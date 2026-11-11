// The simulation engine: one continuous-time run of a population (D024),
// from validated R inputs to ordinary R result tables. Every rate process is
// cell scheduled with one exponential clock (D128); births use cell-accounted
// fertility (D101). All engine state belongs to one Simulation object, so runs
// are independent. Public times are years and rates are per person-year.
#include <Rcpp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "event_queues.h"
#include "fertility_cells.h"
#include "kinship.h"
#include "matching.h"
#include "rate_bands.h"
#include "rate_cells.h"
#include "weight_calibration.h"

using namespace Rcpp;
using namespace hammel;

namespace {

// SHARED DEFINITIONS ------

// Previous-union categories, indexed by their engine codes 0, 1 and 2.
const std::array<const char*, 3> kPreviousLabels = {"none", "separated",
                                                    "widowed"};
// Union-type labels for weight-model member tables, indexed by person value 2.
const std::array<const char*, 3> kUnionTypeLabels = {"unpartnered", "cohab",
                                                     "marriage"};
// Largest integer a double represents exactly: the bound on generated IDs.
constexpr double kMaxExactId = 9007199254740991.0;

// One person's history fields plus working state. `died` is NA while alive;
// `last_birth` enforces birth spacing. Engine code refers to people by their
// zero-based position in Simulation::people, never by public person_id.
struct Person {
  double id, born, died, mother, father, last_birth;
  bool female;
  int group;
};

// Event processes, in their same-time priority order.
enum class EventType {
  death,
  birth,
  dissolution,
  conversion,
  search,
  transition
};

// One union episode. A conversion closes the cohabitation and opens a new
// marriage record linked by `converted_from`.
struct UnionRecord {
  double id, start, end;
  int wife, husband;
  bool wife_died = false;
  bool dissolved = false;
  bool cohab = false, converted = false;
  double converted_from = NA_REAL;
};

// How a process combines group dimensions into strata (D134). A stratum index
// is the sum of category * stride over `dims`; for union processes, `roles`
// says whether each dimension is the wife's (0) or the husband's (1).
struct StrataSpec {
  std::vector<int> dims, roles, strides;
  int count = 1;
};

StrataSpec parse_spec(List x) {
  StrataSpec spec;
  IntegerVector d = x["dims"], r = x["roles"], st = x["strides"];
  spec.dims.assign(d.begin(), d.end());
  spec.roles.assign(r.begin(), r.end());
  spec.strides.assign(st.begin(), st.end());
  spec.count = std::max(1, as<int>(x["count"]));
  return spec;
}

// The active schedules of one mortality/search/conversion stratum.
struct GroupSchedules {
  std::vector<Band> male, female, search, conversion;
  bool duration_conversion = false;
  std::array<std::vector<Band>, 3> previous_search;
  std::array<bool, 3> duration_search = {false, false, false};
};

// Newborn assignment in one group dimension, with that dimension's block of
// the ancestry vector (size zero when ancestry is not tracked). `method` codes
// follow the R constructor: 0 mother, 1 father, 2 random parent, 3 father for
// sons, 4 father for daughters, 5 ancestry draw, 6 ranked, 7 custom, 8 fixed.
struct InheritanceRule {
  int method = 0, count = 1, fixed = -1;
  double passing = 0;
  IntegerVector ranks, ancestry_groups;
  SEXP custom = R_NilValue;
  size_t offset = 0, size = 0;
};

// Read the parts of a rule that a calendar schedule can switch.
void read_rule(List x, InheritanceRule& rule) {
  rule.method = as<int>(x["method"]);
  rule.fixed = as<int>(x["fixed"]);
  rule.passing = as<double>(x["passing"]);
  rule.ranks = as<IntegerVector>(x["ranks"]);
  rule.ancestry_groups = as<IntegerVector>(x["ancestry_groups"]);
  rule.custom = x["custom"];
}

// One outgoing group-transition schedule of an origin stratum: destination,
// optional sex (-1 any, 0 male, 1 female) and union type ("" any). `code` is
// the type as a union-status code: -1 any, 0 unpartnered, 1 cohab, 2
// marriage, 3 partnered (either), 4 anything else (never matches).
struct TransitionSchedule {
  int dest, sex;
  std::string type;
  std::vector<Band> bands;
  int code = -1;
};

int union_type_code(const std::string& type) {
  if (type.empty()) {
    return -1;
  }
  const std::array<const char*, 4> labels = {"unpartnered", "cohab", "marriage",
                                             "partnered"};
  for (int k = 0; k < 4; ++k) {
    if (type == labels[k]) {
      return k;
    }
  }
  return 4;
}

// Transitions in one group dimension. Strata combine the origin category with
// any other dimensions the rate table uses. Household moves (D140) add strata
// [count, 2 * count) that mirror the schedules for dependents, with no risk
// below the dependent age.
struct TransitionProcess {
  int dim = 0, own = 0;  // the dimension, and its position in spec
  StrataSpec spec;
  DataFrame all, active;
  bool duration = false;
  std::vector<std::vector<TransitionSchedule>> schedules;  // per stratum
  std::vector<double> entry;  // entry time into the current category
  std::string weight_key;
  bool household = false;
  double dependent_age = 18;
};

// A weight model (D129) gives each unit a log weight eta; selection within a
// chosen cell is proportional to exp(eta). Formula models hold one factor
// column per term (variable code, level code or -1 for a numeric value,
// ancestry component); custom models call an R function on cell members.
struct WeightSpec {
  bool custom = false, warned_na = false;
  std::string process;
  std::vector<double> coefficients;
  std::vector<std::vector<std::array<int, 3>>> factors;
  SEXP fun = R_NilValue;
};

WeightSpec parse_weight_spec(List x, const std::string& name) {
  WeightSpec spec;
  spec.process = name;
  spec.custom = as<bool>(x["custom"]);
  if (spec.custom) {
    spec.fun = x["fun"];
    return spec;
  }
  spec.coefficients = as<std::vector<double>>(x["coefficients"]);
  List factors = x["factors"];
  for (int k = 0; k < factors.size(); ++k) {
    IntegerMatrix m = factors[k];
    std::vector<std::array<int, 3>> column;
    for (int r = 0; r < m.nrow(); ++r) {
      column.push_back({m(r, 0), m(r, 1), m(r, 2)});
    }
    spec.factors.push_back(column);
  }
  return spec;
}

// Calendar schedule of weight models: intervals and their models. The active
// model is copied into Simulation::weight_specs, so pointers to it stay valid.
struct WeightSchedule {
  std::vector<double> start, end;
  std::vector<WeightSpec> models;
  int active = -1;
};

int weight_interval(const WeightSchedule& w, double now) {
  for (size_t k = 0; k < w.start.size(); ++k) {
    if (w.start[k] <= now && now < w.end[k]) {
      return static_cast<int>(k);
    }
  }
  return -1;
}

// The calendar period of a rate table around now (NA bounds when static).
std::pair<double, double> rate_period(DataFrame rates, double now) {
  if (!rates.containsElementNamed("time_start")) {
    return std::make_pair(NA_REAL, NA_REAL);
  }
  NumericVector lo = rates["time_start"], hi = rates["time_end"];
  double a = R_NegInf, b = R_PosInf;
  for (int r = 0; r < lo.size(); ++r) {
    for (double t : {lo[r], hi[r]}) {
      if (t <= now) {
        a = std::max(a, t);
      } else {
        b = std::min(b, t);
      }
    }
  }
  return std::make_pair(a, b);
}

// Whether a fertility table distinguishes cohabiting from married women.
bool has_cohab_fertility(DataFrame rates) {
  if (!rates.containsElementNamed("union_type")) {
    return false;
  }
  CharacterVector types = rates["union_type"];
  for (auto type : types) {
    if (as<std::string>(type) == "cohab") {
      return true;
    }
  }
  return false;
}

// Whether group weights can change a model's choices. When a model prevents
// all exogamy, every woman sees only her own group's men.
bool applies_weights(const List& model) {
  return !as<bool>(model["prevent_exogamy_all"]);
}

// An Exp(1) cumulative-hazard threshold from R's RNG, which shares the seed
// policy of run_simulation(). Zero is redrawn to keep waiting times positive.
double draw_threshold() {
  double value;
  do {
    value = R::rexp(1.0);
  } while (value <= 0);
  return value;
}

// The strata of each non-fertility process from the active schedules.
struct ProcessStrata {
  std::vector<std::vector<Band>> mortality, search, dissolution, conversion;
  // Per outgoing-schedule slot.
  std::vector<std::vector<std::vector<Band>>> transitions;
  std::vector<std::vector<std::tuple<int, int, std::string>>> transition_layout;
  // (process, slot) of each transition process.
  std::vector<std::pair<int, size_t>> transition_slots;
};

// One event source: a cell-scheduled process, or the fertility cells when
// `cells` is null. `slot` indexes transition processes.
struct EventSource {
  EventType type;
  RateCells* cells;
  int slot;
};

// SIMULATION ------

// All engine state of one run. Construction reads the inputs and places every
// person in their processes; run() processes events until the horizon and
// returns the result tables. Not copyable: processes keep pointers into it.
class Simulation {
 public:
  Simulation(
    DataFrame persons,
    DataFrame mortality_rates,
    DataFrame fertility_rates,
    double start,
    double years,
    double birth_interval,
    DataFrame unions,
    List matching,
    double min_male_age,
    int kinship_generations,
    bool rate_initiation,
    DataFrame union_rates,
    List father_matching,
    bool cancel_birth_if_no_match,
    bool diagnostics,
    NumericMatrix initial_ancestry,
    List inheritance,
    double ffab,
    DataFrame dissolution_rates,
    DataFrame cohab_prob,
    DataFrame conversion_rates,
    List cohab_matching,
    bool counts_cohab,
    bool check_cells,
    Nullable<List> hetfert,
    Nullable<List> transition_input,
    Nullable<List> model_schedule_input,
    Nullable<List> weight_input,
    Nullable<List> group_input
  );
  Simulation(const Simulation&) = delete;
  Simulation& operator=(const Simulation&) = delete;
  List run();

 private:
  // Setup.
  void add_calendar_time(double time);
  void add_calendar_times(const DataFrame& table);
  void slice_rates(double now);
  void read_groups(Nullable<List> group_input, int persons);
  std::map<double, int>
  read_people(DataFrame persons, NumericMatrix initial_ancestry);
  void read_unions(const std::map<double, int>& index);
  void setup_calibration();
  void read_transitions(Nullable<List> transition_input);
  void read_weights(Nullable<List> weight_input);
  void load_probabilities();
  void load_schedules();
  void load_dissolution();
  void load_transitions();
  void load_cells(double now);
  ProcessStrata process_strata();
  void load_processes(double now);
  bool update_process_rates(double now);
  bool update_models(double now);
  bool switch_weights(double now);
  void collect_sources();

  // Group strata and labels.
  int category(int i, int d) const {
    return cats[static_cast<size_t>(i) * dimensions + d];
  }
  int person_stratum(const StrataSpec& spec, int i) const;
  int union_stratum(const StrataSpec& spec, int wife, int husband) const;
  SEXP dim_label(int d, int category) const;

  // Pools, matching and union formation.
  const std::vector<Band>* probability_schedule(int woman) {
    return &probability_schedules
      [person_stratum(probability_spec, woman) * 3 + previous_union[woman]];
  }
  double cohab_probability(int woman, double now);
  void remove_male(int man);
  void add_male(int man);
  void add_counted_woman(int woman, double now);
  void remove_counted_woman(int woman);
  const DrawModel* get_scheduled_draw(int id);
  const std::vector<double>& model_delta(const List& model);
  const LagTable* get_scheduled_lag(int id);
  int find_candidate(
    int woman,
    double now,
    const List& choice_model,
    const DrawModel& draw,
    bool father,
    std::vector<double>* applied,
    bool cohab_role = false
  );
  int form_union(int woman, double now);
  void check_union_capacity() const;
  void
  end_partnership(int wife, int husband, bool cohab, int ending, double now);

  // Weight models.
  WeightSpec* weight_spec(const std::string& name);
  void strata_bin(
    List& bin,
    const StrataSpec& spec,
    int k,
    bool union_level,
    int skip = -1
  );
  double person_value(int var, int component, int i, int dest);
  double union_value(int var, int component, int u);
  void warn_weight_na(WeightSpec& spec, const std::string& what);
  double formula_eta(WeightSpec& spec, int unit, bool is_union, int dest);
  List weight_members(
    const std::vector<int>& units,
    bool is_union,
    int dest,
    const std::string& dest_label = ""
  );
  void custom_eta(
    WeightSpec& spec,
    List members,
    List bin,
    std::vector<double>& eta
  );
  List band_list(
    List bin,
    bool duration,
    double lo,
    double hi,
    DataFrame all,
    double now
  );
  std::function<size_t(int, const std::vector<int>&, double)>
  make_chooser(WeightSpec* spec, RateCells* process, int slot);
  void install_weights(WeightSpec* spec, RateCells* process, int slot);
  void stored_weight(
    WeightSpec* spec,
    RateCells* process,
    int unit,
    bool is_union,
    int dest
  );
  List fertility_bin(int c);
  double search_weight(int woman, int c);
  FertilityCells::State fertility_state(int w);

  // Process membership and household moves.
  void refresh_union(int u, double now);
  void refresh_person(int i, double now);
  void refresh_living(double now);
  int attached_parent(int c) const;
  bool is_dependent(const TransitionProcess& tp, int c, double now) const;
  std::vector<int> household_party(
    const TransitionProcess& tp,
    std::vector<int> party,
    double now
  );
  void move_party(
    TransitionProcess& tp,
    const std::vector<int>& party,
    int destination,
    double initiator,
    double now
  );
  void rejoin(int c, double now);

  // Events.
  double source_clock(const EventSource& x) const {
    return x.cells ? x.cells->next_clock() : cells->next_clock();
  }
  double source_boundary(const EventSource& x) {
    return x.cells ? x.cells->next_boundary() : cells->next_boundary();
  }
  bool apply_calendar_boundary();
  void on_transition(int i, int slot, double now);
  void on_death(int i, double now);
  void on_dissolution(int i, int unit, double now);
  void on_conversion(int i, int unit, double now);
  void on_search(int i, double now);
  void on_birth(int i, double now);
  void append_ancestry(int child, int mother, int father);
  int newborn_category(
    int d,
    int child,
    int mother,
    int father,
    bool female,
    bool partnered
  );
  List results();

  // Inputs that persist through the run. Rate tables hold the active
  // calendar slice; all_* tables hold every period.
  DataFrame mortality_rates, fertility_rates, union_rates, dissolution_rates,
    cohab_prob, conversion_rates, unions;
  DataFrame all_mortality, all_fertility, all_search, all_dissolution,
    all_conversion, all_probability;
  List matching, father_matching, cohab_matching;
  double start, years, end, birth_interval, min_male_age, ffab;
  bool rate_initiation, cancel_birth_if_no_match, diagnostics, counts_cohab,
    check_cells;
  std::set<double> calendar_boundaries;

  // Calendar schedules of matching and inheritance models (D123).
  bool scheduled_models = false;
  List matching_bank, inheritance_bank;
  NumericVector model_times;
  IntegerVector general_route, father_route, cohab_route;
  std::vector<IntegerVector> inheritance_routes;
  int model_epoch = 0, general_id = 0, father_id = 0, cohab_id = 0;

  // Fertility heterogeneity: persistent weights for women (NA for men).
  bool heterogeneous = false;
  double fertility_alpha = 0, fertility_variance = 0;
  std::vector<double> fertility_weights;
  bool cohab_enabled = false, parity_specific = false;

  // Group dimensions (D134): one category per person per dimension, stored
  // row-major. The matching category combines the matching dimensions and
  // carries pools, pair effects, calibration and transitions.
  int dimensions = 0;
  std::vector<int> cats;
  std::vector<std::vector<std::string>> dim_labels;
  std::vector<std::string> dim_names;
  StrataSpec matching_spec, mortality_spec, fertility_spec, search_spec,
    probability_spec, conversion_spec, dissolution_spec, report_spec;
  int group_count = 1;
  std::vector<char> matching_dim;  // whether dimension d is a matching one

  // People and their links. Ancestry shares are row-major; each tracked
  // dimension's block ends with its unknown share.
  std::vector<Person> people;
  double next_id = 0, min_birth = 0;
  std::vector<int> mother_of, father_of, parity;
  std::vector<std::vector<int>> kids;
  size_t ancestry_components = 0;
  std::vector<double> ancestry;
  std::vector<char> unknown_component;
  std::vector<InheritanceRule> inheritance_rules;
  std::unique_ptr<Kinship> kinship;

  // Unions and union status.
  std::vector<UnionRecord> union_history;
  std::vector<std::vector<int>> person_unions;
  std::vector<int> current_partner;
  bool previous_specific = false;
  std::vector<int> previous_union, legal_previous;
  std::vector<double> previous_union_end;
  double next_union_id = 0;
  // Log group weights applied at each union-forming search (empty = unit).
  std::vector<std::vector<double>> diagnostic_log_weights;

  // Exact partner draws (D072): living, unpartnered men by group and bin.
  MalePool pool;
  DrawModel union_draw, father_draw_storage, cohab_draw_storage;
  const DrawModel* union_draw_current = &union_draw;
  const DrawModel* father_draw = &union_draw;
  const DrawModel* cohab_draw = &union_draw;
  std::map<int, DrawModel> scheduled_draws;
  DrawWorkspace draw_workspace;
  double draw_approximate = 0;

  // weighted_choice calibration (D066), present only for grouped runs with
  // weighted matching. Counted women are living and unpartnered.
  std::unique_ptr<WeightCalibrator> calibrator;
  std::vector<char> counted_woman;
  std::unique_ptr<ProbabilityTracker> probability_tracker;
  SearchCalibration search_calibration;
  bool search_calibrated = false;
  LagTable union_table, father_table_storage, cohab_table;
  int union_table_id = 0, father_table_id = 0, cohab_table_id = 0;
  bool union_weighted = false, father_weighted = false, cohab_weighted = false;
  std::map<int, LagTable> scheduled_lags;
  std::vector<double> union_warm, father_warm, cohab_warm;
  std::map<SEXP, std::vector<double>> deltas;  // target log shifts by model
  double weighted_infeasible = 0;

  // Active schedules.
  std::vector<std::vector<Band>> probability_schedules;
  std::vector<GroupSchedules> schedules;
  int schedule_count = 0;
  bool status_fertility = false, typed_fertility = false;
  bool fertility_grouped = false;  // fertility cells stratify by group
  bool has_dissolution = false;
  std::vector<std::vector<Band>> dissolution_schedules,
    cohab_dissolution_schedules;
  // Fertility union-type labels: overall, unpartnered/partnered, or three.
  std::vector<std::string> fertility_types;

  // Rate processes. Search, dissolution and conversion are null when not
  // configured; there is one transition process per outgoing-schedule slot.
  std::unique_ptr<FertilityCells> cells;
  std::unique_ptr<RateCells> mortality, search, dissolution, conversion;
  std::vector<std::unique_ptr<RateCells>> transitions;
  std::vector<std::pair<int, size_t>> transition_slots;
  std::vector<std::vector<std::tuple<int, int, std::string>>> transition_layout;
  std::vector<EventSource> sources;

  // Group transitions and household moves, with their recorded history.
  std::vector<TransitionProcess> tprocs;
  bool any_household = false;
  std::vector<double> tx_person, tx_time, tx_move, tx_initiator;
  std::vector<int> tx_dim, tx_origin, tx_dest;
  double move_count = 0;

  // Weight models.
  std::map<std::string, WeightSpec> weight_specs;
  std::map<std::string, WeightSchedule> weight_schedules;
  std::vector<std::string> weight_ancestry_labels;
  WeightSpec* fertility_weights_spec = nullptr;
  WeightSpec* search_weights_spec = nullptr;
  WeightSpec* mortality_weights_spec = nullptr;
  WeightSpec* dissolution_weights_spec = nullptr;
  WeightSpec* conversion_weights_spec = nullptr;
  std::vector<WeightSpec*> transition_weight_specs;
  bool fertility_refresh_weights = false;
  double weight_now = 0;  // event time, for custom-function rate periods

  // Unrealized births: (period, group, union type, parity start/end, age
  // start/end) -> count.
  std::map<
    std::tuple<long long, int, std::string, double, double, double, double>,
    double>
    unrealized_counts;
  double birth_count = 0, death_count = 0, dissolution_count = 0,
         conversion_count = 0;
  unsigned long long processed = 0;
  bool cells_consistent = true;
};

// SETUP ------

// Constructor parameters that persist are copied into members of the same
// name; the body reads the inputs in a fixed order, which also fixes the
// order of any random draws and warnings during setup.
Simulation::Simulation(
  DataFrame persons,
  DataFrame mortality_rates,
  DataFrame fertility_rates,
  double start,
  double years,
  double birth_interval,
  DataFrame unions,
  List matching,
  double min_male_age,
  int kinship_generations,
  bool rate_initiation,
  DataFrame union_rates,
  List father_matching,
  bool cancel_birth_if_no_match,
  bool diagnostics,
  NumericMatrix initial_ancestry,
  List inheritance,
  double ffab,
  DataFrame dissolution_rates,
  DataFrame cohab_prob,
  DataFrame conversion_rates,
  List cohab_matching,
  bool counts_cohab,
  bool check_cells,
  Nullable<List> hetfert,
  Nullable<List> transition_input,
  Nullable<List> model_schedule_input,
  Nullable<List> weight_input,
  Nullable<List> group_input
)
    : unions(unions),
      all_mortality(mortality_rates),
      all_fertility(fertility_rates),
      all_search(union_rates),
      all_dissolution(dissolution_rates),
      all_conversion(conversion_rates),
      all_probability(cohab_prob),
      matching(matching),
      father_matching(father_matching),
      cohab_matching(cohab_matching),
      start(start),
      years(years),
      end(start + years),
      birth_interval(birth_interval),
      min_male_age(min_male_age),
      ffab(ffab),
      rate_initiation(rate_initiation),
      cancel_birth_if_no_match(cancel_birth_if_no_match),
      diagnostics(diagnostics),
      counts_cohab(counts_cohab),
      check_cells(check_cells) {
  // R performs input validation and seed handling.
  scheduled_models = model_schedule_input.isNotNull();
  if (scheduled_models) {
    List config(model_schedule_input.get());
    matching_bank = config["matching_models"];
    inheritance_bank = config["inheritance_models"];
    DataFrame epochs = as<DataFrame>(config["epochs"]);
    model_times = epochs["time_start"];
    general_route = epochs["matching"];
    father_route = epochs["father_matching"];
    cohab_route = epochs["cohab_matching"];
    List routes = config["inheritance_routes"];
    for (int d = 0; d < routes.size(); ++d) {
      inheritance_routes.push_back(as<IntegerVector>(routes[d]));
    }
    general_id = general_route[0];
    father_id = father_route[0];
    cohab_id = cohab_route[0];
  }
  for (const DataFrame& table :
       {all_mortality, all_fertility, all_search, all_dissolution,
        all_conversion, all_probability}) {
    add_calendar_times(table);
  }
  if (scheduled_models) {
    for (double time : model_times) {
      add_calendar_time(time);
    }
  }
  slice_rates(start);
  this->cohab_prob = calendar_slice(all_probability, start);

  const List heterogeneity =
    hetfert.isNotNull() ? List(hetfert.get()) : List::create();
  heterogeneous = heterogeneity.size() > 0;
  const double fertility_cv =
    heterogeneous ? as<double>(heterogeneity["cv"]) : 0;
  fertility_alpha = heterogeneous ? as<double>(heterogeneity["alpha"]) : 0;
  fertility_variance =
    fertility_cv <= 1 ? std::log1p(fertility_cv * fertility_cv)
                      : 2 * std::log(fertility_cv) +
                          std::log1p((1 / fertility_cv) * (1 / fertility_cv));
  if (heterogeneous) {
    NumericVector supplied = persons["hetfert"];
    fertility_weights.assign(supplied.begin(), supplied.end());
  }
  cohab_enabled = this->cohab_prob.nrows() > 0;
  parity_specific = this->fertility_rates.containsElementNamed("parity_start");
  if (!R_finite(ffab) || ffab < 0 || ffab > 1) {
    stop("ffab must be a finite probability in [0, 1].");
  }
  read_groups(group_input, persons.nrows());
  const std::map<double, int> index = read_people(persons, initial_ancestry);

  // Draw models for the union, father and cohabitation matching models. An
  // empty population never draws, so custom-formula grids are skipped then.
  if (!people.empty() && !scheduled_models) {
    union_draw = make_draw_model(matching, group_count, end - min_birth);
  }
  if (!scheduled_models && !people.empty() &&
      SEXP(father_matching) != SEXP(matching)) {
    father_draw_storage =
      make_draw_model(father_matching, group_count, end - min_birth);
    father_draw = &father_draw_storage;
  }
  if (!scheduled_models && cohab_enabled && !people.empty() &&
      SEXP(cohab_matching) != SEXP(matching)) {
    cohab_draw_storage =
      make_draw_model(cohab_matching, group_count, end - min_birth);
    cohab_draw = &cohab_draw_storage;
  }
  if (scheduled_models) {
    union_draw_current = get_scheduled_draw(general_id);
    father_draw = get_scheduled_draw(father_id);
    cohab_draw = get_scheduled_draw(cohab_id);
  }
  std::vector<std::array<int, 2>> parent_links(people.size());
  for (size_t i = 0; i < people.size(); ++i) {
    parent_links[i] = {mother_of[i], father_of[i]};
  }
  kinship.reset(new Kinship(std::move(parent_links), kinship_generations));

  person_unions.resize(people.size());
  current_partner.assign(people.size(), -1);
  previous_specific = this->union_rates.containsElementNamed("previous_union");
  previous_union.assign(people.size(), 0);
  legal_previous.assign(people.size(), 0);
  previous_union_end.assign(people.size(), NA_REAL);
  if (previous_specific || cohab_enabled) {
    previous_union = as<std::vector<int>>(persons[".hammel_previous_union"]);
    previous_union_end =
      as<std::vector<double>>(persons[".hammel_previous_union_end"]);
    legal_previous = as<std::vector<int>>(persons[".hammel_legal_previous"]);
  }
  probability_schedules.resize(probability_spec.count * 3);
  load_probabilities();
  read_unions(index);

  union_weighted = applies_weights(matching);
  father_weighted = applies_weights(father_matching);
  cohab_weighted = cohab_enabled && applies_weights(cohab_matching);
  for (int d = 0; d < dimensions; ++d) {
    List x = inheritance[d];
    read_rule(x, inheritance_rules[d]);
    inheritance_rules[d].count = as<int>(x["group_count"]);
    inheritance_rules[d].offset =
      static_cast<size_t>(as<int>(x["ancestry_offset"]));
    inheritance_rules[d].size =
      static_cast<size_t>(as<int>(x["ancestry_size"]));
  }
  unknown_component.assign(ancestry_components, 0);
  for (const auto& rule : inheritance_rules) {
    if (rule.size) {
      unknown_component[rule.offset + rule.size - 1] = 1;
    }
  }
  schedule_count =
    std::max({mortality_spec.count, search_spec.count, conversion_spec.count});
  schedules.resize(schedule_count);
  status_fertility = this->fertility_rates.containsElementNamed("union_type");
  load_schedules();
  setup_calibration();
  has_dissolution = this->dissolution_rates.nrows() > 0;
  load_dissolution();
  read_transitions(transition_input);
  load_transitions();
  read_weights(weight_input);
  weight_now = start;
  load_cells(start);
  load_processes(start);
  for (int i = 0; i < static_cast<int>(people.size()); ++i) {
    if (i % 1024 == 0) {
      checkUserInterrupt();
    }
    refresh_person(i, start);
  }
  collect_sources();
}

void Simulation::add_calendar_time(double time) {
  if (start < time && time < end) {
    calendar_boundaries.insert(time);
  }
}

// Record a rate table's calendar band edges inside the run as boundaries.
void Simulation::add_calendar_times(const DataFrame& table) {
  if (!table.containsElementNamed("time_start")) {
    return;
  }
  for (const char* column : {"time_start", "time_end"}) {
    NumericVector times = table[column];
    for (double time : times) {
      add_calendar_time(time);
    }
  }
}

// Select the calendar slices of the demographic rate tables active at now.
// Cohabitation probabilities and transitions are sliced separately.
void Simulation::slice_rates(double now) {
  mortality_rates = calendar_slice(all_mortality, now);
  fertility_rates = calendar_slice(all_fertility, now);
  union_rates = calendar_slice(all_search, now);
  dissolution_rates = calendar_slice(all_dissolution, now);
  conversion_rates = calendar_slice(all_conversion, now);
}

// Group dimensions, their labels, each person's categories, and each
// process's strata.
void Simulation::read_groups(Nullable<List> group_input, int persons) {
  const List groups =
    group_input.isNotNull() ? List(group_input.get()) : List::create();
  auto read_spec = [&](const char* name) {
    if (!groups.size() || !groups.containsElementNamed(name)) {
      return StrataSpec();
    }
    return parse_spec(as<List>(groups[name]));
  };
  dimensions = groups.size() ? as<int>(groups["dimensions"]) : 0;
  cats.resize(static_cast<size_t>(persons) * dimensions);
  if (dimensions) {
    IntegerMatrix input = groups["categories"];
    for (int i = 0; i < persons; ++i) {
      for (int d = 0; d < dimensions; ++d) {
        cats[static_cast<size_t>(i) * dimensions + d] = input(i, d);
      }
    }
  }
  if (groups.size()) {
    List labels = groups["labels"];
    for (int d = 0; d < labels.size(); ++d) {
      dim_labels.push_back(as<std::vector<std::string>>(labels[d]));
    }
    dim_names = as<std::vector<std::string>>(groups["names"]);
  }
  matching_spec = read_spec("matching");
  group_count = matching_spec.count;
  matching_dim.assign(dimensions, 0);
  for (int d : matching_spec.dims) {
    matching_dim[d] = 1;
  }
  mortality_spec = read_spec("mortality");
  fertility_spec = read_spec("fertility");
  search_spec = read_spec("search");
  probability_spec = read_spec("probability");
  conversion_spec = read_spec("conversion");
  dissolution_spec = read_spec("dissolution");
  report_spec = read_spec("report");
  inheritance_rules.resize(dimensions);
}

// People, parent links, ancestry, parity and the male pool from the input.
// Returns the map from public IDs to positions for reading union histories.
std::map<double, int>
Simulation::read_people(DataFrame persons, NumericMatrix initial_ancestry) {
  NumericVector ids = persons["person_id"], births = persons["birth_time"];
  NumericVector deaths = persons["death_time"], moms = persons["mother_id"];
  NumericVector dads = persons["father_id"];
  CharacterVector sexes = persons["sex"];
  // Row-major ancestry shares; zero components means ancestry is not tracked.
  ancestry_components = initial_ancestry.ncol();
  if (initial_ancestry.nrow() != ids.size()) {
    stop("Ancestry rows must match persons.");
  }
  ancestry.resize(static_cast<size_t>(ids.size()) * ancestry_components);
  for (int i = 0; i < ids.size(); ++i) {
    for (size_t g = 0; g < ancestry_components; ++g) {
      ancestry[static_cast<size_t>(i) * ancestry_components + g] =
        initial_ancestry(i, g);
    }
  }
  std::map<double, int> index;
  min_birth = start;
  for (int i = 0; i < ids.size(); ++i) {
    const bool female = as<std::string>(sexes[i]) == "f";
    people.push_back(
      {ids[i], births[i], deaths[i], moms[i], dads[i],
       -std::numeric_limits<double>::infinity(), female,
       person_stratum(matching_spec, i)}
    );
    index[ids[i]] = i;
    next_id = std::max(next_id, ids[i]);
    min_birth = std::min(min_birth, births[i]);
  }
  // Parent positions (-1 unknown) and each person's children, for kinship
  // and household moves (D140). Validation guarantees recorded parents exist.
  mother_of.assign(people.size(), -1);
  father_of.assign(people.size(), -1);
  kids.resize(people.size());
  for (int i = 0; i < static_cast<int>(people.size()); ++i) {
    if (!NumericVector::is_na(people[i].mother)) {
      mother_of[i] = index.at(people[i].mother);
      kids[mother_of[i]].push_back(i);
    }
    if (!NumericVector::is_na(people[i].father)) {
      father_of[i] = index.at(people[i].father);
      kids[father_of[i]].push_back(i);
    }
  }
  // Parity selects fertility cells when schedules are parity-specific; it is
  // always tracked because weight models may use it. On restart, the latest
  // recorded birth, including to children who died, restores spacing.
  parity.assign(people.size(), 0);
  for (int i = 0; i < static_cast<int>(people.size()); ++i) {
    const int m = mother_of[i];
    if (m < 0) {
      continue;
    }
    people[m].last_birth = std::max(people[m].last_birth, people[i].born);
    if (people[m].female) {
      ++parity[m];
    }
  }
  pool = MalePool(min_birth, end, group_count, kWeightBinWidth);
  for (int i = 0; i < static_cast<int>(people.size()); ++i) {
    if (!people[i].female && NumericVector::is_na(people[i].died)) {
      pool.add(i, people[i].group, people[i].born);
    }
  }
  return index;
}

// Input union histories. Ongoing unions set current partners and remove the
// husband from the pool.
void Simulation::read_unions(const std::map<double, int>& index) {
  NumericVector uids = unions["union_id"], wives = unions["wife_id"];
  NumericVector husbands = unions["husband_id"], starts = unions["start_time"];
  NumericVector ends = unions["end_time"];
  const bool typed = unions.containsElementNamed("union_type");
  const bool linked = unions.containsElementNamed("previous_cohab_id");
  CharacterVector types =
    typed ? CharacterVector(unions["union_type"]) : CharacterVector();
  NumericVector links =
    linked ? NumericVector(unions["previous_cohab_id"]) : NumericVector();
  for (int u = 0; u < uids.size(); ++u) {
    const int woman = index.at(wives[u]), man = index.at(husbands[u]);
    union_history.push_back({uids[u], starts[u], ends[u], woman, man});
    if (typed) {
      union_history.back().cohab = as<std::string>(types[u]) == "cohab";
    }
    if (linked) {
      union_history.back().converted_from = links[u];
    }
    next_union_id = std::max(next_union_id, uids[u]);
    if (!NumericVector::is_na(ends[u])) {
      continue;
    }
    if (!NumericVector::is_na(people[woman].died) ||
        !NumericVector::is_na(people[man].died) || !people[woman].female ||
        people[man].female) {
      stop(
        "Unions require living female/male partners in active input unions."
      );
    }
    if (current_partner[woman] >= 0) {
      stop("Union input cannot have multiple ongoing unions for a woman.");
    }
    if (current_partner[man] >= 0) {
      stop("Union input cannot have multiple ongoing unions for a man.");
    }
    current_partner[woman] = man;
    current_partner[man] = woman;
    remove_male(man);
    person_unions[woman].push_back(u);
    person_unions[man].push_back(u);
  }
}

// weighted_choice calibration (D066) needs schedules and group indices, so it
// starts after them. One group means unit weights: no calibration.
void Simulation::setup_calibration() {
  bool any_weighted = union_weighted || father_weighted || cohab_weighted;
  if (scheduled_models) {
    for (auto model : matching_bank) {
      any_weighted = any_weighted || applies_weights(as<List>(model));
    }
  }
  if (!any_weighted || group_count <= 1) {
    return;
  }
  calibrator.reset(new WeightCalibrator(
    min_birth, end, group_count, kWeightBinWidth, min_male_age
  ));
  if (scheduled_models) {
    union_table_id = calibrator->add_table(get_scheduled_lag(general_id));
    father_table_id = calibrator->add_table(get_scheduled_lag(father_id));
    cohab_table_id = calibrator->add_table(get_scheduled_lag(cohab_id));
  } else {
    // A NULL father model shares the general model and calibration table.
    if (union_weighted) {
      union_table = build_lag_table(
        matching, group_count, calibrator->bins(), kWeightBinWidth
      );
      union_table_id = father_table_id = calibrator->add_table(&union_table);
    }
    if (father_weighted && SEXP(father_matching) != SEXP(matching)) {
      father_table_storage = build_lag_table(
        father_matching, group_count, calibrator->bins(), kWeightBinWidth
      );
      father_table_id = calibrator->add_table(&father_table_storage);
    }
    if (cohab_weighted) {
      if (SEXP(cohab_matching) == SEXP(matching)) {
        cohab_table_id = union_table_id;
      } else {
        cohab_table = build_lag_table(
          cohab_matching, group_count, calibrator->bins(), kWeightBinWidth
        );
        cohab_table_id = calibrator->add_table(&cohab_table);
      }
    }
  }
  counted_woman.assign(people.size(), 0);
  // Rate-mode searches take their intensities from the search cells
  // (SearchCalibration) and cohabitation probabilities from the tracker.
  search_calibrated = rate_initiation;
  if (rate_initiation && cohab_enabled) {
    probability_tracker.reset(
      new ProbabilityTracker([this](int woman, double p) {
        search_calibration.set_probability(woman, p);
      })
    );
  }
  for (int man = 0; man < static_cast<int>(people.size()); ++man) {
    if (pool.contains(man)) {
      calibrator->add_man(people[man].group, people[man].born, 1);
    }
  }
  for (int woman = 0; woman < static_cast<int>(people.size()); ++woman) {
    add_counted_woman(woman, start);
  }
}

// One transition process per dimension with transitions, in declaration
// order.
void Simulation::read_transitions(Nullable<List> transition_input) {
  if (transition_input.isNotNull()) {
    List input(transition_input.get());
    for (int t = 0; t < input.size(); ++t) {
      List x = input[t];
      TransitionProcess tp;
      tp.dim = as<int>(x["dimension"]);
      tp.spec = parse_spec(as<List>(x["spec"]));
      for (size_t j = 0; j < tp.spec.dims.size(); ++j) {
        if (tp.spec.dims[j] == tp.dim) {
          tp.own = static_cast<int>(j);
        }
      }
      tp.all = as<DataFrame>(x["rates"]);
      tp.active = calendar_slice(tp.all, start);
      tp.entry = as<std::vector<double>>(x["entry"]);
      tp.weight_key = "transition_" + std::to_string(t + 1);
      if (x.containsElementNamed("household")) {
        tp.household = as<bool>(x["household"]);
        tp.dependent_age = as<double>(x["dependent_age"]);
      }
      tp.schedules.resize(
        static_cast<size_t>(tp.spec.count) * (tp.household ? 2 : 1)
      );
      add_calendar_times(tp.all);
      tprocs.push_back(tp);
    }
  }
  for (const auto& tp : tprocs) {
    any_household = any_household || tp.household;
  }
}

// Weight models by process, with their calendar schedules. Formula models
// whose coefficients are all zero are unweighted.
void Simulation::read_weights(Nullable<List> weight_input) {
  if (weight_input.isNotNull()) {
    List input(weight_input.get());
    weight_ancestry_labels =
      as<std::vector<std::string>>(input["ancestry_labels"]);
    std::vector<std::pair<std::string, std::string>> keys;  // input, process
    for (const char* name :
         {"mortality", "fertility", "dissolution", "conversion", "search"}) {
      keys.push_back({name, name});
    }
    for (const auto& tp : tprocs) {
      keys.push_back({tp.weight_key, "transition"});
    }
    for (const auto& key : keys) {
      const std::string& name = key.first;
      if (!input.containsElementNamed(name.c_str())) {
        continue;
      }
      List x = input[name];
      if (x.containsElementNamed("schedule")) {
        WeightSchedule w;
        w.start = as<std::vector<double>>(x["time_start"]);
        w.end = as<std::vector<double>>(x["time_end"]);
        List models = x["models"];
        for (int k = 0; k < models.size(); ++k) {
          w.models.push_back(parse_weight_spec(models[k], key.second));
        }
        w.active = weight_interval(w, start);
        if (w.active < 0) {
          stop("No " + key.second + " weight model covers the run start.");
        }
        // Scheduled models stay installed even when an interval is
        // unweighted.
        weight_specs[name] = w.models[w.active];
        for (size_t k = 0; k < w.start.size(); ++k) {
          add_calendar_time(w.start[k]);
          add_calendar_time(w.end[k]);
        }
        weight_schedules[name] = w;
        continue;
      }
      WeightSpec spec = parse_weight_spec(x, key.second);
      if (!spec.custom && spec.coefficients.empty()) {
        continue;
      }
      weight_specs[name] = spec;
    }
  }
  // weight_specs is complete: entries are replaced in place, never added.
  fertility_weights_spec = weight_spec("fertility");
  // Formula weights change only at events, so they are stored at refresh
  // and members are chosen by rejection sampling. Custom functions are
  // evaluated when a cell is chosen, except fertility under matching
  // calibration, which stores them at refresh too.
  fertility_refresh_weights =
    fertility_weights_spec && (calibrator || !fertility_weights_spec->custom);
  search_weights_spec = weight_spec("search");
  mortality_weights_spec = weight_spec("mortality");
  dissolution_weights_spec = weight_spec("dissolution");
  conversion_weights_spec = weight_spec("conversion");
  for (const auto& tp : tprocs) {
    transition_weight_specs.push_back(weight_spec(tp.weight_key));
  }
  if (search_calibrated && search_weights_spec && search_weights_spec->custom) {
    Rcpp::warning(
      "With weighted matching calibration, the custom search weight function "
      "is evaluated at every refresh of every searching woman, which can be "
      "slow."
    );
  }
  if (fertility_refresh_weights && fertility_weights_spec->custom &&
      calibrator) {
    Rcpp::warning(
      "With weighted matching calibration, the custom fertility weight "
      "function is evaluated at every refresh of every woman, which can be "
      "slow."
    );
  }
}

void Simulation::load_probabilities() {
  if (!cohab_enabled) {
    return;
  }
  const bool previous_prob = cohab_prob.containsElementNamed("previous_union");
  for (int g = 0; g < probability_spec.count; ++g) {
    for (int p = 0; p < 3; ++p) {
      probability_schedules[g * 3 + p] = read_bands(
        cohab_prob, "", "", g, previous_prob ? kPreviousLabels[p] : ""
      );
    }
  }
}

void Simulation::load_schedules() {
  for (int group = 0; group < schedule_count; ++group) {
    auto& rates = schedules[group];
    rates.male = read_bands(mortality_rates, "m", "", group);
    rates.female = read_bands(mortality_rates, "f", "", group);
    if (conversion_rates.nrows()) {
      NumericVector ages = conversion_rates["age_start"];
      rates.duration_conversion = NumericVector::is_na(ages[0]);
      rates.conversion = read_bands(
        conversion_rates, "", "", group, "", rates.duration_conversion
      );
    }
    if (!rate_initiation) {
      continue;
    }
    if (!previous_specific) {
      rates.search = read_bands(union_rates, "", "", group);
      continue;
    }
    for (int p = 0; p < 3; ++p) {
      rates.duration_search[p] =
        search_uses_duration(union_rates, kPreviousLabels[p]);
      rates.previous_search[p] = read_bands(
        union_rates, "", "", group, kPreviousLabels[p], rates.duration_search[p]
      );
    }
  }
}

void Simulation::load_dissolution() {
  dissolution_schedules.clear();
  cohab_dissolution_schedules.clear();
  if (!has_dissolution) {
    return;
  }
  const bool typed = dissolution_rates.containsElementNamed("union_type");
  for (int pair = 0; pair < dissolution_spec.count; ++pair) {
    dissolution_schedules.push_back(
      read_bands(dissolution_rates, "", typed ? "marriage" : "", pair)
    );
    cohab_dissolution_schedules.push_back(
      read_bands(dissolution_rates, "", typed ? "cohab" : "", pair)
    );
  }
}

// Each stratum lists its outgoing schedules (destination, optional sex and
// union type), built from the active transition table.
void Simulation::load_transitions() {
  for (auto& tp : tprocs) {
    for (auto& stratum : tp.schedules) {
      stratum.clear();
    }
    DataFrame input = tp.active;
    tp.duration = input.containsElementNamed("dur_start");
    IntegerVector strata = input[".hammel_group"],
                  destinations = input[".dest"];
    NumericVector lo = input[tp.duration ? "dur_start" : "age_start"];
    NumericVector hi = input[tp.duration ? "dur_end" : "age_end"],
                  values = input["rate"];
    const bool sexed = input.containsElementNamed("sex"),
               typed = input.containsElementNamed("union_type");
    CharacterVector sexes =
      sexed ? as<CharacterVector>(input["sex"]) : CharacterVector();
    CharacterVector types =
      typed ? as<CharacterVector>(input["union_type"]) : CharacterVector();
    for (int row = 0; row < input.nrows(); ++row) {
      auto& outgoing = tp.schedules[strata[row]];
      const int sex =
        !sexed ? -1 : (as<std::string>(sexes[row]) == "f" ? 1 : 0);
      const std::string type = typed ? as<std::string>(types[row]) : "";
      auto found = std::find_if(
        outgoing.begin(), outgoing.end(), [&](const TransitionSchedule& x) {
          return x.dest == destinations[row] && x.sex == sex && x.type == type;
        }
      );
      if (found == outgoing.end()) {
        outgoing.push_back({destinations[row], sex, type, {}});
        outgoing.back().code = union_type_code(type);
        found = std::prev(outgoing.end());
      }
      found->bands.push_back({lo[row], hi[row], values[row]});
    }
    for (auto& stratum : tp.schedules) {
      for (auto& schedule : stratum) {
        std::sort(
          schedule.bands.begin(), schedule.bands.end(),
          [](const Band& a, const Band& b) { return a.start < b.start; }
        );
      }
    }
    if (!tp.household) {
      continue;
    }
    // Dependents' mirror strata: the same schedules from the dependent age.
    for (int st = 0; st < tp.spec.count; ++st) {
      auto& mirror = tp.schedules[st + tp.spec.count];
      mirror = tp.schedules[st];
      for (auto& schedule : mirror) {
        std::vector<Band> kept;
        for (Band band : schedule.bands) {
          if (band.end <= tp.dependent_age) {
            continue;
          }
          band.start = std::max(band.start, tp.dependent_age);
          kept.push_back(band);
        }
        schedule.bands = kept;
      }
    }
  }
}

// (Re)build the fertility cells from the active fertility table.
void Simulation::load_cells(double now) {
  typed_fertility = has_cohab_fertility(fertility_rates);
  fertility_types =
    !status_fertility
      ? std::vector<std::string>{""}
      : (typed_fertility
           ? std::vector<std::string>{"unpartnered", "marriage", "cohab"}
           : std::vector<std::string>{"unpartnered", "partnered"});
  double oldest_birth = end;
  for (const Person& person : people) {
    oldest_birth = std::min(oldest_birth, person.born);
  }
  fertility_grouped = fertility_rates.containsElementNamed(".hammel_group");
  const int fertility_groups = fertility_spec.count;
  cells.reset(new FertilityCells(
    fertility_strata(fertility_rates, fertility_groups, fertility_types),
    fertility_groups, static_cast<int>(fertility_types.size()), parity_specific,
    end - oldest_birth, birth_interval,
    heterogeneous || fertility_refresh_weights,
    [this](int w) { return fertility_state(w); }, draw_threshold,
    calibrator.get(), now, end
  ));
  // Unrealized births are attributed by members' categories in every
  // dimension.
  cells->set_pooled_groups(report_spec.count > 1);
  cells->set_unrealized_hook([this](int c, int group) {
    // Record one unrealized birth in its one-year period from the run start.
    long long k = static_cast<long long>(std::floor(weight_now - start));
    if (start + static_cast<double>(k + 1) <= weight_now) {
      ++k;  // roundoff at a period boundary
    }
    const auto info = cells->cell_info(c);
    ++unrealized_counts[std::make_tuple(
      k, group, fertility_types[info.type], info.parity_start, info.parity_end,
      info.age_start, info.age_end
    )];
  });
  WeightSpec* spec = fertility_weights_spec;
  if (!spec) {
    return;
  }
  if (fertility_refresh_weights) {
    cells->set_refresh_eta([this, spec](int woman, int c) {
      if (!spec->custom) {
        return formula_eta(*spec, woman, false, -1);
      }
      std::vector<double> eta;
      custom_eta(
        *spec, weight_members({woman}, false, -1), fertility_bin(c), eta
      );
      return eta[0];
    });
    return;
  }
  cells->set_select_eta(
    [this,
     spec](int c, const std::vector<int>& women, std::vector<double>& eta) {
      if (spec->custom) {
        custom_eta(
          *spec, weight_members(women, false, -1), fertility_bin(c), eta
        );
        return;
      }
      eta.resize(women.size());
      for (size_t k = 0; k < women.size(); ++k) {
        eta[k] = formula_eta(*spec, women[k], false, -1);
      }
    }
  );
}

// Strata of the non-fertility processes. Mortality: group x sex. Search:
// group (x previous-union category). Dissolution: group pair x union type.
// Conversion: group. Transition slot k has one stratum per origin stratum
// (its k-th outgoing schedule).
ProcessStrata Simulation::process_strata() {
  ProcessStrata x;
  for (int g = 0; g < mortality_spec.count; ++g) {
    x.mortality.push_back(schedules[g].male);
    x.mortality.push_back(schedules[g].female);
  }
  if (rate_initiation) {
    for (int g = 0; g < search_spec.count; ++g) {
      if (!previous_specific) {
        x.search.push_back(schedules[g].search);
        continue;
      }
      for (int k = 0; k < 3; ++k) {
        x.search.push_back(schedules[g].previous_search[k]);
      }
    }
  }
  if (has_dissolution) {
    for (size_t pair = 0; pair < dissolution_schedules.size(); ++pair) {
      x.dissolution.push_back(dissolution_schedules[pair]);
      x.dissolution.push_back(cohab_dissolution_schedules[pair]);
    }
  }
  if (cohab_enabled) {
    for (int g = 0; g < conversion_spec.count; ++g) {
      x.conversion.push_back(schedules[g].conversion);
    }
  }
  for (size_t p = 0; p < tprocs.size(); ++p) {
    const TransitionProcess& tp = tprocs[p];
    size_t slots = 0;
    for (const auto& stratum : tp.schedules) {
      slots = std::max(slots, stratum.size());
    }
    for (size_t k = 0; k < slots; ++k) {
      std::vector<std::vector<Band>> strata;
      for (const auto& stratum : tp.schedules) {
        strata.push_back(
          k < stratum.size() ? stratum[k].bands : std::vector<Band>()
        );
      }
      x.transitions.push_back(strata);
      x.transition_slots.push_back({static_cast<int>(p), k});
    }
    for (const auto& stratum : tp.schedules) {
      std::vector<std::tuple<int, int, std::string>> layout;
      for (const auto& schedule : stratum) {
        layout.emplace_back(schedule.dest, schedule.sex, schedule.type);
      }
      x.transition_layout.push_back(layout);
    }
  }
  return x;
}

// (Re)build the non-fertility rate processes from the active schedules, and
// install their weight models.
void Simulation::load_processes(double now) {
  auto make = [&](const std::vector<std::vector<Band>>& strata) {
    return std::unique_ptr<RateCells>(
      new RateCells(strata, draw_threshold, now, end)
    );
  };
  const ProcessStrata layout = process_strata();
  transition_layout = layout.transition_layout;
  transition_slots = layout.transition_slots;
  mortality = make(layout.mortality);
  search.reset();
  if (rate_initiation) {
    search = make(layout.search);
    if (search_calibrated) {
      search_calibration.reset(search.get(), calibrator->cells());
      // Stored weights (if any) pick the searcher; bin moves update the
      // calibration sums. Without a weight model every weight is one.
      if (search_weights_spec) {
        search->enable_stored_weights();
      }
      search->set_move_hook([this](int woman, int, int cell) {
        if (cell < 0) {
          search_calibration.update(woman, -1, 0, 1);
          return;
        }
        double w = 1;
        if (search_weights_spec) {
          w = search_weights_spec->custom || !search_calibration.active(woman)
                ? search_weight(woman, cell)
                : search_calibration.weight(woman);
          search->set_weight(woman, w);
        }
        search_calibration.update(
          woman, cell,
          calibrator->bin_index(people[woman].group, people[woman].born), w
        );
      });
    } else if (search_weights_spec) {
      install_weights(search_weights_spec, search.get(), -1);
    }
  }
  dissolution.reset();
  if (has_dissolution) {
    dissolution = make(layout.dissolution);
  }
  conversion.reset();
  if (cohab_enabled) {
    conversion = make(layout.conversion);
  }
  transitions.clear();
  for (size_t k = 0; k < layout.transitions.size(); ++k) {
    transitions.push_back(make(layout.transitions[k]));
    if (WeightSpec* spec =
          transition_weight_specs[layout.transition_slots[k].first]) {
      install_weights(spec, transitions.back().get(), static_cast<int>(k));
    }
  }
  if (mortality_weights_spec) {
    install_weights(mortality_weights_spec, mortality.get(), -1);
  }
  if (dissolution_weights_spec && dissolution) {
    install_weights(dissolution_weights_spec, dissolution.get(), -1);
  }
  if (conversion_weights_spec && conversion) {
    install_weights(conversion_weights_spec, conversion.get(), -1);
  }
}

// Calendar rate change: update every process's rates in place when all keep
// their cells; otherwise return false so the caller rebuilds them.
bool Simulation::update_process_rates(double now) {
  const ProcessStrata layout = process_strata();
  if (layout.transition_layout != transition_layout ||
      !mortality->same_layout(layout.mortality) ||
      (search && !search->same_layout(layout.search)) ||
      (dissolution && !dissolution->same_layout(layout.dissolution)) ||
      (conversion && !conversion->same_layout(layout.conversion)) ||
      layout.transitions.size() != transitions.size()) {
    return false;
  }
  for (size_t k = 0; k < transitions.size(); ++k) {
    if (!transitions[k]->same_layout(layout.transitions[k])) {
      return false;
    }
  }
  mortality->set_rates(layout.mortality, now);
  if (search) {
    search->set_rates(layout.search, now);
    if (search_calibrated) {
      search_calibration.mark_all_dirty();
    }
  }
  if (dissolution) {
    dissolution->set_rates(layout.dissolution, now);
  }
  if (conversion) {
    conversion->set_rates(layout.conversion, now);
  }
  for (size_t k = 0; k < transitions.size(); ++k) {
    transitions[k]->set_rates(layout.transitions[k], now);
  }
  return true;
}

// Switch matching and inheritance models at a schedule boundary (D123).
// Returns whether a matching model changed.
bool Simulation::update_models(double now) {
  if (!scheduled_models) {
    return false;
  }
  while (model_epoch + 1 < model_times.size() &&
         model_times[model_epoch + 1] <= now) {
    ++model_epoch;
  }
  const bool matching_changed = general_id != general_route[model_epoch] ||
                                father_id != father_route[model_epoch] ||
                                cohab_id != cohab_route[model_epoch];
  if (matching_changed) {
    general_id = general_route[model_epoch];
    father_id = father_route[model_epoch];
    cohab_id = cohab_route[model_epoch];
    matching = as<List>(matching_bank[general_id]);
    father_matching = as<List>(matching_bank[father_id]);
    cohab_matching = as<List>(matching_bank[cohab_id]);
    union_draw_current = get_scheduled_draw(general_id);
    father_draw = get_scheduled_draw(father_id);
    cohab_draw = get_scheduled_draw(cohab_id);
    union_weighted = applies_weights(matching);
    father_weighted = applies_weights(father_matching);
    cohab_weighted = cohab_enabled && applies_weights(cohab_matching);
    union_warm.clear();
    father_warm.clear();
    cohab_warm.clear();
    if (calibrator) {
      calibrator->set_table(union_table_id, get_scheduled_lag(general_id));
      calibrator->set_table(father_table_id, get_scheduled_lag(father_id));
      calibrator->set_table(cohab_table_id, get_scheduled_lag(cohab_id));
      calibrator->refresh();
    }
  }
  for (int d = 0; d < dimensions; ++d) {
    read_rule(
      as<List>(inheritance_bank[inheritance_routes[d][model_epoch]]),
      inheritance_rules[d]
    );
  }
  return matching_changed;
}

// Switch scheduled weight models at a calendar boundary; true if any changed.
bool Simulation::switch_weights(double now) {
  bool changed = false;
  for (auto& entry : weight_schedules) {
    WeightSchedule& w = entry.second;
    const int k = weight_interval(w, now);
    if (k < 0 || k == w.active) {
      continue;
    }
    w.active = k;
    weight_specs[entry.first] = w.models[k];
    changed = true;
  }
  return changed;
}

// Event sources in same-time priority order: death, birth, dissolution,
// conversion, search, then transition slots.
void Simulation::collect_sources() {
  sources.clear();
  sources.push_back({EventType::death, mortality.get(), -1});
  sources.push_back({EventType::birth, nullptr, -1});
  if (dissolution) {
    sources.push_back({EventType::dissolution, dissolution.get(), -1});
  }
  if (conversion) {
    sources.push_back({EventType::conversion, conversion.get(), -1});
  }
  if (search) {
    sources.push_back({EventType::search, search.get(), -1});
  }
  for (size_t k = 0; k < transitions.size(); ++k) {
    sources.push_back(
      {EventType::transition, transitions[k].get(), static_cast<int>(k)}
    );
  }
}

// GROUP STRATA AND LABELS ------

// A person's stratum in a process: their categories in its dimensions.
int Simulation::person_stratum(const StrataSpec& spec, int i) const {
  int k = 0;
  for (size_t j = 0; j < spec.dims.size(); ++j) {
    k += category(i, spec.dims[j]) * spec.strides[j];
  }
  return k;
}

// A union's stratum: the wife's (role 0) and husband's (role 1) categories.
int Simulation::union_stratum(
  const StrataSpec& spec,
  int wife,
  int husband
) const {
  int k = 0;
  for (size_t j = 0; j < spec.dims.size(); ++j) {
    k +=
      category(spec.roles[j] ? husband : wife, spec.dims[j]) * spec.strides[j];
  }
  return k;
}

// A category label of dimension d (NA when out of range).
SEXP Simulation::dim_label(int d, int category) const {
  if (d < 0 || d >= static_cast<int>(dim_labels.size()) || category < 0 ||
      category >= static_cast<int>(dim_labels[d].size())) {
    return NA_STRING;
  }
  return Rf_mkChar(dim_labels[d][category].c_str());
}

// POOLS, MATCHING AND UNION FORMATION ------

// A woman's current cohabitation probability (zero without cohabitation).
double Simulation::cohab_probability(int woman, double now) {
  if (!cohab_enabled) {
    return 0;
  }
  return current_hazard(probability_schedule(woman), people[woman].born, now)
    .first;
}

// Pool updates are constant-time (swap removal within a cell).
void Simulation::remove_male(int man) {
  if (!pool.contains(man)) {
    return;  // partnered men are already absent when they die
  }
  pool.remove(man);
  if (calibrator) {
    calibrator->add_man(people[man].group, people[man].born, -1);
  }
}

void Simulation::add_male(int man) {
  if (people[man].female || !NumericVector::is_na(people[man].died) ||
      current_partner[man] >= 0 || pool.contains(man)) {
    return;
  }
  pool.add(man, people[man].group, people[man].born);
  if (calibrator) {
    calibrator->add_man(people[man].group, people[man].born, 1);
  }
}

// weighted_choice counts living, unpartnered women. Their search and
// fertility intensities come from the cells; with cohabitation, their current
// cohabitation probabilities come from the probability tracker.
void Simulation::add_counted_woman(int woman, double now) {
  if (!calibrator || !people[woman].female ||
      !NumericVector::is_na(people[woman].died) ||
      current_partner[woman] >= 0 || counted_woman[woman]) {
    return;
  }
  counted_woman[woman] = 1;
  calibrator->add_woman(people[woman].group, people[woman].born, 1);
  if (probability_tracker) {
    probability_tracker->activate(
      woman, people[woman].born, probability_schedule(woman), now
    );
  }
}

void Simulation::remove_counted_woman(int woman) {
  if (!calibrator || !counted_woman[woman]) {
    return;
  }
  counted_woman[woman] = 0;
  if (probability_tracker) {
    probability_tracker->remove(woman);
  }
  calibrator->add_woman(people[woman].group, people[woman].born, -1);
}

// Draw models and lag tables of scheduled matching models, built on first
// use.
const DrawModel* Simulation::get_scheduled_draw(int id) {
  auto found = scheduled_draws.find(id);
  if (found == scheduled_draws.end()) {
    DrawModel draw;
    if (!people.empty()) {
      draw = make_draw_model(
        as<List>(matching_bank[id]), group_count, end - min_birth
      );
    }
    found = scheduled_draws.emplace(id, draw).first;
  }
  return &found->second;
}

// A matching model's calibration target shifts (delta), read once per model.
const std::vector<double>& Simulation::model_delta(const List& model) {
  auto found = deltas.find(SEXP(model));
  if (found == deltas.end()) {
    found = deltas.emplace(SEXP(model), as<std::vector<double>>(model["delta"]))
              .first;
  }
  return found->second;
}

const LagTable* Simulation::get_scheduled_lag(int id) {
  auto found = scheduled_lags.find(id);
  if (found == scheduled_lags.end()) {
    found = scheduled_lags
              .emplace(
                id, build_lag_table(
                      as<List>(matching_bank[id]), group_count,
                      calibrator->bins(), kWeightBinWidth
                    )
              )
              .first;
  }
  return &found->second;
}

// Partner and father searches share eligibility and the exact draw; only
// form_union changes partner state or removes the selected man. Searches and
// on-demand unions use the union-formation calibration; rate-mode nonmarital
// fathers have their own, driven by fertility. With check_cells, calibration
// totals are compared with explicit scans of the women.
int Simulation::find_candidate(
  int woman,
  double now,
  const List& choice_model,
  const DrawModel& draw,
  bool father,
  std::vector<double>* applied,
  bool cohab_role
) {
  if (pool.empty()) {
    return -1;
  }
  std::vector<double> log_weights(group_count, 0.0);
  const bool weighted =
    cohab_role ? cohab_weighted : (father ? father_weighted : union_weighted);
  if (calibrator && weighted) {
    if (probability_tracker) {
      probability_tracker->advance(now);
    }
    const bool fertility_driven = !rate_initiation || father;
    if (cells && !cohab_role && fertility_driven) {
      const std::vector<double> filled = cells->fill_calibration(*calibrator);
      if (check_cells) {
        std::vector<double> scan(filled.size(), 0.0);
        for (int w = 0; w < static_cast<int>(counted_woman.size()); ++w) {
          if (counted_woman[w]) {
            scan[calibrator->bin_index(people[w].group, people[w].born)] +=
              cells->hazard(w);
          }
        }
        for (size_t k = 0; k < scan.size(); ++k) {
          if (std::abs(scan[k] - filled[k]) >
              1e-9 * std::max(1.0, std::abs(scan[k]))) {
            stop(
              "Internal error: cell calibration totals differ from an "
              "explicit scan."
            );
          }
        }
      }
    }
    if (search_calibrated && !fertility_driven) {
      // Searches: marriage and cohabitation channels from the search cells.
      std::vector<double> marriage, cohab;
      search_calibration.apply(marriage, cohab);
      if (check_cells) {
        std::vector<double> scan_m(marriage.size(), 0.0),
          scan_c(cohab.size(), 0.0);
        for (int w = 0; w < static_cast<int>(people.size()); ++w) {
          if (!search_calibration.active(w)) {
            continue;
          }
          // Independent of stored values: p from the woman's schedule now.
          const double h = search_calibration.hazard(w);
          const double p = cohab_probability(w, now);
          const size_t b = search_calibration.bin(w);
          scan_m[b] += h * (1 - p);
          scan_c[b] += h * p;
        }
        for (size_t k = 0; k < scan_m.size(); ++k) {
          if (std::abs(scan_m[k] - marriage[k]) >
                1e-9 * std::max(1.0, std::abs(scan_m[k])) ||
              std::abs(scan_c[k] - cohab[k]) >
                1e-9 * std::max(1.0, std::abs(scan_c[k]))) {
            stop(
              "Internal error: weighted search calibration differs from an "
              "explicit scan."
            );
          }
        }
      }
      calibrator->assign_channel(0, marriage);
      calibrator->assign_channel(2, cohab);
    }
    const WeightSolution weights = calibrator->solve(
      now, cohab_role ? 2 : (fertility_driven ? 1 : 0),
      cohab_role ? cohab_table_id : (father ? father_table_id : union_table_id),
      model_delta(choice_model),
      cohab_role ? cohab_warm : (father ? father_warm : union_warm)
    );
    if (weights.infeasible) {
      ++weighted_infeasible;
    }
    log_weights = weights.log_weights;
  }
  if (applied) {
    *applied = log_weights;
  }
  return draw_partner(
    pool, draw, people[woman].group, people[woman].born, now, min_male_age,
    log_weights, [this](int person) { return people[person].born; },
    [this, woman](int man) {
      return kinship->generations() > 0 && kinship->known_related(woman, man);
    },
    draw_workspace, draw_approximate
  );
}

// Form a union for an unpartnered woman, choosing its type by her current
// cohabitation probability. Returns her partner, or -1 when none is found.
int Simulation::form_union(int woman, double now) {
  int partner = current_partner[woman];
  if (partner >= 0) {
    return partner;
  }
  std::vector<double> applied;
  const double p = cohab_probability(woman, now);
  const bool is_cohab = cohab_enabled && R::runif(0, 1) < p;
  partner = find_candidate(
    woman, now, is_cohab ? cohab_matching : matching,
    is_cohab ? *cohab_draw : *union_draw_current, false,
    diagnostics ? &applied : nullptr, is_cohab
  );
  if (partner < 0) {
    return partner;
  }
  check_union_capacity();
  const int u = static_cast<int>(union_history.size());
  union_history.push_back({++next_union_id, now, NA_REAL, woman, partner});
  union_history.back().cohab = is_cohab;
  if (diagnostics) {
    diagnostic_log_weights.push_back(applied);
  }
  person_unions[woman].push_back(u);
  person_unions[partner].push_back(u);
  current_partner[woman] = partner;
  current_partner[partner] = woman;
  remove_male(partner);
  remove_counted_woman(woman);
  return partner;
}

void Simulation::check_union_capacity() const {
  if (next_union_id >= kMaxExactId ||
      union_history.size() >=
        static_cast<size_t>(std::numeric_limits<int>::max())) {
    stop("Union history exceeds supported ID or row limits.");
  }
}

// A partnership ends (`ending` 1 dissolution, 2 widowhood): neither partner
// has a current partner, both record the ending as their previous union, and
// both return to the unpartnered pools.
void Simulation::end_partnership(
  int wife,
  int husband,
  bool cohab,
  int ending,
  double now
) {
  current_partner[wife] = -1;
  current_partner[husband] = -1;
  if (!cohab) {
    legal_previous[wife] = legal_previous[husband] = ending;
  }
  previous_union[wife] = cohab && !counts_cohab ? legal_previous[wife] : ending;
  previous_union[husband] =
    cohab && !counts_cohab ? legal_previous[husband] : ending;
  previous_union_end[wife] = previous_union_end[husband] = now;
  add_male(husband);
  add_counted_woman(wife, now);
}

// WEIGHT MODELS ------

WeightSpec* Simulation::weight_spec(const std::string& name) {
  auto found = weight_specs.find(name);
  return found == weight_specs.end() ? nullptr : &found->second;
}

// Describe a process stratum k by its dimensions: each declared dimension
// (with `_w`/`_h` for union strata) is its category's label, or NA when the
// process does not stratify by it. `skip` omits one dimension.
void Simulation::strata_bin(
  List& bin,
  const StrataSpec& spec,
  int k,
  bool union_level,
  int skip
) {
  for (int d = 0; d < dimensions; ++d) {
    if (d == skip) {
      continue;
    }
    std::vector<std::string> names =
      union_level
        ? std::vector<std::string>{dim_names[d] + "_w", dim_names[d] + "_h"}
        : std::vector<std::string>{dim_names[d]};
    for (size_t role = 0; role < names.size(); ++role) {
      SEXP value = NA_STRING;
      for (size_t j = 0; j < spec.dims.size(); ++j) {
        if (spec.dims[j] == d && static_cast<size_t>(spec.roles[j]) == role) {
          value = dim_label(
            d, (k / spec.strides[j]) % static_cast<int>(dim_labels[d].size())
          );
        }
      }
      bin[names[role]] = Rf_ScalarString(value);
    }
  }
}

// Person variable values, coded as in R's weight_variable_codes (NA is NaN).
// `dest` is a transition destination for the `destination` variable.
double Simulation::person_value(int var, int component, int i, int dest) {
  const Person& p = people[i];
  const int partner = current_partner[i];
  switch (var) {
    case 0:
      return p.female ? 1 : 0;
    case 1:
      return category(i, component);
    case 2:
      if (partner < 0) {
        return 0;
      }
      return union_history[person_unions[i].back()].cohab ? 1 : 2;
    case 3:
      return partner >= 0 ? 1 : 0;
    case 4:
      return parity[i];
    case 5:
      return previous_union[i];
    case 6:
      return partner < 0 ? R_NaN : category(partner, component);
    case 7:
      return ancestry[static_cast<size_t>(i) * ancestry_components + component];
    case 8:
      return dest;
  }
  stop("Internal error: unknown person weight variable.");
}

// Union variable values, coded as in R's weight_variable_codes.
double Simulation::union_value(int var, int component, int u) {
  const UnionRecord& x = union_history[u];
  const size_t wife = static_cast<size_t>(x.wife) * ancestry_components;
  const size_t husband = static_cast<size_t>(x.husband) * ancestry_components;
  switch (var) {
    case 2:
      return x.cohab ? 1 : 2;
    case 9:
      return category(x.wife, component);
    case 10:
      return category(x.husband, component);
    case 11:
      return parity[x.wife];
    case 12:
      return previous_union[x.wife];
    case 13:
      return previous_union[x.husband];
    case 14:
      return ancestry[wife + component];
    case 15:
      return ancestry[husband + component];
    case 16:
      return category(x.wife, component) == category(x.husband, component);
    case 17:
      // Husband's age minus wife's age.
      return people[x.wife].born - people[x.husband].born;
  }
  stop("Internal error: unknown union weight variable.");
}

// Warn once per model that missing values were treated as zero.
void Simulation::warn_weight_na(WeightSpec& spec, const std::string& what) {
  if (spec.warned_na) {
    return;
  }
  spec.warned_na = true;
  Rcpp::warning(
    "Missing values in the " + spec.process + " weight " + what +
    " were treated as a zero log-weight contribution."
  );
}

// A formula's log weight; a term with a missing variable contributes zero.
double
Simulation::formula_eta(WeightSpec& spec, int unit, bool is_union, int dest) {
  double eta = 0;
  for (size_t k = 0; k < spec.coefficients.size(); ++k) {
    double product = 1;
    bool missing = false;
    for (const auto& f : spec.factors[k]) {
      const double v = is_union ? union_value(f[0], f[2], unit)
                                : person_value(f[0], f[2], unit, dest);
      if (std::isnan(v)) {
        missing = true;
        break;
      }
      product *= f[1] == -1 ? v : (f[1] >= 0 && v == f[1] ? 1.0 : 0.0);
    }
    if (missing) {
      warn_weight_na(spec, "formula");
      continue;
    }
    eta += spec.coefficients[k] * product;
  }
  return eta;
}

// The members tibble a custom weight function receives.
List Simulation::weight_members(
  const std::vector<int>& units,
  bool is_union,
  int dest,
  const std::string& dest_label
) {
  const int n = static_cast<int>(units.size());
  List columns;
  CharacterVector names;
  auto add = [&](const std::string& name, SEXP value) {
    columns.push_back(value);
    names.push_back(name);
  };
  if (!is_union) {
    NumericVector id(n);
    CharacterVector sex(n), type(n), previous(n);
    IntegerVector partnered(n), births(n);
    std::vector<CharacterVector> own(dimensions), partners(dimensions);
    for (int d = 0; d < dimensions; ++d) {
      own[d] = CharacterVector(n);
      partners[d] = CharacterVector(n);
    }
    for (int r = 0; r < n; ++r) {
      const int i = units[r];
      const int partner = current_partner[i];
      id[r] = people[i].id;
      sex[r] = people[i].female ? "f" : "m";
      type[r] = kUnionTypeLabels[static_cast<int>(person_value(2, 0, i, -1))];
      partnered[r] = partner >= 0;
      births[r] = parity[i];
      previous[r] = kPreviousLabels[previous_union[i]];
      for (int d = 0; d < dimensions; ++d) {
        own[d][r] = dim_label(d, category(i, d));
        partners[d][r] =
          partner < 0 ? NA_STRING : dim_label(d, category(partner, d));
      }
    }
    add("person_id", id);
    add("sex", sex);
    for (int d = 0; d < dimensions; ++d) {
      add(dim_names[d], own[d]);
    }
    add("union_type", type);
    add("partnered", partnered);
    add("parity", births);
    add("previous_union", previous);
    for (int d = 0; d < dimensions; ++d) {
      add("partner_" + dim_names[d], partners[d]);
    }
    for (size_t k = 0; k < ancestry_components; ++k) {
      NumericVector share(n);
      for (int r = 0; r < n; ++r) {
        share[r] = person_value(7, static_cast<int>(k), units[r], -1);
      }
      add(weight_ancestry_labels[k], share);
    }
    if (dest >= 0) {
      CharacterVector destination(n, Rf_mkChar(dest_label.c_str()));
      add("destination", destination);
    }
  } else {
    NumericVector id(n), wife(n), husband(n), gap(n);
    CharacterVector previous_w(n), previous_h(n), type(n);
    IntegerVector parity_w(n);
    std::vector<CharacterVector> wives(dimensions), husbands(dimensions);
    std::vector<IntegerVector> same(dimensions);
    for (int d = 0; d < dimensions; ++d) {
      wives[d] = CharacterVector(n);
      husbands[d] = CharacterVector(n);
      same[d] = IntegerVector(n);
    }
    for (int r = 0; r < n; ++r) {
      const UnionRecord& x = union_history[units[r]];
      id[r] = x.id;
      wife[r] = people[x.wife].id;
      husband[r] = people[x.husband].id;
      for (int d = 0; d < dimensions; ++d) {
        wives[d][r] = dim_label(d, category(x.wife, d));
        husbands[d][r] = dim_label(d, category(x.husband, d));
        same[d][r] = category(x.wife, d) == category(x.husband, d);
      }
      parity_w[r] = parity[x.wife];
      previous_w[r] = kPreviousLabels[previous_union[x.wife]];
      previous_h[r] = kPreviousLabels[previous_union[x.husband]];
      gap[r] = people[x.wife].born - people[x.husband].born;
      type[r] = x.cohab ? "cohab" : "marriage";
    }
    add("union_id", id);
    add("wife_id", wife);
    add("husband_id", husband);
    for (int d = 0; d < dimensions; ++d) {
      add(dim_names[d] + "_w", wives[d]);
      add(dim_names[d] + "_h", husbands[d]);
    }
    add("parity_w", parity_w);
    add("previous_union_w", previous_w);
    add("previous_union_h", previous_h);
    for (const std::string suffix : {"_w", "_h"}) {
      for (size_t k = 0; k < ancestry_components; ++k) {
        NumericVector share(n);
        for (int r = 0; r < n; ++r) {
          share[r] = union_value(
            suffix == "_w" ? 14 : 15, static_cast<int>(k), units[r]
          );
        }
        add(weight_ancestry_labels[k] + suffix, share);
      }
    }
    for (int d = 0; d < dimensions; ++d) {
      add("same_" + dim_names[d], same[d]);
    }
    add("age_difference", gap);
    add("union_type", type);
  }
  columns.attr("names") = names;
  columns.attr("class") =
    CharacterVector::create("tbl_df", "tbl", "data.frame");
  columns.attr("row.names") = IntegerVector::create(NA_INTEGER, -n);
  return columns;
}

// Call a custom weight function and check its log weights. R's RNG state is
// synchronized around the call, so R code sees and continues the native
// stream.
void Simulation::custom_eta(
  WeightSpec& spec,
  List members,
  List bin,
  std::vector<double>& eta
) {
  Function f(spec.fun);
  PutRNGstate();
  NumericVector result;
  try {
    result = as<NumericVector>(f(members, bin));
  } catch (...) {
    GetRNGstate();
    throw;
  }
  GetRNGstate();
  const int n = Rf_xlength(members[0]);
  if (result.size() != n) {
    stop(
      "A custom " + spec.process + " weight function returned " +
      std::to_string(result.size()) + " values for " + std::to_string(n) +
      " members."
    );
  }
  eta.assign(n, 0.0);
  for (int r = 0; r < n; ++r) {
    if (NumericVector::is_na(result[r])) {
      warn_weight_na(spec, "function result");
      continue;
    }
    if (!R_finite(result[r])) {
      stop(
        "A custom " + spec.process +
        " weight function returned an infinite value."
      );
    }
    eta[r] = result[r];
  }
}

// Add a cell's band edges and calendar period to its description.
List Simulation::band_list(
  List bin,
  bool duration,
  double lo,
  double hi,
  DataFrame all,
  double now
) {
  bin[duration ? "dur_start" : "age_start"] = lo;
  bin[duration ? "dur_end" : "age_end"] = hi;
  const auto period = rate_period(all, now);
  bin["time_start"] = period.first;
  bin["time_end"] = period.second;
  return bin;
}

// Selection within a chosen RateCells cell by a custom weight function.
std::function<size_t(int, const std::vector<int>&, double)>
Simulation::make_chooser(WeightSpec* spec, RateCells* process, int slot) {
  return [this, spec, process, slot](
           int c, const std::vector<int>& members, double uniform
         ) -> size_t {
    const std::string& name = spec->process;
    const bool is_union = name == "dissolution" || name == "conversion";
    const int stratum = process->cell_stratum(c);
    int dest = -1;
    std::string dest_label;
    if (name == "transition") {
      const auto& at = transition_slots[slot];
      dest = tprocs[at.first].schedules[stratum][at.second].dest;
      dest_label = dim_labels[tprocs[at.first].dim][dest];
    }
    // Only custom functions choose this way; formula weights are stored.
    std::vector<double> eta(members.size());
    List bin = List::create(_["process"] = name);
    bool duration = false;
    DataFrame all = all_mortality;
    if (name == "mortality") {
      bin["sex"] = stratum % 2 ? "f" : "m";
      strata_bin(bin, mortality_spec, stratum / 2, false);
    } else if (name == "dissolution") {
      strata_bin(bin, dissolution_spec, stratum / 2, true);
      bin["union_type"] = stratum % 2 ? "cohab" : "marriage";
      duration = true;
      all = all_dissolution;
    } else if (name == "conversion") {
      strata_bin(bin, conversion_spec, stratum, true);
      duration = schedules[stratum].duration_conversion;
      all = all_conversion;
    } else if (name == "search") {
      const int g = previous_specific ? stratum / 3 : stratum;
      const int previous = previous_specific ? stratum % 3 : 0;
      strata_bin(bin, search_spec, g, false);
      bin["previous_union"] = Rf_ScalarString(
        previous_specific ? Rf_mkChar(kPreviousLabels[previous]) : NA_STRING
      );
      duration = previous_specific && schedules[g].duration_search[previous];
      all = all_search;
    } else {
      const auto& at = transition_slots[slot];
      const TransitionProcess& tp = tprocs[at.first];
      const auto& x = tp.schedules[stratum][at.second];
      const auto& labels = dim_labels[tp.dim];
      const int base = stratum % tp.spec.count;  // dependents' mirror strata
      const int origin =
        (base / tp.spec.strides[tp.own]) % static_cast<int>(labels.size());
      bin["origin"] = labels[origin];
      bin["destination"] = labels[x.dest];
      strata_bin(bin, tp.spec, base, false, tp.dim);
      bin["sex"] =
        Rf_ScalarString(x.sex < 0 ? NA_STRING : Rf_mkChar(x.sex ? "f" : "m"));
      bin["union_type"] =
        Rf_ScalarString(x.type.empty() ? NA_STRING : Rf_mkChar(x.type.c_str()));
      duration = tp.duration;
      all = tp.all;
    }
    custom_eta(
      *spec, weight_members(members, is_union, dest, dest_label),
      band_list(
        bin, duration, process->cell_start(c), process->cell_end(c), all,
        weight_now
      ),
      eta
    );
    return pick_log_weights(eta, uniform);
  };
}

// Formula models store weights (set by stored_weight at each refresh);
// custom functions choose within the fired cell.
void Simulation::install_weights(
  WeightSpec* spec,
  RateCells* process,
  int slot
) {
  if (spec->custom) {
    process->set_chooser(make_chooser(spec, process, slot));
  } else {
    process->enable_stored_weights();
  }
}

void Simulation::stored_weight(
  WeightSpec* spec,
  RateCells* process,
  int unit,
  bool is_union,
  int dest
) {
  if (!spec || spec->custom || !process) {
    return;
  }
  const double w = std::exp(formula_eta(*spec, unit, is_union, dest));
  if (!(w > 0) || !R_finite(w)) {
    stop("A " + spec->process + " weight is not finite and positive.");
  }
  process->set_weight(unit, w);
}

// The description of fertility cell c that custom functions receive.
List Simulation::fertility_bin(int c) {
  const auto info = cells->cell_info(c);
  List bin = List::create(_["process"] = "fertility");
  strata_bin(bin, fertility_spec, info.group, false);
  const std::string& type = fertility_types[info.type];
  bin["union_type"] =
    Rf_ScalarString(type.empty() ? NA_STRING : Rf_mkChar(type.c_str()));
  bin["parity_start"] = info.parity_start;
  bin["parity_end"] = info.parity_end;
  return band_list(
    bin, false, info.age_start, info.age_end, all_fertility, weight_now
  );
}

// A searching woman's stored weight in search cell c (calibration mode).
double Simulation::search_weight(int woman, int c) {
  double eta;
  if (!search_weights_spec->custom) {
    eta = formula_eta(*search_weights_spec, woman, false, -1);
  } else {
    std::vector<double> values;
    List bin = List::create(_["process"] = "search");
    const int stratum = search->cell_stratum(c);
    const int g = previous_specific ? stratum / 3 : stratum;
    const int previous = previous_specific ? stratum % 3 : 0;
    strata_bin(bin, search_spec, g, false);
    bin["previous_union"] = Rf_ScalarString(
      previous_specific ? Rf_mkChar(kPreviousLabels[previous]) : NA_STRING
    );
    custom_eta(
      *search_weights_spec, weight_members({woman}, false, -1),
      band_list(
        bin, previous_specific && schedules[g].duration_search[previous],
        search->cell_start(c), search->cell_end(c), all_search, weight_now
      ),
      values
    );
    eta = values[0];
  }
  const double w = std::exp(eta);
  if (!(w > 0) || !R_finite(w)) {
    stop("A search weight is not finite and positive.");
  }
  return w;
}

// A woman's fertility-cell state; men and the dead have no stratum.
FertilityCells::State Simulation::fertility_state(int w) {
  FertilityCells::State st;
  const Person& p = people[w];
  if (!p.female || !NumericVector::is_na(p.died)) {
    return st;
  }
  int type = 0;
  if (status_fertility && current_partner[w] >= 0) {
    type =
      typed_fertility && union_history[person_unions[w].back()].cohab ? 2 : 1;
  }
  st.stratum = cells->stratum(
    person_stratum(fertility_spec, w), type, parity_specific ? parity[w] : 0
  );
  st.group = p.group;
  st.report = person_stratum(report_spec, w);
  if (st.stratum < 0) {
    stop("No fertility schedule covers a woman's current parity.");
  }
  st.weight = heterogeneous ? fertility_weights[w] : 1;
  st.born = p.born;
  st.ready = p.last_birth + birth_interval;
  st.unpartnered = current_partner[w] < 0;
  return st;
}

// PROCESS MEMBERSHIP AND HOUSEHOLD MOVES ------

// Recompute union u's membership in the union processes (dissolution and
// conversion), whose units are unions identified by their history index.
// Called when a union starts or ends and when a partner changes group.
void Simulation::refresh_union(int u, double now) {
  const UnionRecord& x = union_history[u];
  const bool ongoing = NumericVector::is_na(x.end);
  if (dissolution) {
    // Directional strata from the wife's and husband's categories.
    const int pair = union_stratum(dissolution_spec, x.wife, x.husband);
    dissolution->refresh(
      u, ongoing ? 2 * pair + (x.cohab ? 1 : 0) : -1, x.start, now
    );
    if (ongoing) {
      stored_weight(dissolution_weights_spec, dissolution.get(), u, true, -1);
    }
  }
  if (conversion) {
    // Conversion uses the partners' categories, by duration or her age.
    const int g = union_stratum(conversion_spec, x.wife, x.husband);
    const double origin =
      schedules[g].duration_conversion ? x.start : people[x.wife].born;
    conversion->refresh(u, ongoing && x.cohab ? g : -1, origin, now);
    if (ongoing && x.cohab) {
      stored_weight(conversion_weights_spec, conversion.get(), u, true, -1);
    }
  }
}

// Recompute a person's membership in every rate process from their current
// state. Called for each person whose state an event changes.
void Simulation::refresh_person(int i, double now) {
  const Person& p = people[i];
  const bool alive = NumericVector::is_na(p.died);
  const int partner = current_partner[i];
  if (p.female) {
    cells->refresh(i, now);
  }
  mortality->refresh(
    i, alive ? 2 * person_stratum(mortality_spec, i) + (p.female ? 1 : 0) : -1,
    p.born, now
  );
  if (alive) {
    stored_weight(mortality_weights_spec, mortality.get(), i, false, -1);
  }
  if (search) {
    int stratum = -1;
    double origin = p.born;
    if (alive && p.female && partner < 0) {
      const int previous = previous_specific ? previous_union[i] : 0;
      const int g = person_stratum(search_spec, i);
      stratum = previous_specific ? 3 * g + previous : g;
      if (previous_specific && schedules[g].duration_search[previous]) {
        origin = previous_union_end[i];
        if (!R_finite(origin)) {
          stop(
            "Duration search risk requires a recorded previous union ending."
          );
        }
      }
    }
    search->refresh(i, stratum, origin, now);
    if (!search_calibrated && stratum >= 0) {
      stored_weight(search_weights_spec, search.get(), i, false, -1);
    }
    if (search_calibrated) {
      const int cell = search->cell_of(i);
      if (cell < 0) {
        search_calibration.update(i, -1, 0, 1);
      } else {
        search_calibration.set_probability(i, cohab_probability(i, now));
        double w = 1;
        if (search_weights_spec) {
          w = search_weight(i, cell);
          search->set_weight(i, w);
        }
        search_calibration.update(
          i, cell, calibrator->bin_index(p.group, p.born), w
        );
      }
    }
  }
  if ((dissolution || conversion) && partner >= 0) {
    refresh_union(person_unions[i].back(), now);
  }
  if (transitions.empty()) {
    return;
  }
  // Union-status code as in TransitionSchedule::code.
  const int status =
    partner < 0 ? 0 : (union_history[person_unions[i].back()].cohab ? 1 : 2);
  for (size_t f = 0; f < transitions.size(); ++f) {
    const auto& at = transition_slots[f];
    const TransitionProcess& tp = tprocs[at.first];
    const size_t k = at.second;
    // Household moves: one clock per household, carried by the woman in a
    // union; dependents use the mirror strata (no risk before the dependent
    // age).
    const int stratum =
      person_stratum(tp.spec, i) +
      (tp.household && is_dependent(tp, i, now) ? tp.spec.count : 0);
    const auto& outgoing = tp.schedules[stratum];
    bool at_risk = alive && k < outgoing.size() &&
                   !(tp.household && !p.female && partner >= 0);
    if (at_risk) {
      const auto& x = outgoing[k];
      if (x.sex >= 0 && x.sex != (p.female ? 1 : 0)) {
        at_risk = false;
      }
      if (x.code >= 0 && x.code != status && !(x.code == 3 && partner >= 0)) {
        at_risk = false;
      }
    }
    transitions[f]->refresh(
      i, at_risk ? stratum : -1, tp.duration ? tp.entry[i] : p.born, now
    );
    if (at_risk) {
      stored_weight(
        transition_weight_specs[at.first], transitions[f].get(), i, false,
        outgoing[k].dest
      );
    }
  }
}

// Refresh every living person (after processes are rebuilt or weight models
// change). The dead hold no process membership.
void Simulation::refresh_living(double now) {
  for (int i = 0; i < static_cast<int>(people.size()); ++i) {
    if (i % 1024 == 0) {
      checkUserInterrupt();
    }
    if (NumericVector::is_na(people[i].died)) {
      refresh_person(i, now);
    }
  }
}

// A living child's attached parent for household moves: the living mother,
// otherwise the living father (-1: independent).
int Simulation::attached_parent(int c) const {
  const int m = mother_of[c], f = father_of[c];
  if (m >= 0 && NumericVector::is_na(people[m].died)) {
    return m;
  }
  if (f >= 0 && NumericVector::is_na(people[f].died)) {
    return f;
  }
  return -1;
}

// A dependent: alive, unpartnered, below the process's dependent age and
// attached to a living parent.
bool Simulation::is_dependent(
  const TransitionProcess& tp,
  int c,
  double now
) const {
  return NumericVector::is_na(people[c].died) && current_partner[c] < 0 &&
         people[c].born + tp.dependent_age > now && attached_parent(c) >= 0;
}

// Household moves (D140). A party grows from its heads by adding the
// dependents attached to any member, recursively (a dependent's own child
// follows them).
std::vector<int> Simulation::household_party(
  const TransitionProcess& tp,
  std::vector<int> party,
  double now
) {
  for (size_t h = 0; h < party.size(); ++h) {
    for (int c : kids[party[h]]) {
      if (is_dependent(tp, c, now) && attached_parent(c) == party[h] &&
          std::find(party.begin(), party.end(), c) == party.end()) {
        party.push_back(c);
      }
    }
  }
  return party;
}

// Move a party to a destination in tp's dimension as one recorded move;
// members already there stay and are not recorded.
void Simulation::move_party(
  TransitionProcess& tp,
  const std::vector<int>& party,
  int destination,
  double initiator,
  double now
) {
  ++move_count;
  std::vector<int> moved;
  for (int m : party) {
    const size_t cell = static_cast<size_t>(m) * dimensions + tp.dim;
    const int origin = cats[cell];
    if (origin == destination) {
      continue;
    }
    if (people[m].female) {
      remove_counted_woman(m);
    } else {
      remove_male(m);
    }
    tx_person.push_back(people[m].id);
    tx_time.push_back(now);
    tx_dim.push_back(tp.dim);
    tx_origin.push_back(origin);
    tx_dest.push_back(destination);
    tx_move.push_back(move_count);
    tx_initiator.push_back(initiator);
    cats[cell] = destination;
    if (matching_dim[tp.dim]) {
      people[m].group = person_stratum(matching_spec, m);
    }
    tp.entry[m] = now;
    if (people[m].female) {
      add_counted_woman(m, now);
    } else {
      add_male(m);
    }
    moved.push_back(m);
  }
  // A move changes union strata and the partner's partner variables.
  for (int m : moved) {
    const int partner = current_partner[m];
    refresh_person(m, now);
    if (partner >= 0) {
      refresh_person(partner, now);
    }
  }
}

// A dependent whose attachment changed (a mother's death, the end of their
// own union) joins the parent they now depend on, with their own dependents,
// recorded as a move initiated by that parent.
void Simulation::rejoin(int c, double now) {
  if (!any_household || !NumericVector::is_na(people[c].died)) {
    return;
  }
  for (auto& tp : tprocs) {
    if (!tp.household || !is_dependent(tp, c, now)) {
      continue;
    }
    const int parent = attached_parent(c);
    const int destination = category(parent, tp.dim);
    if (category(c, tp.dim) == destination) {
      continue;
    }
    move_party(
      tp, household_party(tp, {c}, now), destination, people[parent].id, now
    );
  }
}

// EVENTS ------

// Process events in chronological order until no scheduled event remains
// before the horizon.
List Simulation::run() {
  while (true) {
    if (apply_calendar_boundary()) {
      continue;
    }
    // Cell intensities depend on current membership, so bookkeeping
    // boundaries are processed in time order before any event at or after
    // them. Exact event ties follow the source priority order.
    double boundary = R_PosInf, next_event = R_PosInf;
    size_t boundary_source = 0, event_source = 0;
    for (size_t k = 0; k < sources.size(); ++k) {
      const double b = source_boundary(sources[k]);
      const double c = source_clock(sources[k]);
      if (b < boundary) {
        boundary = b;
        boundary_source = k;
      }
      if (c < next_event) {
        next_event = c;
        event_source = k;
      }
    }
    if (boundary < end && boundary <= next_event) {
      weight_now = boundary;
      const EventSource& x = sources[boundary_source];
      if (x.cells) {
        x.cells->advance_boundary();
      } else {
        cells->advance_boundary();
      }
      continue;
    }
    if (!(next_event < end)) {
      break;
    }
    if (++processed % 1024 == 0) {
      checkUserInterrupt();
    }
    const double now = next_event;
    weight_now = now;
    const EventSource fired = sources[event_source];
    // Union processes return a union index; their events act on its wife.
    const int unit =
      fired.cells ? fired.cells->fire_next(now) : cells->fire_next(now);
    if (unit < 0) {
      continue;  // an unrealized birth: no event follows
    }
    const bool union_event = fired.type == EventType::dissolution ||
                             fired.type == EventType::conversion;
    const int i = union_event ? union_history[unit].wife : unit;
    if (probability_tracker) {
      probability_tracker->advance(now);
    }
    if (check_cells && cells && !cells->consistent()) {
      stop(
        "Internal error: fertility cell counts differ from an explicit scan."
      );
    }
    if (check_cells) {
      for (const EventSource& x : sources) {
        if (x.cells && !x.cells->consistent()) {
          stop(
            "Internal error: rate cell counts differ from an explicit scan."
          );
        }
      }
    }
    switch (fired.type) {
      case EventType::transition:
        on_transition(i, fired.slot, now);
        break;
      case EventType::death:
        on_death(i, now);
        break;
      case EventType::dissolution:
        on_dissolution(i, unit, now);
        break;
      case EventType::conversion:
        on_conversion(i, unit, now);
        break;
      case EventType::search:
        on_search(i, now);
        break;
      case EventType::birth:
        on_birth(i, now);
        break;
    }
  }
  return results();
}

// Apply the next calendar boundary if it precedes every pending boundary and
// event; calendar changes win exact ties. Returns whether one was applied.
// Conditional Exp(1) resampling is valid after survival to the boundary, and
// candidates under superseded rates never execute beyond it.
bool Simulation::apply_calendar_boundary() {
  if (calendar_boundaries.empty()) {
    return false;
  }
  const double now = *calendar_boundaries.begin();
  for (const EventSource& x : sources) {
    if (!(now <= source_boundary(x) && now <= source_clock(x))) {
      return false;
    }
  }
  calendar_boundaries.erase(calendar_boundaries.begin());
  weight_now = now;
  checkUserInterrupt();
  const bool matching_changed = update_models(now);
  const bool weights_changed = switch_weights(now);
  const bool probability_changed = calendar_changes_at(all_probability, now);
  const bool fertility_change = calendar_changes_at(all_fertility, now);
  const bool rate_change =
    calendar_changes_at(all_mortality, now) || fertility_change ||
    calendar_changes_at(all_search, now) ||
    calendar_changes_at(all_dissolution, now) ||
    calendar_changes_at(all_conversion, now) ||
    std::any_of(tprocs.begin(), tprocs.end(), [&](const TransitionProcess& tp) {
      return calendar_changes_at(tp.all, now);
    });
  if (!rate_change && !matching_changed && !probability_changed &&
      !weights_changed) {
    return true;
  }
  for (int w = 0; w < static_cast<int>(counted_woman.size()); ++w) {
    remove_counted_woman(w);
  }
  if (probability_changed) {
    cohab_prob = calendar_slice(all_probability, now);
    load_probabilities();
  }
  if (!rate_change) {
    // Model/probability changes retain all demographic event clocks. Role
    // preferences and intensities use the newly active inputs.
    for (int w = 0; w < static_cast<int>(people.size()); ++w) {
      if (w % 1024 == 0) {
        checkUserInterrupt();
      }
      add_counted_woman(w, now);
    }
    // New weight models: recompute every stored weight. Intensities and
    // event clocks are unchanged, so no draws are needed.
    if (weights_changed) {
      refresh_living(now);
    }
    return true;
  }
  if (fertility_change) {
    cells_consistent = cells_consistent && cells->consistent();
  }
  for (const EventSource& x : sources) {
    if (x.cells) {
      cells_consistent = cells_consistent && x.cells->consistent();
    }
  }
  slice_rates(now);
  for (auto& tp : tprocs) {
    tp.active = calendar_slice(tp.all, now);
  }
  load_schedules();
  load_dissolution();
  load_transitions();
  // Rates that keep every cell change in place: clocks continue and nobody
  // moves. A changed layout rebuilds that part, and its new thresholds are
  // valid after survival to the boundary (memoryless).
  bool rebuilt = false;
  if (fertility_change) {
    const bool same = has_cohab_fertility(fertility_rates) == typed_fertility &&
                      fertility_rates.containsElementNamed(".hammel_group") ==
                        fertility_grouped &&
                      cells->update_rates(
                        fertility_strata(
                          fertility_rates, fertility_spec.count, fertility_types
                        ),
                        now
                      );
    if (!same) {
      load_cells(now);
      rebuilt = true;
    }
  }
  if (!update_process_rates(now)) {
    load_processes(now);
    collect_sources();
    rebuilt = true;
  }
  // All table objects referenced by the tracker keep their addresses. Clock
  // selection is fixed across calendar periods by R validation.
  for (int w = 0; w < static_cast<int>(people.size()); ++w) {
    add_counted_woman(w, now);
  }
  // Rebuilt processes start empty; only the living (and their ongoing
  // unions) re-enter. A new weight model also needs every stored weight
  // recomputed, even when rates changed in place.
  if (rebuilt || weights_changed) {
    refresh_living(now);
  }
  return true;
}

// A group transition moves the person, or for household moves (D140) their
// partner and, recursively, the dependents attached to party members.
void Simulation::on_transition(int i, int slot, double now) {
  const auto& at = transition_slots[slot];
  TransitionProcess& tp = tprocs[at.first];
  const int destination =
    tp.schedules[person_stratum(tp.spec, i)][at.second].dest;
  if (destination == category(i, tp.dim)) {
    stop("Invalid transition destination.");
  }
  std::vector<int> heads = {i};
  if (tp.household && current_partner[i] >= 0) {
    heads.push_back(current_partner[i]);
  }
  move_party(
    tp, tp.household ? household_party(tp, heads, now) : heads, destination,
    people[i].id, now
  );
}

// Retain the person's history, record death, widow any partner, and leave
// every process.
void Simulation::on_death(int i, double now) {
  people[i].died = now;
  std::vector<int> widowed, ended;
  for (int u : person_unions[i]) {
    UnionRecord& relationship = union_history[u];
    if (!NumericVector::is_na(relationship.end)) {
      continue;
    }
    relationship.end = now;
    relationship.wife_died = i == relationship.wife;
    end_partnership(
      relationship.wife, relationship.husband, relationship.cohab, 2, now
    );
    widowed.push_back(
      i == relationship.wife ? relationship.husband : relationship.wife
    );
    ended.push_back(u);
  }
  if (!people[i].female) {
    remove_male(i);
  } else {
    remove_counted_woman(i);
  }
  refresh_person(i, now);
  for (int survivor : widowed) {
    refresh_person(survivor, now);
  }
  for (int u : ended) {
    refresh_union(u, now);
  }
  // Children may now be attached to the other parent (joining them), or
  // independent; a widowed survivor under the age of independence may be a
  // dependent again.
  if (any_household) {
    for (int c : kids[i]) {
      if (NumericVector::is_na(people[c].died)) {
        refresh_person(c, now);
        rejoin(c, now);
      }
    }
    for (int survivor : widowed) {
      rejoin(survivor, now);
    }
  }
  ++death_count;
}

void Simulation::on_dissolution(int i, int unit, double now) {
  const int partner = current_partner[i];
  if (partner < 0) {
    stop("Dissolution requires an ongoing union.");
  }
  UnionRecord& relationship = union_history[unit];
  relationship.end = now;
  relationship.dissolved = true;
  end_partnership(i, partner, relationship.cohab, 1, now);
  ++dissolution_count;
  refresh_union(unit, now);
  refresh_person(i, now);
  refresh_person(partner, now);
  // A partner under the age of independence may be a dependent again.
  rejoin(i, now);
  rejoin(partner, now);
}

// A cohabitation converts to marriage: the episode closes and a linked
// marriage record opens at the same time.
void Simulation::on_conversion(int i, int unit, double now) {
  const int partner = current_partner[i];
  if (partner < 0 || !union_history[unit].cohab) {
    stop("Conversion requires an ongoing cohabitation.");
  }
  union_history[unit].end = now;
  union_history[unit].converted = true;
  const double prior_id = union_history[unit].id;
  check_union_capacity();
  const int next = static_cast<int>(union_history.size());
  union_history.push_back({++next_union_id, now, NA_REAL, i, partner});
  union_history.back().converted_from = prior_id;
  person_unions[i].push_back(next);
  person_unions[partner].push_back(next);
  if (diagnostics) {
    diagnostic_log_weights.emplace_back(group_count, NA_REAL);
  }
  ++conversion_count;
  refresh_union(unit, now);  // the cohabitation closes; the marriage enters
  refresh_person(i, now);
  refresh_person(partner, now);
}

void Simulation::on_search(int i, double now) {
  const int partner = form_union(i, now);
  refresh_person(i, now);
  if (partner >= 0) {
    refresh_person(partner, now);
  }
}

// A birth to woman i. On-demand unions form first; a rate-mode nonmarital
// father is drawn without forming a union and stays available. Without a
// father the birth keeps an unknown father, unless births are cancelled.
void Simulation::on_birth(int i, double now) {
  double father = NA_REAL;
  int partner = current_partner[i];
  if (!rate_initiation && partner < 0) {
    partner = form_union(i, now);
    if (partner >= 0) {
      refresh_person(partner, now);
    }
  }
  int father_position = partner;
  if (partner < 0 && rate_initiation) {
    father_position =
      find_candidate(i, now, father_matching, *father_draw, true, nullptr);
  }
  if (father_position < 0 && cancel_birth_if_no_match) {
    // No child was born: retain last_birth; her union state may have changed.
    refresh_person(i, now);
    return;
  }
  if (father_position >= 0) {
    father = people[father_position].id;
  }
  if (next_id >= kMaxExactId ||
      people.size() >= static_cast<size_t>(std::numeric_limits<int>::max())) {
    stop("Population exceeds supported ID or row limits.");
  }
  // IDs increase beyond existing histories. Union formation above happens
  // before the newborn exists.
  const double mother = people[i].id;
  people[i].last_birth = now;
  ++parity[i];
  // Default ffab is 0.4886 - its the signature of Berkeley training!
  const bool female = R::runif(0, 1) < ffab;
  const int child = static_cast<int>(people.size());
  people.push_back(
    {++next_id, now, NA_REAL, mother, father,
     -std::numeric_limits<double>::infinity(), female, 0}
  );
  mother_of.push_back(i);
  father_of.push_back(father_position);
  kids.emplace_back();
  kids[i].push_back(child);
  if (father_position >= 0) {
    kids[father_position].push_back(child);
  }
  parity.push_back(0);
  if (heterogeneous) {
    fertility_weights.push_back(
      female ? inherit_fertility_weight(
                 fertility_weights[i], fertility_variance, fertility_alpha
               )
             : NA_REAL
    );
  }
  append_ancestry(child, i, father_position);
  // Sex and ancestry are now available; assign each dimension's category, in
  // declaration order, before any newborn hazard or matching lookup.
  cats.resize(static_cast<size_t>(child + 1) * dimensions);
  for (int d = 0; d < dimensions; ++d) {
    cats[static_cast<size_t>(child) * dimensions + d] =
      newborn_category(d, child, i, father_position, female, partner >= 0);
  }
  const int child_group = person_stratum(matching_spec, child);
  people[child].group = child_group;
  kinship->append(i, father_position);
  if (!female) {
    pool.add(child, child_group, now);
    if (calibrator) {
      calibrator->add_man(child_group, now, 1);
    }
  }
  person_unions.emplace_back();
  current_partner.push_back(-1);
  previous_union.push_back(0);
  legal_previous.push_back(0);
  previous_union_end.push_back(NA_REAL);
  if (calibrator) {
    counted_woman.push_back(0);
    add_counted_woman(child, now);
  }
  for (auto& tp : tprocs) {
    tp.entry.push_back(now);
  }
  ++birth_count;
  // The newborn enters its cells; the mother's spacing (and possibly her
  // union state, in on-demand mode) has changed.
  refresh_person(child, now);
  refresh_person(i, now);
}

// A newborn's ancestry is the average of the parents'. A missing father
// contributes his half to each dimension's unknown share; supplied parental
// unknown shares transmit through the same rule.
void Simulation::append_ancestry(int child, int mother, int father) {
  const size_t offset = static_cast<size_t>(child) * ancestry_components;
  const size_t maternal = static_cast<size_t>(mother) * ancestry_components;
  ancestry.resize(offset + ancestry_components);
  for (size_t g = 0; g < ancestry_components; ++g) {
    const double paternal =
      father >= 0
        ? ancestry[static_cast<size_t>(father) * ancestry_components + g]
        : (unknown_component[g] ? 1.0 : 0.0);
    ancestry[offset + g] = 0.5 * ancestry[maternal + g] + 0.5 * paternal;
  }
}

// A newborn's category in dimension d under that dimension's active rule.
// Transitions separately update existing identities.
int Simulation::newborn_category(
  int d,
  int child,
  int mother,
  int father,
  bool female,
  bool partnered
) {
  const InheritanceRule& rule = inheritance_rules[d];
  const size_t base = rule.offset, size = rule.size;
  const double* shares =
    ancestry.data() + static_cast<size_t>(child) * ancestry_components + base;
  const int maternal_group = category(mother, d);
  const int paternal_group = father < 0 ? -1 : category(father, d);
  const bool mixed = paternal_group >= 0 && paternal_group != maternal_group;
  int child_group = maternal_group;
  if (rule.method == 8) {
    child_group = rule.fixed;
  } else if (rule.method == 1 && paternal_group >= 0) {
    child_group = paternal_group;
  } else if (rule.method == 2 && mixed) {
    child_group = R::runif(0, 1) < 0.5 ? maternal_group : paternal_group;
  } else if (rule.method == 3 && !female && paternal_group >= 0) {
    child_group = paternal_group;
  } else if (rule.method == 4 && female && paternal_group >= 0) {
    child_group = paternal_group;
  } else if (rule.method == 5 && size > 0 && shares[size - 1] < 1.0) {
    // Ancestry draw: a group with probability equal to its known share.
    double total = 0;
    int positive = 0, sole_group = maternal_group;
    for (size_t g = 0; g + 1 < size; ++g) {
      if (shares[g] > 0) {
        if (rule.ancestry_groups[g] < 0) {
          stop("Positive ancestry has no declared identity group.");
        }
        total += shares[g];
        ++positive;
        sole_group = rule.ancestry_groups[g];
      }
    }
    if (positive == 1) {
      child_group = sole_group;
    } else if (positive > 1) {
      const double draw = R::runif(0, total);
      double cumulative = 0;
      child_group = sole_group;  // protect the final cell against rounding
      for (size_t g = 0; g + 1 < size; ++g) {
        cumulative += shares[g];
        if (draw < cumulative) {
          child_group = rule.ancestry_groups[g];
          break;
        }
      }
    }
  } else if (rule.method == 6 && mixed) {
    // Ranked: the higher-ranked parent's group, reversed with probability
    // `passing`.
    const int maternal_rank = rule.ranks[maternal_group];
    const int paternal_rank = rule.ranks[paternal_group];
    if (maternal_rank >= 0 && paternal_rank >= 0) {
      const bool reverse = rule.passing == 1 ||
                           (rule.passing > 0 && R::runif(0, 1) < rule.passing);
      const bool maternal_first = maternal_rank < paternal_rank;
      child_group = maternal_first != reverse ? maternal_group : paternal_group;
    }
  } else if (rule.method == 7) {
    NumericVector maternal(size), paternal(size), child_shares(size);
    for (size_t g = 0; g < size; ++g) {
      maternal[g] =
        ancestry[static_cast<size_t>(mother) * ancestry_components + base + g];
      if (father >= 0) {
        paternal[g] = ancestry
          [static_cast<size_t>(father) * ancestry_components + base + g];
      }
      child_shares[g] = shares[g];
    }
    // R's random functions must see all preceding native draws, and native
    // scheduling must resume from the state left by the R callback.
    PutRNGstate();
    try {
      Function inherit_custom(rule.custom);
      child_group = as<int>(inherit_custom(
        maternal_group, paternal_group, female ? "f" : "m", maternal, paternal,
        child_shares, partnered
      ));
    } catch (...) {
      GetRNGstate();
      throw;
    }
    GetRNGstate();
  }
  if (child_group < 0 || child_group >= rule.count) {
    stop("Invalid newborn group index.");
  }
  return child_group;
}

// RESULTS ------

// Materialize ordinary R columns from the engine records. The R wrapper
// restores custom columns, closes existing unions affected by deaths, and
// assembles the final population and run metadata. Pending events are not
// saved as a checkpoint.
List Simulation::results() {
  // Unrealized births by one-year period from the run start and cell.
  const size_t nu = unrealized_counts.size();
  NumericVector period_start(nu), period_end(nu), parity_start(nu),
    parity_end(nu), age_start(nu), age_end(nu), unrealized_births(nu);
  IntegerVector unrealized_group(nu);
  CharacterVector unrealized_type(nu);
  size_t r = 0;
  for (const auto& entry : unrealized_counts) {
    const auto& key = entry.first;
    const double grid_start = start + static_cast<double>(std::get<0>(key));
    period_start[r] = std::max(start, grid_start);
    period_end[r] = std::min(end, grid_start + 1);
    unrealized_group[r] = std::get<1>(key);
    unrealized_type[r] = std::get<2>(key);
    parity_start[r] = std::get<3>(key);
    parity_end[r] = std::get<4>(key);
    age_start[r] = std::get<5>(key);
    age_end[r] = std::get<6>(key);
    unrealized_births[r] = entry.second;
    ++r;
  }
  DataFrame unrealized = DataFrame::create(
    _["period_start"] = period_start, _["period_end"] = period_end,
    _["group"] = unrealized_group, _["union_type"] = unrealized_type,
    _["parity_start"] = parity_start, _["parity_end"] = parity_end,
    _["age_start"] = age_start, _["age_end"] = age_end,
    _["unrealized_births"] = unrealized_births
  );
  if (cells) {
    cells_consistent = cells_consistent && cells->consistent();
  }

  const size_t n = people.size();
  NumericVector out_id(n), out_birth(n), out_death(n), out_mother(n),
    out_father(n);
  CharacterVector out_sex(n);
  IntegerVector out_group(n);
  IntegerMatrix out_categories(n, dimensions);
  NumericMatrix out_ancestry(n, ancestry_components);
  for (size_t i = 0; i < n; ++i) {
    const Person& p = people[i];
    out_id[i] = p.id;
    out_birth[i] = p.born;
    out_death[i] = p.died;
    out_mother[i] = p.mother;
    out_father[i] = p.father;
    out_sex[i] = p.female ? "f" : "m";
    out_group[i] = p.group;
    for (int d = 0; d < dimensions; ++d) {
      out_categories(i, d) = category(static_cast<int>(i), d);
    }
    for (size_t g = 0; g < ancestry_components; ++g) {
      out_ancestry(i, g) = ancestry[i * ancestry_components + g];
    }
  }

  const size_t m = union_history.size();
  NumericVector out_uid(m), out_wife(m), out_husband(m), out_start(m),
    out_end(m), out_from(m);
  CharacterVector out_reason(m), out_type(m);
  CharacterVector input_reasons = unions["end_reason"];
  NumericVector input_ends = unions["end_time"];
  for (size_t u = 0; u < m; ++u) {
    const UnionRecord& x = union_history[u];
    out_type[u] = x.cohab ? "cohab" : "marriage";
    out_from[u] = x.converted_from;
    out_uid[u] = x.id;
    out_wife[u] = people[x.wife].id;
    out_husband[u] = people[x.husband].id;
    out_start[u] = x.start;
    out_end[u] = x.end;
    if (NumericVector::is_na(x.end)) {
      out_reason[u] = NA_STRING;
    } else if (x.converted) {
      out_reason[u] = "conversion";
    } else if (x.dissolved) {
      out_reason[u] = "dissolution";
    } else if (u < static_cast<size_t>(unions.nrows()) &&
               !NumericVector::is_na(input_ends[u])) {
      out_reason[u] = input_reasons[u];  // ended in the input history
    } else {
      out_reason[u] = x.wife_died ? "death_wife" : "death_husband";
    }
  }
  DataFrame out_unions =
    cohab_enabled || unions.containsElementNamed("union_type")
      ? DataFrame::create(
          _["union_id"] = out_uid, _["wife_id"] = out_wife,
          _["husband_id"] = out_husband, _["start_time"] = out_start,
          _["end_time"] = out_end, _["end_reason"] = out_reason,
          _["union_type"] = out_type, _["previous_cohab_id"] = out_from
        )
      : DataFrame::create(
          _["union_id"] = out_uid, _["wife_id"] = out_wife,
          _["husband_id"] = out_husband, _["start_time"] = out_start,
          _["end_time"] = out_end, _["end_reason"] = out_reason
        );

  // Rows are new unions in order; columns are group indices. Unions formed
  // without weights (unit weights) record zeros.
  NumericMatrix log_weights(
    static_cast<int>(diagnostic_log_weights.size()), group_count
  );
  for (size_t u = 0; u < diagnostic_log_weights.size(); ++u) {
    for (int g = 0; g < group_count; ++g) {
      log_weights(u, g) =
        diagnostic_log_weights[u].empty() ? 0 : diagnostic_log_weights[u][g];
    }
  }
  return List::create(
    _["persons"] = DataFrame::create(
      _["person_id"] = out_id, _["sex"] = out_sex, _["birth_time"] = out_birth,
      _["death_time"] = out_death, _["mother_id"] = out_mother,
      _["father_id"] = out_father
    ),
    _["transitions"] = DataFrame::create(
      _["person_id"] = wrap(tx_person), _["time"] = wrap(tx_time),
      _["dimension"] = wrap(tx_dim), _["origin"] = wrap(tx_origin),
      _["destination"] = wrap(tx_dest), _["move_id"] = wrap(tx_move),
      _["initiator_id"] = wrap(tx_initiator)
    ),
    _["groups"] = out_group, _["categories"] = out_categories,
    _["unions"] = out_unions, _["ancestry"] = out_ancestry,
    _["births"] = birth_count, _["deaths"] = death_count,
    _["dissolutions"] = dissolution_count, _["conversions"] = conversion_count,
    _["weighted_infeasible"] = weighted_infeasible,
    _["draw_approximate"] = draw_approximate, _["unrealized"] = unrealized,
    _["cells_consistent"] = cells_consistent,
    _["hetfert"] = wrap(fertility_weights),
    _["diagnostic_log_weights"] = log_weights
  );
}

}  // namespace

// ENTRY POINT ------

// Run one simulation. run_simulation() validates and prepares the inputs in
// R and handles the seed.
// [[Rcpp::export]]
List simulate_cpp(
  DataFrame persons,
  DataFrame mortality_rates,
  DataFrame fertility_rates,
  double start,
  double years,
  double birth_interval,
  DataFrame unions,
  List matching,
  double min_male_age,
  int kinship_generations,
  bool rate_initiation,
  DataFrame union_rates,
  List father_matching,
  bool cancel_birth_if_no_match,
  bool diagnostics,
  NumericMatrix initial_ancestry,
  List inheritance,
  double ffab,
  DataFrame dissolution_rates,
  DataFrame cohab_prob,
  DataFrame conversion_rates,
  List cohab_matching,
  bool counts_cohab,
  bool check_cells = false,
  Nullable<List> hetfert = R_NilValue,
  Nullable<List> transition_input = R_NilValue,
  Nullable<List> model_schedule_input = R_NilValue,
  Nullable<List> weight_input = R_NilValue,
  Nullable<List> group_input = R_NilValue
) {
  Simulation simulation(
    persons, mortality_rates, fertility_rates, start, years, birth_interval,
    unions, matching, min_male_age, kinship_generations, rate_initiation,
    union_rates, father_matching, cancel_birth_if_no_match, diagnostics,
    initial_ancestry, inheritance, ffab, dissolution_rates, cohab_prob,
    conversion_rates, cohab_matching, counts_cohab, check_cells, hetfert,
    transition_input, model_schedule_input, weight_input, group_input
  );
  return simulation.run();
}

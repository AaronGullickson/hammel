# Analytical investigation of column-power matching adjustments

Status: exploratory calculation authorized October 22, 2026. This does not
change Hammel's engine or accept a demographic matching mechanism. D062's
additive selection remains implemented; D063's geometric experiment was removed.

## Calculation specification

Use a fixed group-only market, with positive female and male group shares
f[i] and m[j]. Each sex's shares sum to one. A woman in group i initiates
with probability f[i], equivalent to equal initiation per woman. There are no
ages, event times, demographic transitions, sampled markets, or hard exclusions.
Pair utilities U[i,j] are dimensionless log weights, with female-first group
orientation. Their exponentials W[i,j] are positive matching weights.

Define each man's weighted female availability as
C[j] = sum_i f[i] W[i,j]. The proposed rule is
P[j | i] proportional to m[j] W[i,j] / C[j]^alpha.
Normalize separately for each female group, then calculate expected pairing
shares T[i,j] = f[i] P[j | i]. These are exact expectations for one choice;
there is no Monte Carlo simulation. Pool-size factors common to every column
cancel. No target male margins are imposed.

Compare alpha = 0 (female choice), 0.5 (geometric), and one constant exponent
calibrated to the previous symmetric 75/25, cross-utility -3 example. Separately
find the balancing alpha within each symmetric two-group case to investigate
whether it has a general form. Never calibrate directional equality for
asymmetric sex compositions or directed preferences.

Validation: probability sums and female row margins; fixed-table cross-product
log odds equal U[A,B]+U[B,A]-U[A,A]-U[B,B] for every exponent; agreement with
expanded individual-level matrices; reproduction of the original additive
example; numerical/closed-form agreement for symmetric two-group roots.
Equal group shares or zero preference make the balancing exponent unidentified:
all exponents balance in those cases. Three-group cases assess all pairwise
flows, not just their total. Directed and sex-asymmetric examples are diagnostic
outputs; unequal flows are permitted, not treated as failures.

The analysis cannot establish recovery in age-structured or changing markets,
or determine which male participation assumptions are demographically right.

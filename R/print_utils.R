# Functions useful for printing methods

# HELPER FUNCTIONS ------------------------------------------------------------

# format the printing of a stat
print_stat <- function(label, value) {
  cat(sprintf("%-28s %12s\n", label, value))
}


with_database_reader <- function(name, replacement, code) {
  bindings <- setNames(list(replacement), name)
  do.call(local_mocked_bindings, c(bindings, list(.package = "HMDHFDplus", .env = environment())))
  force(code)
}

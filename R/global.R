.scrattch.io_env <- new.env(parent = emptyenv())

.scrattch.io_env$overwrite <- FALSE

#' Set global overwrite parameter
#'
#' This parameter will affect the default behavior of all write_tome functions, but will be overridden by parameters passed to each function.
#'
#' @param overwrite logical, the setting of the default overwrite parameter for all scrattch.io functions.
#' Default = FALSE.
#'
set_scrattch.io_global_overwrite <- function(overwrite = FALSE) {
  .scrattch.io_env$overwrite <- overwrite
}

.scrattch.io_env$verbosity <- 2

#' Set global verbosity parameter
#'
#' This parameter will affect the default behavior of all scrattch.io read or write functions that return print or cat statements
#'
#' @param verbosity integer, the verbosity level to use.
#'
#' There are 3 levels:
#' 0 removes all messages.
#' 1 provides TRUE or FALSE for success or failure of reads/writes.
#' 2 provides useful feedback for interactive use.
#' Default is 2.
#'
set_scrattch.io_global_verbosity <- function(verbosity = 2) {
  .scrattch.io_env$verbosity <- verbosity
}

#' Detects the number of cores that can be used for benchmarking.
#' 
#' Uses all available resources by default
#' 
#' @importFrom parallel detectCores
#' @export
detect_cores <- function() {
  n_cores <- detectCores()
  if (n_cores < 2) stop("At least two cores are required for benchmarking")
  c(floor(n_cores / 2), 2)
}

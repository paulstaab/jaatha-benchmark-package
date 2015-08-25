#' Excutes Jaathas self test
#' @param seed The Random Seed to use.
#' @param cores The number of cores to use. This is a numeric vector of length
#'   two, where the first entry is the number of Jaatha runs execute in parallel,
#'   and the second is the number of cores that are used in each run.
#' @param folder The folder where the simulation results are saved. Must not exists
#'   when calling the function.
#' @param test_data If specified, use the test data generated with 
#'   \code{\link{createTestData}} instead of generating it in the beginning.
#'   The arguments \code{n.points} and \code{reps} are ignored then.
#' @param ... Additional argument for \code{jaatha}.
#' @inheritParams createTestData
#' 
#' @return Nothing.
#' @importFrom jaatha jaatha create_jaatha_model create_jaatha_data
#' @importFrom evaluate try_capture_stack
#' @export
#' @examples
#' library('jaatha')
#' dm <- coala:::model_theta_tau()
#' testJaatha(dm, 2, 1, folder = tempfile(), scaling_factor=10)
#' 
#' test_data <- createTestData(dm, 2, 1, grid.pars=2)
#' testJaatha(dm, test_data = test_data, folder=tempfile())
testJaatha <- function(dm, n.points=2, reps=1, seed=12523, cores=detect_cores(), 
                       folder=".", grid.pars='all', test_data=NULL, ...) {
  
  # Setup parallization backend
  mc.opt <- list(preschedule = FALSE, set.seed = FALSE)
  registerDoMC(cores[1])
  
  # Set up directories
  folder.results <- file.path(folder, "results")
  folder.logs  <- file.path(folder, "logs")
  if ( file.exists(folder.results) || file.exists(folder.logs) ) {
    stop("Folders already exists")
  }
  dir.create(folder.results, recursive = TRUE)
  dir.create(folder.logs, recursive = TRUE)
  if ( !(file.exists(folder.results) && file.exists(folder.logs)) ) {
    stop("Failed to create folders")
  }
  
  # Create test data sets if not provided
  if (is.null(test_data)) {
    set.seed(seed + 1)
    test_data <- createTestData(dm, n.points, reps, grid.pars, cores = prod(cores))
    save(test_data, file = file.path(folder.results, 'test_datasets.Rda'))
  }
  
  faulty <- vapply(test_data$data, inherits, logical(1), what = "simpleError")
  if (any(faulty)) {
    print(test_data$data[faulty])
    stop("Error in test data")
  }
  
  # Sample seeds for each run
  set.seed(seed)
  seeds <- sample.int(2 ^ 20, nrow(test_data$par_grid))
  
  # Write some information about this run to disk
  envir <- c(folder = folder,
             jaatha.version = as.character(packageVersion("jaatha")),
             test_jaatha.version = as.character(packageVersion("testJaatha")),
             hostname = Sys.info()["nodename"],
             seed = seed)
  print(envir)
  write.table(envir, file = paste0(folder.logs, "/envir.txt"), col.names = FALSE)
  
  n <- nrow(test_data$par_grid)
  # The actual simulation
  i <- 0
  results <- foreach(i = 1:n, .combine = rbind, .options.multicore = mc.opt) %dopar% { 
    cat("Run", i, "of", n, "\n")
    log <- file(file.path(folder.logs, paste0("run_", i, ".txt")))
    sink(log)
    sink(log, type = "message")
    seed <- seeds[i]

    cat("----------------------------------------------------------------------\n")
    cat("Run", i, "of" , n, "\n")
    cat("Real parameters:", test_data$par_grid[i,], "\n")
    cat("Seed:", seed, "\n")
    cat("----------------------------------------------------------------------\n")
    res <- try_capture_stack({
      jaatha_model <- create_jaatha_model(dm)
      jaatha_data <- create_jaatha_data(test_data$data[[i]], jaatha_model)
      reps <- 2
      save(jaatha_model, jaatha_data, seed, reps,
           file = file.path(folder.logs, paste0("run_", i, "_setup.Rda")))
      
      set.seed(seed)
      runtimes <- system.time(
        result <- jaatha(jaatha_model, jaatha_data, 
                         repetitions = 2, cores = cores[2], ...)
      )
      
      save(result, runtimes,
           file = file.path(folder.logs, paste0("run_", i, "_result.Rda")))
      estimates <- result$param
      sink(NULL)
      sink(NULL)
      c(runtimes, estimates)
    }, new.env())
    
    if (inherits(res, "simpleError")) {
      cat("Error:", res$message, "\n")
      print(res$call)
      sink(NULL)
      sink(NULL)
      cat("Run", i, " Error:", res$message, "\n")
      return(rep(NA, ncol(test_data$par_grid) + 6))
    }
    
    res
  }
  
  estimates <- results[, -(1:5)]
  runtimes <- results[, 1:5]
  
  write.table(estimates, file = paste0(folder.results, "/estimates.txt"), row.names = FALSE)
  write.table(test_data$par_grid,  file = paste0(folder.results, "/true_values.txt"), row.names = FALSE)
  write.table(runtimes,  file = paste0(folder.results, "/runtimes.txt"), row.names = FALSE)
}


#' Creates test datasets
#' 
#' @param dm The demographic model to test
#' @param n.points The number of test points for each parameter.
#' @param reps The number of repetitions of each parameter combination
#' @param grid.pars The index of parameter that are used for builing the grid
#'   of true values. Use 'all' (default) for all parameters.
#' @param cores The number of cores on which the simulations are distributed.
#' @param grid.values The values that are used for the parameters. A named list.
#' 
#' @export
#' @examples
#' library('jaatha')
#' dm <- coala:::model_theta_tau()
#' test_data <- createTestData(dm, 2, 2)
createTestData <- function(dm, n.points=2, reps=1, grid.pars='all', 
                           cores=prod(detect_cores()),
                           grid.values = NULL) {
  
  test_data <- list()
  
  # Create the parameter grid for true values
  test_data$par_grid <- createParGrid(dm, n.points, reps, grid.pars, grid.values)
  
  # Simulate test data sets
  seeds <- sample(10000000, nrow(test_data$par_grid))
  test_data$data <- mclapply(1:nrow(test_data$par_grid), function(x) {
    try_capture_stack({
      set.seed(seeds[x])
      simulate(dm, pars = test_data$par_grid[x,])
    }, new.env())
  }, mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE)

  test_data
}

#' Create a parameter grid
#' 
#' @importFrom coala get_parameter_table
#' @inheritParams createTestData
#' @examples
#' dm <- coala:::model_theta_tau()
#' testJaatha:::createParGrid(dm, 2, 1)
#' testJaatha:::createParGrid(dm, 2, 1, grid.pars=1)
#' testJaatha:::createParGrid(dm, 2, 1, grid.pars=2)
#' testJaatha:::createParGrid(dm, 2, 1, grid.pars=1:2)
#' testJaatha:::createParGrid(dm, 2, 1, grid.pars=1, 
#'                            grid.values = data.frame(tau = 1:3))
createParGrid <- function(dm, n.points, reps, grid.pars='all', grid.values = NULL){
  par.ranges <- get_parameter_table(dm)[,2:3]
  rownames(par.ranges)  <- get_parameter_table(dm)[,1]
  n.dim <- nrow(par.ranges)
  
  if (any(grid.pars != 'all')) {
    grid.par.mask = 1:n.dim %in% grid.pars
  } else {
    grid.par.mask = rep(TRUE, n.dim)
  }
  
  if (!is.null(grid.values)) {
    par.values <- grid.values
  } else {
    par.values <- data.frame(apply(par.ranges[grid.par.mask, , drop = FALSE], 1, function(x) {
      seq(x[1], x[2], length = n.points + 2)[-c(1, n.points + 2)]
    }))
  }

  # Create the grid
  par.grid <- expand.grid(par.values)
  
  # Add middle values for the non-grid parameters
  if (any(grid.pars != 'all')) {
    non.grid.pars <- t(as.matrix(apply(par.ranges[!grid.par.mask,], 1, function(x) {
      mean(as.numeric(x[1:2]))
    })))
    par.grid = cbind(par.grid, non.grid.pars)
    par.grid = par.grid[ , rownames(par.ranges)] # Sort again
  }
  
  # And repeat
  apply(par.grid, 2, rep, reps)
}

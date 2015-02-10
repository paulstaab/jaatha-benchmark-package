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
#' @param Additional argument for \code{Jaatha.initialize}.
#' @inheritParams createTestData
#' 
#' @return Nothing.
#' @export
#' @examples
#' library('jaatha')
#' dm <- coalsimr:::model_theta_tau()
#' testJaatha(dm, 2, 1, cores=c(2,1), folder=tempfile(), scaling.factor=10)
#' 
#' test_data <- createTestData(dm, 2, 1)
#' testJaatha(dm, test_data = test_data, cores=c(1,1), folder=tempfile())
testJaatha <- function(dm, n.points=2, reps=1, seed=12523, cores=c(16,2), 
                       folder=".", grid.pars='all', test_data=NULL, ...) {
  
  # Setup parallization backend
  mc.opt <- list(preschedule=FALSE, set.seed=FALSE)
  registerDoMC(cores[1])
  
  # Set up directories
  folder.results <- paste(folder, "results", sep="/")
  folder.logs  <- paste(folder, "logs", sep="/")
  if ( file.exists(folder.results) || file.exists(folder.logs) ) {
    stop("Folders already exists")
  }
  dir.create(folder.results, recursive=T)
  dir.create(folder.logs, recursive=T)
  
  # Create test data sets if not provided
  if (is.null(test_data)) {
    set.seed(seed+1)
    test_data <- createTestData(dm, n.points, reps, grid.pars, cores = prod(cores))
    save(test_data, file = paste(folder.results, 'test_datasets.Rda' , sep="/"))
  }
  
  # Sample seeds for each run
  set.seed(seed)
  seeds <- sample.int(2^20, nrow(test_data$par_grid))
  
  # Write some information about this run to disk
  envir <- c(folder=folder,
             jaatha.version=as.character(packageVersion("jaatha")),
             test_jaatha.version=as.character(packageVersion("testJaatha")),
             hostname=Sys.info()["nodename"],
             seed=seed)
  print(envir)
  write.table(envir, file=paste0(folder.logs, "/envir.txt"), col.names=F)
  
  n <- nrow(test_data$par_grid)
  # The actual simulation
  results <- foreach(i=1:n, .combine=rbind, .options.multicore=mc.opt) %dopar% { 
                       
    cat("Run",i,"of",n,"\n")
    log <- file(paste(folder.logs, "/run_", i, ".txt", sep=""))
    sink(log)
    sink(log, type = "message")
    set.seed(seeds[i])
    cat("----------------------------------------------------------------------\n")
    cat("Run",i,"of",n,"\n")
    cat("Real parameters:",test_data$par_grid[i,],"\n")
    cat("----------------------------------------------------------------------\n")
    tryCatch({
      jaatha <- Jaatha.initialize(data = test_data$data[[i]], 
                                  model = dm,
                                  cores = cores[2],
                                  ...)
      save(jaatha, file=paste(folder.logs, "/run_", i, ".Rda", sep=""))
      
      runtimes <- rep(0, 6)
      names(runtimes) <-
        c('init.user','init.system','init.elapsed','ref.user','ref.system','ref.elapsed')
      
      runtimes[1:3] <- system.time(
        jaatha <- Jaatha.initialSearch(jaatha)
      )
      save(jaatha, file=paste(folder.logs, "/run_", i, ".Rda", sep=""))
      
      runtimes[4:6] <- system.time(
        jaatha <- Jaatha.refinedSearch(jaatha, 2)
      )
      save(jaatha, file=paste(folder.logs, "/run_", i, ".Rda", sep=""))
      estimates <-  Jaatha.getLikelihoods(jaatha)[1,-(1:2)]
      sink(NULL)
      sink(NULL)
      return(c(runtimes, estimates))
    }, error = function(e) {
      cat("Error:", e$message,
          file=paste(folder.logs, "/run_", i, ".txt", sep=""), append=TRUE)
      sink(NULL)
      sink(NULL)
      cat("Run", i, " ERROR:", e$msg, "\n")
      print(e)
      return(NA)
    })
  }
  
  estimates <- results[, -(1:6)]
  runtimes <- results[, 1:6]
  
  write.table(estimates, file=paste(folder.results, "estimates.txt", sep="/"), row.names=F)
  write.table(test_data$par_grid,  file=paste(folder.results, "true_values.txt", sep="/"), row.names=F)
  write.table(runtimes,  file=paste(folder.results, "runtimes.txt", sep="/"), row.names=F)
}

#' Creates test datasets
#' 
#' @param dm The demographic model to test
#' @param n.points The number of test points for each parameter.
#' @param reps The number of repetitions of each parameter combination
#' @param grid.pars The index of parameter that are used for builing the grid
#'   of true values. Use 'all' (default) for all parameters.
#' @param cores The number of cores on which the simulations are distributed.
#' 
#' @importFrom jaatha dm.simSumStats
#' 
#' @export
#' @examples
#' library('jaatha')
#' dm <- coalsimr:::model_theta_tau()
#' test_data <- createTestData(dm, 2, 2)
createTestData <- function(dm, n.points=2, reps=1, grid.pars='all', cores=2) {
  test_data <- list()
  
  dm <- dm + sumstat_seg_sites()
  
  # Create the parameter grid for true values
  test_data$par_grid <- createParGrid(dm, n.points, reps, grid.pars)
  
  # Simulate test data sets
  seeds <- sample(10000000, nrow(test_data$par_grid))
  test_data$data <- mclapply(1:nrow(test_data$par_grid), function(x) {
    set.seed(seeds[x])
    simulate(dm, pars = test_data$par_grid[x,])
  }, mc.cores=cores, mc.preschedule=FALSE, mc.set.seed=FALSE)

  test_data
}


createParGrid <- function(dm, n.points, reps, grid.pars='all'){
  par.ranges <- get_parameter_table(dm)
  n.dim <- nrow(par.ranges)
  
  if (any(grid.pars != 'all')) {
    grid.pars.mask = 1:n.dim %in%grid.pars
    par.ranges = par.ranges[grid.pars.mask,]
  }
  
  par.values <- data.frame(apply(par.ranges, 1, function(x) {
    seq(as.numeric(x[2]), as.numeric(x[3]),
               length=n.points+2)[-c(1,n.points+2)]
  }))
  
  par.grid <- expand.grid(par.values)
  
  if (any(grid.pars != 'all')) {
    non.grid.pars <- t(as.matrix(apply(get_parameter_table(dm)[!grid.pars.mask,], 1, function(x) {
      mean(as.numeric(x[2:3]))
    })))
    par.grid = cbind(par.grid, non.grid.pars)
  }
  
  par.grid <- apply(par.grid, 2, rep, reps)
  colnames(par.grid) <- get_parameter_table(dm)[,1]
  par.grid
}

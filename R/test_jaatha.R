#' Excutes Jaathas self test
#' @export
#' @examples
#' library('jaatha')
#' dm <- dm.createThetaTauModel(10:11, 100)
#' testJaatha(dm, 2, 1, cores=c(1,2), folder=tempfile())
#' 
#' test_data <- createTestData(dm, 2, 1)
#' testJaatha(dm, test_data = test_data, cores=c(1,2), folder=tempfile())
testJaatha <- function(dm, n.points=2, reps=1, seed=12523, smoothing=FALSE, 
                       cores=c(16,2), folder=".", grid.pars='all', 
                       test_data=NULL, scaling_factor=1) {
  
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
                                  smoothing = smoothing,
                                  scaling.factor =  scaling_factor)
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
#' @export
#' @examples
#' library('jaatha')
#' dm <- dm.createThetaTauModel(10:11, 100)
#' dm <- dm.addSummaryStatistic(dm, 'fpc')
#' test_data <- createTestData(dm, 2, 2)
createTestData <- function(dm, n.points=2, reps=1, grid.pars='all', cores=2) {
  test_data <- list()
  
  # Set Summary Statistics to JSFS + seg.sites
  dm@sum.stats <- dm.createDemographicModel(5:6, 100)@sum.stats
  dm <- dm.addSummaryStatistic(dm, 'seg.sites')
  
  # Create the parameter grid for true values
  test_data$par_grid <- createParGrid(dm, n.points, reps, grid.pars)
  
  # Simulate test data sets
  seeds <- sample(10000000, nrow(test_data$par_grid))
  test_data$data <- mclapply(1:nrow(test_data$par_grid), function(x) {
    set.seed(seeds[x])
    dm.simSumStats(dm.addSummaryStatistic(dm, 'seg.sites'), test_data$par_grid[x,])
  }, mc.cores=cores, mc.preschedule=FALSE, mc.set.seed=FALSE)

  test_data
}


createParGrid <- function(dm, n.points, reps, grid.pars='all'){
  par.ranges <- jaatha:::dm.getParRanges(dm)
  n.dim <- length(jaatha:::dm.getParameters(dm))
  
  if (any(grid.pars != 'all')) {
    grid.pars.mask = 1:n.dim %in%grid.pars
    par.ranges = par.ranges[grid.pars.mask,]
  }
  
  par.values <- data.frame(apply(par.ranges, 1, function(x) seq(x[1], x[2],length=n.points+2)[-c(1,n.points+2)]))
  par.grid <- expand.grid(par.values)
  
  if (any(grid.pars != 'all')) {
    non.grid.pars <- t(as.matrix(apply(jaatha:::dm.getParRanges(dm)[!grid.pars.mask,], 1, mean)))
    par.grid = cbind(par.grid, non.grid.pars)
    par.grid = par.grid[,row.names(jaatha:::dm.getParRanges(dm))]
  }
  
  par.grid <- apply(par.grid, 2, rep, reps)
  return(par.grid)
}

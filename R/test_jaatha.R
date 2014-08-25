#' @export
testJaatha <- function(dm, n.points=2, reps=1, seed=12523, smoothing=FALSE, 
                       cores=c(16,2), folder=".", fpc=FALSE, grid.pars='all') {
  
  # Set up directories
  folder.results <- paste(folder, "results", sep="/")
  folder.logs  <- paste(folder, "logs", sep="/")
  if ( file.exists(folder.results) || file.exists(folder.logs) ) {
    stop("Folders already exists")
  }
  dir.create(folder.results, recursive=T)
  dir.create(folder.logs, recursive=T)
  
  # Write some information about this run to disk
  envir <- c(folder=folder,
             jaatha.version=as.character(packageVersion("jaatha")),
             test_jaatha.version=as.character(packageVersion("testJaatha")),
             hostname=Sys.info()["nodename"],
             seed=seed)
  print(envir)
  write.table(envir, file=paste0(folder.logs, "/envir.txt"), col.names=F)
  
  # Create the parameter grid for true values
  par.grid <- createParGrid(dm, n.points, reps, grid.pars)
  n <- dim(par.grid)[1]
  set.seed(seed)
  seeds <- sample.int(2^20, n)
  
  # Setup parallization backend
  registerDoMC(cores[1])
  mc.opt <- list(preschedule=FALSE, set.seed=FALSE)
  
  # The actual simulation
  results <- foreach(i=1:n, .combine=rbind, .options.multicore=mc.opt) %dopar% { 
    cat("Run",i,"of",n,"\n")
    log <- file(paste(folder.logs, "/run_", i, ".txt", sep=""))
    sink(log)
    sink(log, type = "message")
    set.seed(seeds[i])
    cat("----------------------------------------------------------------------\n")
    cat("Run",i,"of",n,"\n")
    cat("Real parameters:",par.grid[i,],"\n")
    cat("----------------------------------------------------------------------\n")
    tryCatch({
      if (fpc) {
        jsfs <- dm.simSumStats(jaatha:::dm.addSummaryStatistic(dm, 'seg.sites'), par.grid[i, ])
        jaatha <- Jaatha.initialize(dm, jsfs=jsfs, 
                                    cores=cores[2],
                                    smoothing=smoothing)
      } else {
        jsfs <- dm.simSumStats(dm, par.grid[i, ])
        jaatha <- Jaatha.initialize(dm, jsfs=jsfs, 
                                    cores=cores[2],
                                    smoothing=smoothing)
      }
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
  write.table(par.grid,  file=paste(folder.results, "true_values.txt", sep="/"), row.names=F)
  write.table(runtimes,  file=paste(folder.results, "runtimes.txt", sep="/"), row.names=F)
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

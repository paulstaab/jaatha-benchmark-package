#' @export
testJaatha <- function(dm, n.points=2, reps=1, seed=12523, smoothing=FALSE, 
                       cores=c(8,4), folder=".", fpc=FALSE) {
  # Set up directories
  folder.results <- paste(folder, "results", sep="/")
  folder.logs  <- paste(folder, "logs", sep="/")
  if ( file.exists(folder.results) || file.exists(folder.logs) ) {
    stop("Folders already exists")
  }
  dir.create(folder.results, recursive=T)
  dir.create(folder.logs, recursive=T)
  
  envir <- c(jaatha.version=as.character(packageVersion("jaatha")),
             test_jaatha.version=as.character(packageVersion("testJaatha")),
             hostname=Sys.info()["nodename"],
             seed=seed)
  print(envir)
  write.table(envir, file=paste0(folder.logs, "/envir.txt"), col.names=F)

  set.seed(seed)
  par.grid <- createParGrid(dm, n.points, reps)
  n <- dim(par.grid)[1]
  seeds <- sample.int(2^20, n)

  registerDoMC(cores[1])

  results <- foreach(i=1:n, .combine=rbind) %dopar% { 
      cat("Run",i,"of",n,"\n")
      sink(file=paste(folder.logs, "/run_", i, ".txt", sep=""))
      set.seed(seeds[i])
      cat("----------------------------------------------------------------------\n")
      cat("Run",i,"of",n,"\n")
      cat("Real parameters:",par.grid[i,],"\n")
      cat("----------------------------------------------------------------------\n")
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

      runtimes <- rep(0, 6)
      names(runtimes) <-
        c('init.user','init.system','init.elapsed','ref.user','ref.system','ref.elapsed')

      runtimes[1:3] <- system.time(
        jaatha <- Jaatha.initialSearch(jaatha)
      )

      runtimes[4:6] <- system.time(
        jaatha <- Jaatha.refinedSearch(jaatha, 2)
      )
      estimates <-  Jaatha.getLikelihoods(jaatha)[1,-(1:2)]
      sink(NULL)
      return(c(runtimes, estimates))
  }

  estimates <- results[, -(1:6)]
  runtimes <- results[, 1:6]

  write.table(estimates, file=paste(folder.results, "estimates.txt", sep="/"), row.names=F)
  write.table(par.grid,  file=paste(folder.results, "true_values.txt", sep="/"), row.names=F)
  write.table(runtimes,  file=paste(folder.results, "runtimes.txt", sep="/"), row.names=F)
}

createParGrid <- function(dm, n.points, reps){
  par.ranges <- jaatha:::dm.getParRanges(dm)
  n.dim <- length(jaatha:::dm.getParameters(dm))
  n.runs <- n.points^n.dim
  par.grid <- matrix(0,n.runs,n.dim)
  
  for (i in 1:n.dim){
    pars.current.dim <- seq(par.ranges[i,1],par.ranges[i,2],length=n.points+2)[-c(1,n.points+2)]
    par.grid[,i] <- rep(pars.current.dim,each=n.points^(i-1))
  }

  par.grid <- apply(par.grid, 2, rep, reps)
  return(par.grid)
}

dm.getParRanges <- function(dm,inklExtTheta=T){
    parMask <- !is.na(dm@features$parameter)
    parRanges <- cbind(lower=dm@features$lowerRange[parMask],upper=dm@features$upperRange[parMask])
    parRanges <- parRanges[sort.list(dm@features$parameter[parMask]),]
    if (!inklExtTheta) parRanges <- parRanges[1:dm.getNPar(dm),]
    return( parRanges )
}

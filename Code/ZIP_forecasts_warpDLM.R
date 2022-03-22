#install.packages(c("doParallel", "foreach", "KFAS", "truncdist", "doSNOW", "tscount", "VGAM", "mvtnorm", "MASS", "TruncatedNormal", "optimParallel"))
#install.packages("remotes")
#remotes::install_github("drkowal/rSTAR")
library(parallel)
library(doParallel)
library(rSTAR)
library(foreach)
library(KFAS)
library(truncdist)
library(doSNOW)
library(VGAM)
library(mvtnorm)
library(TruncatedNormal)
library(MASS)
library(optimParallel)

source("./Code/helper_functions.R")

numCores <- detectCores()

#Run simulations
comb <- function(...) {
  mapply('rbind', ..., SIMPLIFY=FALSE)
}
ptm <- proc.time()
#cl <- makeCluster(numCores)
cl <- makeCluster(numCores, outfile="") # if you want to see everything on the terminal
registerDoSNOW(cl)
on.exit(stopCluster(cl))
numSeries <- 30
iterations <- 50
itSeq <- seq(from=100, by=2, length.out=iterations)
pb <- txtProgressBar(max = numSeries*iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
zipFC <- foreach(j=1:numSeries) %:%
  foreach(i=itSeq, .combine=comb, .packages = c("KFAS", "truncdist", "VGAM", "mvtnorm", "TruncatedNormal", "MASS"),
          .options.snow = opts, .multicombine=TRUE, .errorhandling = "pass") %dopar% {
            #Simulate Data
            set.seed(32*j)
            tslen=200
            lambda <- array(0, tslen)
            lambda[1] <- runif(1,5,15)
            for(k in 2:tslen){
              lambda[k] <- lambda[k-1] + rnorm(1, mean=0, sd=.2)
            }
            y <- pmin(rzipois(tslen,lambda, runif(1, .1,.3)),rep(24,tslen))
            dim(y) <- NULL
            #print(length(unique(y))/length(y))
            #print(table(y)[1])
            #plot(y)
            
            #Run DLM-STAR simulations for forecasts
            init_mod <- KFAS::SSModel(y[1:i] ~ SSMtrend(1, Q=15, a1=0, P1 =3),  H = 15)
            ymax=24
            
            #Square root
            sqrt_gibbs <- DLM_STAR_mcmc_KFAS(y[1:i], init_mod = init_mod, 
                                             transformation = "sqrt", y_max=ymax,
                                             update_mod = updatemod_unif_loclev_KFAS,
                                             nburn=500, nsave=1500, nskip=0)
            a_y = a_j(y[1:i], ymax)
            a_yp1 = a_j(y[1:i] + 1, ymax)
            lower = sqrt_gibbs$g_func(a_y)
            upper= sqrt_gibbs$g_func(a_yp1)
            tempFunc <- function(par){
              locLevLike_log(3, par[1], par[2], lower, upper)
            }
            sqrt_fc <- tryCatch(expr={
              optimResults <- optim(par=c(mean(sqrt_gibbs$V_post), mean(sqrt_gibbs$W_post)), fn=tempFunc, 
                                    lower=c(0, 0), upper=rep(Inf, 2), method="L-BFGS-B",control=list(factr=100))
              dirSampler_locLev(y[1:i], V=optimResults$par[1], W=optimResults$par[2], g=sqrt_gibbs$g_func,
                                g_inv=sqrt_gibbs$g_inv_func, y_max=ymax)
            },
            error=function(cond){
              #message(cond)
              print(paste("Error with sqrt in series", j, "and time point", i))
              return(array(-999,10000))
            })
            
            
            #Identity
            id_gibbs <- DLM_STAR_mcmc_KFAS(y[1:i], init_mod = init_mod, 
                                           transformation = "identity", y_max=ymax,
                                           update_mod = updatemod_unif_loclev_KFAS,
                                           nburn=500, nsave=1500, nskip=0)
            a_y = a_j(y[1:i], ymax)
            a_yp1 = a_j(y[1:i] + 1, ymax)
            lower = id_gibbs$g_func(a_y)
            upper=id_gibbs$g_func(a_yp1)
            tempFunc <- function(par){
              locLevLike_log(3, par[1], par[2], lower, upper)
            }
            id_fc <- tryCatch(expr={
              optimResults <- optim(par=c(mean(id_gibbs$V_post),mean(id_gibbs$W_post)), fn=tempFunc, 
                                    lower=c(0, 0), upper=rep(Inf, 2), method="L-BFGS-B",control=list(factr=100))
              dirSampler_locLev(y[1:i], V=optimResults$par[1], W=optimResults$par[2], g=id_gibbs$g_func,
                                g_inv=id_gibbs$g_inv_func, y_max=ymax)
            },
            error=function(cond){
              #message(cond)
              print(paste("Error with id in series", j, "and time point", i))
              return(array(-999,10000))
            })
            
            
            #Empirical CDF
            np_gibbs <- DLM_STAR_mcmc_KFAS(y[1:i], init_mod = init_mod, 
                                           transformation = "np", y_max=ymax,
                                           update_mod = updatemod_unif_loclev_KFAS,
                                           nburn=500, nsave=1500, nskip=0)
            a_y = a_j(y[1:i], ymax)
            a_yp1 = a_j(y[1:i] + 1, ymax)
            lower = np_gibbs$g_func(a_y)
            upper=np_gibbs$g_func(a_yp1)
            tempFunc <- function(par){
              locLevLike_log(3, par[1], par[2], lower, upper)
            }
            np_fc <- tryCatch(expr={
              optimResults <- optim(par=c(mean(np_gibbs$V_post),mean(np_gibbs$W_post)), fn=tempFunc, 
                                    lower=c(0, 0), upper=rep(Inf, 2), method="L-BFGS-B",control=list(factr=100))
              dirSampler_locLev(y[1:i], V=optimResults$par[1], W=optimResults$par[2], g=np_gibbs$g_func,
                                g_inv=np_gibbs$g_inv_func, y_max=ymax)
            },
            error=function(cond){
              #message(cond)
              print(paste("Error with np in series", j, "and time point", i))
              return(array(-999,10000))
            })
            
            star_results <- list(as.vector(sqrt_fc), 
                                 as.vector(id_fc), 
                                 as.vector(np_fc))
            
            return(c(star_results, list(y)))
          }

proc.time() - ptm
saveRDS(zipFC, file='./Outputs/ModelResults/simulations/ZIP_forecasts_warpDLM_direct.rds')

zipFC <- readRDS(file="./Outputs/ModelResults/simulations/ZIP_forecasts_warpDLM_direct.rds")

errorSearch <- function(fc_dat){
  nSeries <- length(fc_dat)
  nMethods <- length(fc_dat[[1]])-1
  nTimepoints <- nrow(fc_dat[[1]][[1]])
  matrix <- matrix(NA,0,3)
  for(i in 1:nSeries){
    for(j in 1:nMethods){
      for(k in 1:nTimepoints){
        if(fc_dat[[i]][[j]][k,1]==-999){
          matrix = rbind(matrix, c(i,j,k))
        }
      }
    }
  }
  print(paste("There are",nrow(matrix),"errors"))
  return(matrix)
}
errors = errorSearch(zipFC)

done = (nrow(errors)==0)
while(!done){
  for(m in 1:nrow(errors)){
    cur = errors[m,]
    #Get correct data 
    y <- zipFC[[cur[1]]][[4]][1,]
    i = itSeq[cur[3]]
    #Run DLM-STAR simulations for forecasts
    init_mod <- KFAS::SSModel(y[1:i] ~ SSMtrend(1, Q=15, a1=0, P1 =3),  H = 15)
    ymax=24
    if(cur[2]==1){
      #Square root
      sqrt_gibbs <- DLM_STAR_mcmc_KFAS(y[1:i], init_mod = init_mod, 
                                         transformation = "sqrt", y_max=ymax,
                                         update_mod = updatemod_unif_loclev_KFAS,
                                         nburn=500, nsave=1500, nskip=0)
      a_y = a_j(y[1:i], ymax)
      a_yp1 = a_j(y[1:i] + 1, ymax)
      lower = sqrt_gibbs$g_func(a_y)
      upper= sqrt_gibbs$g_func(a_yp1)
      tempFunc <- function(par){
        locLevLike_log(3, par[1], par[2], lower, upper)
      }
      newfc <- tryCatch(expr={
        optimResults <- optim(par=c(mean(sqrt_gibbs$V_post), mean(sqrt_gibbs$W_post)), fn=tempFunc, 
                              lower=c(0, 0), upper=rep(Inf, 2), method="L-BFGS-B",control=list(factr=100))
        dirSampler_locLev(y[1:i], V=optimResults$par[1], W=optimResults$par[2], g=sqrt_gibbs$g_func,
                          g_inv=sqrt_gibbs$g_inv_func, y_max=ymax)
      },
      error=function(cond){
        #message(cond)
        print(paste("Error with sqrt in series", j, "and time point", i))
        return(array(-999,10000))
      })
    }
    else if(cur[2]==2){
      #Identity
      id_gibbs <- DLM_STAR_mcmc_KFAS(y[1:i], init_mod = init_mod, 
                                     transformation = "identity", y_max=ymax,
                                     update_mod = updatemod_unif_loclev_KFAS,
                                     nburn=500, nsave=1500, nskip=0)
      a_y = a_j(y[1:i], ymax)
      a_yp1 = a_j(y[1:i] + 1, ymax)
      lower = id_gibbs$g_func(a_y)
      upper=id_gibbs$g_func(a_yp1)
      tempFunc <- function(par){
        locLevLike_log(3, par[1], par[2], lower, upper)
      }
      newfc <- tryCatch(expr={
        optimResults <- optim(par=c(mean(id_gibbs$V_post),mean(id_gibbs$W_post)), fn=tempFunc, 
                              lower=c(0, 0), upper=rep(Inf, 2), method="L-BFGS-B",control=list(factr=100))
        dirSampler_locLev(y[1:i], V=optimResults$par[1], W=optimResults$par[2], g=id_gibbs$g_func,
                          g_inv=id_gibbs$g_inv_func, y_max=ymax)
      },
      error=function(cond){
        #message(cond)
        print(paste("Error with id in series", j, "and time point", i))
        return(array(-999,10000))
      })
    }
    else if(cur[2]==3){
      #Empirical CDF
      np_gibbs <- DLM_STAR_mcmc_KFAS(y[1:i], init_mod = init_mod, 
                                     transformation = "np", y_max=ymax,
                                     update_mod = updatemod_unif_loclev_KFAS,
                                     nburn=500, nsave=1500, nskip=0)
      a_y = a_j(y[1:i], ymax)
      a_yp1 = a_j(y[1:i] + 1, ymax)
      lower = np_gibbs$g_func(a_y)
      upper=np_gibbs$g_func(a_yp1)
      tempFunc <- function(par){
        locLevLike_log(3, par[1], par[2], lower, upper)
      }
      newfc <- tryCatch(expr={
        optimResults <- optim(par=c(mean(np_gibbs$V_post),mean(np_gibbs$W_post)), fn=tempFunc, 
                              lower=c(0, 0), upper=rep(Inf, 2), method="L-BFGS-B",control=list(factr=100))
        dirSampler_locLev(y[1:i], V=optimResults$par[1], W=optimResults$par[2], g=np_gibbs$g_func,
                          g_inv=np_gibbs$g_inv_func, y_max=ymax)
      },
      error=function(cond){
        #message(cond)
        print(paste("Error with np in series", j, "and time point", i))
        return(array(-999,10000))
      })
    }
    else{
      stop("You shouldn't have reached here")
    }
    zipFC[[cur[1]]][[cur[2]]][cur[3],]=newfc
  }
  print("Loop complete")
  errors = errorSearch(zipFC)
  done = (nrow(errors)==0)
}

saveRDS(zipFC, file='./Outputs/ModelResults/simulations/ZIP_forecasts_warpDLM.rds')

#install.packages(c("doParallel", "foreach", "KFAS", "truncdist", "doSNOW", "tscount"))
#install.packages("remotes")
#remotes::install_github("drkowal/rSTAR")
library(parallel)
library(doParallel)
library(rSTAR)
library(foreach)
library(KFAS)
library(truncdist)
library(doSNOW)
library(tscount)

source("./Code/helper_functions.R")

(numCores <- detectCores())

#Run simulations
comb <- function(...) {
  mapply('rbind', ..., SIMPLIFY=FALSE)
}
ptm <- proc.time()
cl <- makeCluster(numCores)
#cl <- makeCluster(numCores, outfile="") # if you want to see everything on the terminal
registerDoSNOW(cl)
numSeries <- 30
iterations <- 50
itSeq <- seq(from=100, by=2, length.out=iterations)
pb <- txtProgressBar(max = numSeries*iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
inGarchFC <- foreach(j=1:numSeries) %:%
  foreach(i=itSeq, .combine=comb, .packages = c("KFAS", "truncdist", "tscount"),
          .options.snow = opts, .multicombine=TRUE) %dopar% {
            #Simulate Data
            set.seed(32*j)
            y <- tsglm.sim(200, param=list(intercept=.3, past_obs=.6, past_mean=.2),
                           model = list(past_obs=1, past_mean=1))$ts
            #print(length(unique(y))/length(y))
            #print(table(y)[1])
            #plot(y)
            
            #Run DGLM models
            update_model <- function(pars, model) {
              model["Q"] <- pars[1]
              model
            }
            #check that variances are non-negative
            check_model <- function(model) {
              (model["Q"] > 0)
            }
            pois_mod <- SSModel(y[1:i] ~ -1 + SSMtrend(degree=1, Q=matrix(NA)), distribution = "poisson")
            pois_fit <- fitSSM(pois_mod, inits = c(1), method = "BFGS",
                               updatefn = update_model, checkfn = check_model)
            pois_fc <- importSim(pois_fit$model, nsim= 1250, n.ahead=1)
            
            #Negative Binomial 
            nb_mod <- SSModel(y[1:i] ~ -1 + SSMtrend(degree=1, Q=matrix(NA)), 
                              distribution = "negative binomial")
            nb_fit <- fitSSM(nb_mod, inits = c(0.001), method = "BFGS",
                             updatefn = update_model, checkfn = check_model)
            nb_fc <- importSim(nb_fit$model, nsim= 1250, n.ahead=1)
            dglm_results <- list(as.vector(pois_fc), as.vector(nb_fc))
            
            return(c(dglm_results, list(y)))
            
          }
close(pb)
stopCluster(cl)
proc.time() - ptm
saveRDS(inGarchFC, file='./Outputs/ModelResults/simulations/INGARCH_forecasts_dglm.rds')
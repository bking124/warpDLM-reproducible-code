library(tidyverse)
library(ggplot2)
library(wesanderson)
library(ddst)

#Analysis flag
#Currently showing INGARCH (flag "ing")
#Change to "zip" to show zero-inflated Poisson results
flag = "ing"
#flag = "zip"

#Read in data
if(flag == "ing"){
  dglm_dat <- readRDS('Outputs/ModelResults/simulations/INGARCH_forecasts_dglm.rds')
  warpDLM_dat <- readRDS('Outputs/ModelResults/simulations/INGARCH_forecasts_warpDLM.rds')
}

if (flag == "zip"){
  dglm_dat <- readRDS('Outputs/ModelResults/simulations/ZIP_forecasts_dglm.rds')
  warpDLM_dat <- readRDS('Outputs/ModelResults/simulations/ZIP_forecasts_warpDLM.rds')
}

listLen <- length(warpDLM_dat[[1]])
warpDLM_dat <- lapply(warpDLM_dat, function(x){x[[listLen]] <- NULL; return(x)})

full_dat <- mapply(c, warpDLM_dat, dglm_dat, SIMPLIFY = FALSE)
listLen <- length(full_dat[[1]])

#Simulation stores many versions of y; get rid of them
full_dat <- lapply(full_dat, function(x){x[[listLen]] <- x[[listLen]][1,]; return(x)})

score_titles <- c("logarithmic", "quadratic", "spherical", 
                  "rankprob", "dawseb", "normsq", "sqerror")
methods_titles <- c("sqrt-wDLM", "id-wDLM", "np-wDLM", "PoisDGLM", "nbDGLM")
y_all <- sapply(full_dat, function(x){x[[listLen]]})

#Function that takes results from simulation for single time series
#and turns it into a list of scores for each time point and method
#Pass in list and vector of points where data was passed in to the methods 
scoreCalc <- function(sim_results, timePoints){
  y <- sim_results[[listLen]]
  len <- length(sim_results)
  n <- length(timePoints)
  y_test <- y[timePoints+1]
  
  #Calculate scores for STAR-DLM and for Poisson DGLM
  means <- matrix(NA, nrow=n, ncol=(len-1))
  stdevs <- matrix(NA, nrow=n, ncol=(len-1))
  probs_list <- list()
  cdf_list <- list()
  nbins <- max(unlist(sim_results[-len]))
  for(i in 1:(len-1)){
    probs_mat <- matrix(NA, nrow=n, ncol=nbins)
    dat <- sim_results[[i]]
    means[,i] <- rowMeans(dat)
    stdevs[,i] <- apply(dat,1, sd)
    probs_mat <- apply(dat+1, 1, tabulate, nbins=nbins)/ncol(dat)
    cdf_mat <- apply(probs_mat, 2, cumsum)
    cdf_list <- c(cdf_list, list(cdf_mat))
    probs_list <- c(probs_list, list(probs_mat))
  }
  
  #Based on means and standard deviations
  sqerror <- (y_test-means)^2
  normsq <- (sqerror/(stdevs^2))
  dawseb <- normsq + 2*log(stdevs)
  colMeans(sqerror)
  colMeans(normsq)
  colMeans(dawseb)
  
  #Based on mass function
  calc_prob <- function(prob_mat){
    prob <- array(0, dim=n)
    for(i in 1:n){
      prob[i] <- prob_mat[y_test[i]+1,i]
      if(prob[i]==0)
        prob[i]=0.0001
    }
    return(prob)
  }
  pt_yt <- sapply(probs_list, calc_prob)
  log_scores <- -log(pt_yt)
  colMeans(log_scores)
  
  pnorm <- sapply(probs_list, function(x){sqrt(colSums(x^2))})
  quadratic <- pnorm^2 - 2*pt_yt
  colMeans(quadratic)
  
  spherical <- -pt_yt/pnorm
  colMeans(spherical)
  
  #Based on CDF
  calc_rank <- function(cdf_mat){
    rank_score <- array(0, dim=n)
    for(i in 1:n){
      indicator_vec <- (y_test[i] <= 0:(nbins-1))
      rank_score[i] <- sum((cdf_mat[,i]-indicator_vec)^2)
    }
    return(rank_score)
  }
  rank_scores <- sapply(cdf_list, calc_rank)
  colMeans(rank_scores)
  
  #Get mean scores across time points
  dlm_score_list <- list(log_scores, quadratic, spherical, 
                         rank_scores, dawseb, normsq, sqerror)
  score_list = dlm_score_list
  for(index in 1:length(score_list)){
    colnames(score_list[[index]]) <- methods_titles[1:(len-1)]
  }
  names(score_list) <- score_titles
  
  #Change normsq so that min is better
  score_list$normsq <- (score_list$normsq-1)^2
  return(score_list)
}

numSeries <- 30
iterations <- 50
itSeq <- seq(from=100, by=2, length.out=iterations)

full_score_list <- full_dat
for(i in 1:numSeries){
  full_score_list[[i]] <- scoreCalc(full_dat[[i]], itSeq)
}

#Get mean scores for each TS
meanScores_byTS <- lapply(full_score_list, function(x){sapply(x, colMeans)})

##############################Plot log score chart####################################
log_scores <- t(sapply(meanScores_byTS, function(x){x[,1]}))
log_scores_diff <- ((log_scores-log_scores[,4])/log_scores[,4])*100
log_scores_diff <- log_scores_diff[,-4]
df <- pivot_longer(data.frame(log_scores_diff), cols=everything(), names_to="Method")
df$Method <- factor(df$Method , levels=c("nbDGLM","sqrt.wDLM", "id.wDLM", "np.wDLM"))

# Violin chart
png(filename = paste0("Outputs/Figures/logscore_",flag,".png"), width=800)
ggplot(df, aes(x=Method, y=value, fill=Method)) + geom_violin() + 
  geom_hline(yintercept = 0) + ylab("% Log Score Difference vs. Poisson DGLM")+
  scale_fill_manual(values=wes_palette(n=5, name="Darjeeling1")[-1])
dev.off()

##############################Plot calibration chart####################################
#Compute p values for test of uniformity
ddst <- matrix(NA, 30, length(methods_titles))
for(j in 1:30){
  for(k in 1:length(methods_titles)){
    df <- full_dat[[j]][[k]]
    y <- full_dat[[j]][[length(full_dat[[1]])]]
    iterations <- 50
    timePoints <- seq(from=100, by=2, length.out=iterations)
    y_test <- y[timePoints+1]
    pvals = array(NA, 50)
    for(i in 1:50){
      empcdf = ecdf(df[i,])
      pvals[i] = runif(1,empcdf(y_test[i]-1),empcdf(y_test[i]))
    }
    ddst[j,k] <- ddst.uniform.test(pvals, compute.p = T)$p.value
  }
}
colnames(ddst) <- methods_titles
ddst_df <- pivot_longer(data.frame(ddst), cols=everything(), names_to="Method")
ddst_df$Method <- factor(ddst_df$Method , levels=c("PoisDGLM", "nbDGLM","sqrt.wDLM", "id.wDLM", "np.wDLM"))

#Box plot of p values
png(filename = paste0("Outputs/Figures/cal_pval_",flag,".png"), width=800)
ggplot(ddst_df, aes(x=Method, y=value, fill=Method)) + geom_boxplot() +
  geom_hline(yintercept = 0.05) + ylab("P-value")+
  scale_fill_manual(values=wes_palette(n=5, name="Darjeeling1"))
dev.off()
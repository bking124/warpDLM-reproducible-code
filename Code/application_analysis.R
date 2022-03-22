library(dlm)
library(rSTAR)
library(tidyverse)
library(mc2d)
library(bayesplot)
library(TruncatedNormal)
library(spatstat)  #ewcdf function
library(mvnfast)

source("./Code/helper_functions.R")

#Load Cincy Data
load("./Data/heroinCallsCounts")
load("./Data/drugCallsCounts")

df <- data.frame(heroin = heroinCounts$Count, other=otherCounts$Count, Date=heroinCounts$Date)
df <- subset(df, Date >= '2018-01-01')
sub_df <- subset(df, Date < '2020-01-01' & Date >= '2018-01-01')
df_online <- subset(df, Date >='2020-01-01')

dat <- cbind(sub_df$other, sub_df$heroin)
dat_online <- cbind(df_online$other, df_online$heroin)

######################################Offline Analysis###################################
#Can skip this section and instead load model results from save

#Initial model
init_mod <- dlm(FF = matrix(c(1, 0), nrow = 1) %x% diag(2), V = diag(2),
                GG = matrix(c(1, 0, 1, 1), 2, 2) %x% diag(2),
                W = bdiag(diag(2), diag(2)),
                m0 = c(dat[1, 1], dat[1, 2], 0, 0),
                C0 = diag(x = 3, nrow = 4))

#This can take several hours to run
drug_results_np <- SUTSE_mcmc_dlm(dat, init_mod, update_mod = updatemod_invWish_dlm, transformation = "np",
                                  nsave=10000, nburn=5000, nskip=1, nfc=1, particle = T)


save(drug_results_np, file="./Outputs/ModelResults/application/CincyData_offline")

######################################Particle Filtering###################################
#This section also can be skipped and model results can simply be loaded from saved files

load(file="./Outputs/ModelResults/application/CincyData_offline")      #drug_results_np

#Get np g
y <- dat
g_wrap <- apply(y, 2, g_cdf, distribution = "np")
g_np <- function(y){ylist <- split(y, col(y)); mapply(function(x, y) x(y), g_wrap, ylist)}
t_grid <- function(y, y_max){
  sort(unique(c(0:min(2*max(y), y_max),
                seq(0, min(2*max(y), y_max), length.out = 100),
                quantile(unique(y[y!=0]), seq(0, 1, length.out = 100)))))
}
t_grid_list = apply(y,2,t_grid, y_max=Inf)
g_inv_wrap <- mapply(g_inv_approx, g_wrap, t_grid_list)
g_inv_np = function(y){ylist <- split(y, col(y)); mapply(function(x, y) x(y), g_inv_wrap, ylist)}

pf_np_all <- SUTSE_PF(init_mod = init_mod, offline_mod = drug_results_np, 
                       dat_online = dat_online, g=g_np, g_inv=g_inv_np, resample="all")

#Don't need to save draws of z_pred, theta_pred, or theta_samp for analysis
pf_np_all[[3]] <- NULL
pf_np_all[[1]] <- NULL
pf_np_all[[4]] <- NULL

save(pf_np_all, file="./Outputs/ModelResults/application/CincyData_online", compress="xz")

######################################Explore PF Results###################################
load(file="./Outputs/ModelResults/application/CincyData_offline")      #drug_results_np
load(file="./Outputs/ModelResults/application/CincyData_online")       #pf_np_all

y_pred_samples <- pf_np_all$y_pred
weights <- pf_np_all$weights
Neff <-  pf_np_all$Neff
tpl <- pf_np_all$update_time
particles <- pf_np_all$theta_particles

#Get np g
y <- dat
g_wrap <- apply(y, 2, g_cdf, distribution = "np")
g_np <- function(y){ylist <- split(y, col(y)); mapply(function(x, y) x(y), g_wrap, ylist)}
t_grid <- function(y, y_max){
  sort(unique(c(0:min(2*max(y), y_max),
                seq(0, min(2*max(y), y_max), length.out = 100),
                quantile(unique(y[y!=0]), seq(0, 1, length.out = 100)))))
}
t_grid_list = apply(y,2,t_grid, y_max=Inf)
g_inv_wrap <- mapply(g_inv_approx, g_wrap, t_grid_list)
g_inv_np = function(y){ylist <- split(y, col(y)); mapply(function(x, y) x(y), g_inv_wrap, ylist)}

#Draw from 'smoothing' distribution
FF <- matrix(c(1, 0), nrow = 1) %x% diag(2)
offline_mod <- drug_results_np
#Get V and W
W <- apply(offline_mod$W_post, c(2,3), mcmcMean, sd=FALSE)
V <- apply(offline_mod$V_post, c(2,3), mcmcMean, sd=FALSE)

drawSamp <- function(x){
  z <- rmvn(1,FF%*%x, V)
  y <- round_fun(g_inv_np(z))
  return(y)
}

y_post_pf <- apply(particles, c(1,2), drawSamp)
y_smooth_pf <- apply(y_post_pf, c(1,2), median)
y_smooth_offline <- apply(offline_mod$post_pred, c(2,3), median)
y_smooth_all <- cbind(t(y_smooth_offline), y_smooth_pf)

#Produce smoothing plot
df_withsmooth <- data.frame(df, heroin_smooth = y_smooth_all[2,], 
                            other_smooth=y_smooth_all[1,])
long_df <- df_withsmooth %>%  gather(key = "variable", value = "Count", -Date)

png(filename = "Outputs/Figures/ODCountsTS.png", width=1200)
ggplot(long_df, aes(x = Date, y = Count)) + 
  geom_line(aes(color = variable,  alpha=variable), size=1) +
  scale_color_manual(values = c("#00AFBB", "#00AFBB", "#E7B800","#E7B800")) +
  scale_alpha_manual(values = c(.5, 1, .5, 1))+
  legend_none()
dev.off()


#Plot efficiency
#hist(Neff, breaks = "Scott", )
png(filename = "Outputs/Figures/PF_ESS.png", width = 800)
plot(1:length(Neff), Neff, type='l', xlab="Particle Filter Time Step", ylab="ESS")
dev.off()

#Plot Time in Loop
png(filename = "Outputs/Figures/PF_TPL.png", width=800)
plot(1:length(tpl), tpl, xlab="Particle Filter Time Step", ylab="Seconds in Loop Iteration")
dev.off()

#Plot PIT
numTest = length(df_online$other)-1
uniform <- replicate(100,runif(numTest))
y_test <- df_online$other[-1]
y_test2 <- df_online$heroin[-1]
pvals = array(NA, numTest)
pvals2 = array(NA, numTest)
for(i in 1:numTest){
  empcdf = ewcdf(y_pred_samples[i,,1], weights = weights[i,])
  empcdf2 = ewcdf(y_pred_samples[i,,2],  weights = weights[i,])
  pvals[i] = runif(1,empcdf(y_test[i]-1),empcdf(y_test[i]))
  pvals2[i] = runif(1,empcdf2(y_test2[i]-1),empcdf2(y_test2[i]))
}
x <- seq(0,1, length.out=numTest)
png(filename = "Outputs/Figures/calibration.png", width = 1200)
par(pty="s",mfrow=c(1,2))
plot(punif(x), sort(pvals), xlab="rPIT Values", ylab="Uniform Quantiles",
     main="Other OD Calibration", xlim=c(0,1), ylim=c(0,1))
matlines(x, apply(uniform,2,sort), col=alpha(rgb(0,0,0), 0.05), lty=1)
lines(x,x, lty=2)
plot(punif(x), sort(pvals2), xlab="rPIT Values", ylab="Uniform Quantiles",
     main="Heroin OD Calibration", xlim=c(0,1), ylim=c(0,1))
matlines(x, apply(uniform,2,sort), col=alpha(rgb(0,0,0), 0.05), lty=1)
lines(x,x, lty=2)
dev.off()

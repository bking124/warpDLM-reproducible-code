#Works with CDF transformations but right now all series must be the same transformation
SUTSE_mcmc_dlm <- function (y, init_mod, 
                            update_mod=function(mod_level, z_star, theta){return(mod_level)}, 
                            transformation = "log", lambda = NULL, y_max = Inf, 
                            nsave = 5000, nburn = 5000, nskip = 2, nfc=5, particle=FALSE) {
  #Checking for proper input
  if (any(y < 0) || any(y != floor(y))) 
    stop("y must be nonnegative counts")
  if (any(y > y_max)) 
    stop("y must not exceed y_max")
  transformation = tolower(transformation)
  if (!is.element(transformation, c("identity", "log", "sqrt", 
                                    "np", "pois", "neg-bin", "box-cox"))) 
    stop("The transformation must be one of 'identity', 'log', 'sqrt', 'np', 'pois', 'neg-bin', or 'box-cox'")
  transform_family = ifelse(test = is.element(transformation, 
                                              c("identity", "log", "sqrt", "box-cox")), yes = "bc", 
                            no = "cdf")
  n = length(y)
  if (transform_family == "bc") {
    if (transformation == "identity") 
      lambda = 1
    if (transformation == "log") 
      lambda = 0
    if (transformation == "sqrt") 
      lambda = 1/2
    if (transformation == "box-cox") 
      lambda = runif(n = 1)
    g = function(t) g_bc(t, lambda = lambda)
    g_inv = function(s) g_inv_bc(s, lambda = lambda)
    sum_log_deriv = (lambda - 1) * sum(log(y + 1))
  }
  if (transform_family == "cdf") {
    g_wrap <- apply(y, 2, g_cdf, distribution = transformation)
    g <- function(y){ylist <- split(y, col(y)); mapply(function(x, y) x(y), g_wrap, ylist)}
    t_grid <- function(y, y_max){
      sort(unique(c(0:min(2*max(y), y_max),
                    seq(0, min(2*max(y), y_max), length.out = 100),
                    quantile(unique(y[y!=0]), seq(0, 1, length.out = 100)))))
    }
    t_grid_list = apply(y,2,t_grid, y_max=y_max)
    g_inv_wrap <- mapply(g_inv_approx, g_wrap, t_grid_list)
    g_inv = function(y){ylist <- split(y, col(y)); mapply(function(x, y) x(y), g_inv_wrap, ylist)}
    lambda = NULL
  }
  
  n = nrow(y) #Length of time series
  m=ncol(y) #Number of time series
  z_star = g(y + abs(rnorm(n = n)))
  mod_level <- init_mod
  filt <- dlmFilter(z_star, mod_level)
  theta <- dlmBSample(filt)
  p <- ncol(mod_level$W) #Length of state vector
  
  a_y = a_j(y)
  a_yp1 = a_j(y + 1)
  
  nstot = nburn + (nskip + 1) * (nsave)
  skipcount = 0
  isave = 0
  
  #Allocate space to save
  V_save <- array(0, dim = c(nsave, m, m))
  W_save <- array(0, dim = c(nsave, p, p))
  y_fc_save <- array(0, dim = c(nsave, nfc, m))
  pred_save <- array(0, dim = c(nsave, n, m))
  theta_save <- array(0, dim=c(nsave, p))
  
  #MCMC Loop
  ptm <- proc.time()
  pb <- txtProgressBar(0, nstot, style = 3)
  progressBar=interactive()
  for (it in 1 : nstot){
    if (progressBar) 
      setTxtProgressBar(pb, it)
    
    # Draw zstar
    z_star <- y
    muvec <- t(tcrossprod(filt$mod$FF, as.matrix(theta[-1,])))
    mulist <- split(muvec, row(muvec))
    lower <- g(a_y)
    lowerlist <- split(lower, row(lower))
    upper <- g(a_yp1)
    upperlist <- split(upper, row(upper))
    z_star <- t(mapply(function(x,y,z){rtmvnorm(1,mu=x,sigma=filt$mod$V, lb=y, ub=z)},
                       mulist, lowerlist, upperlist))
    
    ## draw the states: FFBS
    filt <- dlmFilter(z_star, mod_level, simplify=T)
    theta <- dlmBSample(filt)

    ## Update model
    mod_level <- update_mod(mod_level, z_star, theta)
    
    #Draw forecast and transform to y scale
    y_fc <- dlmForecast(filt, nAhead=nfc, sampleNew=1)$newObs %>% unlist %>% matrix(.,nfc,m) %>%
      g_inv %>% round_fun
    
    # update and save
    if(it > nburn){
      skipcount = skipcount + 1
      if (skipcount > nskip) {
        isave = isave + 1
        V_save[isave,,] <- mod_level$V
        W_save[isave,,] <- mod_level$W
        theta_save[isave, ] <- theta[nrow(theta),]
        y_fc_save[isave, ,] <- y_fc
        #Draw form posterior predictive
        pred_save[isave,,] <- rmultinormal(n=n, 
                                           mean = t(tcrossprod(filt$mod$FF, 
                                                               as.matrix(theta[-1,]))), 
                                           sigma = as.vector(filt$mod$V)) %>% 
          g_inv %>% round_fun
        skipcount=0
      }
    }
  }
  end <- (proc.time()-ptm)[1]
  print(paste("Time taken: ", round(end, digits=5), " seconds"))
  
  #Return results
  if(particle)
    return(list(V_post = V_save, W_post =W_save, fc_post = y_fc_save, 
                post_pred = pred_save, theta_last = theta_save))
  else
    return(list(V_post = V_save, W_post =W_save, fc_post = y_fc_save, post_pred = pred_save))
  
}

#Update model for multivariate
updatemod_invWish_dlm <-
  function(mod_level, z_star, theta){
    TT= nrow(z_star)
    
    # prior hyperparameters
    delta0 <- delta2 <- 3 ; delta1 <- 100
    V0 <- (delta0 - 2) *diag(c(2^2, 3^2))
    Wmu0 <- (delta1 - 2) * diag(0.01^2, 2)
    Wbeta0 <- (delta2 - 2) * diag(c(1, 2^2))
    
    #V0 <- (delta0 - 2) *diag(c(10^2, 500^2))
    #Wmu0 <- (delta1 - 2) * diag(0.01^2, 2)
    #Wbeta0 <- (delta2 - 2) * diag(c(5^2, 100^2))
    
    # update V
    S <- crossprod(z_star - theta[-1, ] %*% t(mod_level$FF)) + V0
    mod_level$V <- solve(dlm::rwishart(df = delta0 + 1 + TT, p = 2, Sigma = solve(S)))
    
    # update Wmu and Wbeta
    theta.center <- theta[-1, ] - (theta[-(TT + 1), ] %*% t(mod_level$GG))
    SS1 <- crossprod(theta.center)[1 : 2, 1 : 2] + Wmu0
    SS2 <- crossprod(theta.center)[3 : 4, 3 : 4] + Wbeta0
    gibbsWmu <- solve(dlm::rwishart(df =delta1 + 1 + TT, Sigma = solve(SS1)))
    gibbsWbeta <- solve(dlm::rwishart(df = delta2 + 1 + TT, Sigma = solve(SS2)))
    mod_level$W <- bdiag(gibbsWmu, gibbsWbeta)
    
    return (mod_level)
  }

#Particle Filter for multivariate model
SUTSE_PF <- function(init_mod, offline_mod, dat_online, g, g_inv, S=10000, resample = "threshold", g_dyn = -1){
  if(!(resample %in% c("all","threshold")))
    stop("Resampling algorithm must be 'all' or 'threshold'")
  
  #Get V and W
  W <- apply(offline_mod$W_post, c(2,3), mcmcMean, sd=FALSE)
  V <- apply(offline_mod$V_post, c(2,3), mcmcMean, sd=FALSE)
  
  #Outer loop is for number of time points
  theta_samp <- offline_mod$theta_last
  p=ncol(theta_samp)
  n <- ncol(dat_online)
  Tmax = nrow(dat_online)
  Fmat <- init_mod$FF
  Gmat <- init_mod$GG
  
  numResample = 0
  prevWeights = rep(1,S)
  weights = array(0,S)
  Neff_vec = array(0,Tmax)
  timeInLoop = array(0,Tmax)
  theta_samples <- array(NA, c(Tmax, S, p))
  theta_particles <- array(NA, c(Tmax, S, p))
  theta_pred_samples <- array(NA, c(Tmax, S, p))
  z_pred_samples <- array(NA, c(Tmax, S, n))
  y_pred_samples <- array(NA, c(Tmax, S, n))
  weights_samples <- array(NA, c(Tmax, S)) 
  
  #Calculate things that don't depend on time
  Sig_z <-  V + Fmat%*%W%*%t(Fmat)
  Sig_ztheta <- Fmat%*%W
  Sig_z_inv <- solve(Sig_z)
  V_1var <- W - t(Sig_ztheta)%*%Sig_z_inv%*%Sig_ztheta
  for(i in 1:Tmax){
    ptm=proc.time()
    
    #Reestimate g
    if(g_dyn!=-1 & i%%g_dyn==0){
      y <- dat_online[(i-g_dyn+1):i,]
      g_wrap <- apply(y, 2, g_cdf, distribution = 'np')
      g <- function(y){ylist <- split(y, col(y)); mapply(function(x, y) x(y), g_wrap, ylist)}
      t_grid_mat = apply(y,2,t_grid, y_max=Inf)
      if(is.matrix(t_grid_mat))
        t_grid_mat = split(t_grid_mat, col(t_grid_mat))
      g_inv_wrap <- mapply(g_inv_approx, g_wrap, t_grid_mat)
      g_inv = function(y){ylist <- split(y, col(y)); mapply(function(x, y) x(y), g_inv_wrap, ylist)}
    }
    
    #Sample from importance density
    mu_z <- apply(theta_samp, 1, function(x) {Fmat%*%Gmat%*%x})
    mu_theta <- apply(theta_samp, 1, function(x) {Gmat%*%x})
    a_y = a_j(dat_online[i,])
    a_yp1 = a_j(dat_online[i,]+1)
    lower <- g(matrix(a_y, 1,2))
    upper <- g(matrix(a_yp1, 1,2))

    theta_propose = matrix(NA, nrow=S, ncol=p)
    for(j in 1:S){
      #V_1samp <- mvrnorm(1, rep(0,nrow(V_1var)), V_1var)
      V_1samp <- as.vector(rmvn(1, rep(0,nrow(V_1var)), V_1var))
      V_0samp <- TruncatedNormal::rtmvnorm(1, rep(0,nrow(Sig_z)),Sig_z,lower-mu_z[,j], upper-mu_z[,j])
      theta_propose[j,] <- mu_theta[,j] + V_1samp + 
        t(Sig_ztheta)%*%Sig_z_inv%*%as.matrix(V_0samp)
      weights[j] = mvtnorm::pmvnorm(mean = mu_z[,j], sigma=Sig_z, 
                                    lower=as.vector(lower), upper=as.vector(upper))
    }
    weights <- prevWeights*weights
    normWeights = weights/sum(weights)
    Neff = 1/sum(normWeights^2)
    
    #Resample
    if(resample=="threshold"){
      if(Neff < (S/2)){
        index = sample(S, S, replace=TRUE, prob=normWeights)
        theta_samp <- theta_propose[index, ]
        numResample = numResample+1
        prevWeights <- rep(1,S)
      }
      else{
        theta_samp <- theta_propose
        prevWeights <- weights
      }
    }
    else if(resample=="all"){
      index = sample(S, S, replace=TRUE, prob=normWeights)
      theta_samp <- theta_propose[index, ]
      numResample = numResample+1
      prevWeights <- rep(1,S)
    }
    else
      stop("You shouldn't have gotten to this point")

    
    #theta_pred <- apply(theta_samp,1,function(x) mvrnorm(1, Gmat%*%x, W))
    theta_pred <- apply(theta_samp,1,function(x) rmvn(1, Gmat%*%x, W))
    #z_pred <- t(apply(theta_pred,2,function(x) mvrnorm(1, Fmat%*%x, V)))
    z_pred <- t(apply(theta_pred,2,function(x) rmvn(1, Fmat%*%x, V)))
    y_pred <- round_fun(g_inv(z_pred))
    print(paste("Done with time point",i))
    #Store results
    timeInLoop[i] = (proc.time()-ptm)[3]
    theta_samples[i,,] <- theta_samp
    theta_particles[i,,] <- theta_propose
    theta_particles <- theta_particles[,index,]
    theta_pred_samples[i,,] <-  theta_pred
    z_pred_samples[i,,] <-  z_pred
    y_pred_samples[i,,] <- y_pred
    weights_samples[i,] <- normWeights
    Neff_vec[i] = Neff
  }
  
  #return(list(theta_samp=theta_samples, y_pred=y_pred_samples, weights=weights_samples,
  #            Neff=Neff_vec,update_time=timeInLoop, numResample=numResample))
  return(list(theta_samp=theta_samples,theta_particles=theta_particles, theta_pred=theta_pred_samples, weights=weights_samples,
              z_pred=z_pred_samples,y_pred=y_pred_samples,Neff=Neff_vec, 
              update_time=timeInLoop, numResample=numResample))
}
#Rounding Function
a_j <- function (j, y_max = Inf) 
{
  val = j
  val[j <= 0] = -Inf
  val[j >= y_max + 1] = Inf
  return(val)
}

#Transform back to y scale
round_fun <- function (z, y_max = Inf) 
{
  pmin(floor(z) * I(z > 0), y_max)
}


#CDF Estimator function (taken directly from rSTAR)
g_cdf <- function (y, distribution = "np") 
{
  distribution = tolower(distribution)
  if (!is.element(distribution, c("np", "pois", "neg-bin", 
                                  "box-cox"))) 
    stop("The distribution must be one of 'np', 'pois', or 'neg-bin'")
  n = length(y)
  mu_y = mean(y)
  sigma_y = sd(y)
  if (distribution == "np") {
    F_y = function(t) n/(n + 1) * ecdf(y)(t)
  }
  if (distribution == "pois") {
    F_y = function(t) ppois(t, lambda = mu_y)
  }
  if (distribution == "neg-bin") {
    if (mu_y >= sigma_y^2) {
      warning("'neg-bin' not recommended for underdispersed data")
      sigma_y = 1.1 * sqrt(abs(mu_y))
    }
    F_y = function(t) pnbinom(t, size = mu_y^2/(sigma_y^2 - 
                                                  mu_y), prob = mu_y/sigma_y^2)
  }
  t0 = sort(unique(y[y != 0]))
  g0 = mu_y + sigma_y * qnorm(F_y(t0 - 1))
  t0 = t0[which(is.finite(g0))]
  g0 = g0[which(is.finite(g0))]
  splinefun(t0, g0, method = "monoH.FC")
}

#Inverted g for nonparametric CDF transformation
g_inv_approx <- function (g, t_grid) 
{
  g_grid = g(t_grid)
  function(s) {
    sapply(s, function(si) t_grid[which.min(abs(si - g_grid))])
  }
}

#Box-Cox Transformation
g_bc <- function (t, lambda) 
{
  if (lambda == 0) {
    sign(t) * log(abs(t))
  }
  else {
    (sign(t) * abs(t)^lambda - 1)/lambda
  }
}

#Inverse Box-Cox Transformation
g_inv_bc <- function (s, lambda) 
{
  if (lambda == 0) {
    exp(s)
  }
  else {
    sign(lambda * s + 1) * abs(lambda * s + 1)^(1/lambda)
  }
}

#Likelihood for univariate local level model
#Assumes theta0 mean is 0
locLevLike <- function(P0, V, W, lower, upper){
  Tlen = length(lower)
  Vmat <- diag(V, Tlen)
  j=P0+W
  Sigma_theta <- matrix(0, Tlen, Tlen)
  for(i in 1:Tlen){
    Sigma_theta[i,] <- rep(j, Tlen)
    j=j+W
  }
  s <- diag(Sigma_theta)
  Sigma_theta[lower.tri(Sigma_theta)] <- 0
  Sigma_theta <- Sigma_theta + t(Sigma_theta) - diag(s)
  Sigma_z <- Sigma_theta + Vmat
  mu =rep(0,nrow(Sigma_z))
  return(mvtnorm::pmvnorm(as.numeric(lower), as.numeric(upper), sigma=Sigma_z))
}

#Returns negative log likelihood for marginal maximum likelihood estimation
locLevLike_log <- function(P0, V, W, lower, upper, logParams=FALSE){
  if(logParams){
    V=exp(V)
    W=exp(W)
  }
  Tlen = length(lower)
  Vmat <- diag(V, Tlen)
  j=P0+W
  Sigma_theta <- matrix(0, Tlen, Tlen)
  for(i in 1:Tlen){
    Sigma_theta[i,] <- rep(j, Tlen)
    j=j+W
  }
  s <- diag(Sigma_theta)
  Sigma_theta[lower.tri(Sigma_theta)] <- 0
  Sigma_theta <- Sigma_theta + t(Sigma_theta) - diag(s)
  Sigma_z <- Sigma_theta + Vmat
  mu =rep(0,nrow(Sigma_z))
  return(-log(mvtnorm::pmvnorm(as.numeric(lower), as.numeric(upper), sigma=Sigma_z)[1]))
}

#Returns Highest Density Region
getHDR <- function(probs, vals, level){
  if(sum(probs)< .95)
    print("Probabilities sum to less than 95%: consider expanding range")
  if(sum(probs)< level)
    stop("Probabilities given to do not sum to desired level")
  if(length(probs)!= length(vals))
    stop("Values vector must have same length as probabilities vector")
  temp <- data.frame(probs,vals)
  ordered <- temp[order(temp$probs, decreasing=T),]
  index <- which(probs==max(probs))
  sum=0
  count=1
  while(sum < level){
    sum=sum+ordered[count,1]
    count=count+1
  }
  return(ordered[1:(count-1),])
}


#MCMC for univariate models
#Utilizes KFAS
DLM_STAR_mcmc_KFAS <- 
  function (y, init_mod, update_mod=function(fit, z_star, theta){return(fit)}, 
            transformation = "log", lambda = NULL, y_max = Inf, 
            nsave = 5000, nburn = 5000, nskip = 2, nfc=1) {
    #Checking for proper input
    if(!is.SSModel(init_mod))
      stop("initial model is not proper KFAS model")
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
      g = g_cdf(y = y, distribution = transformation)
      t_grid = sort(unique(round(c(seq(0, min(2 * max(y), 
                                              y_max), length.out = 250), quantile(unique(y[y < 
                                                                                             y_max] + 1), seq(0, 1, length.out = 250))), 8)))
      g_inv = g_inv_approx(g = g, t_grid = t_grid)
      lambda = NULL
    }
    
    
    n = length(y)
    z_star = g(y + abs(rnorm(n = n)))
    
    #first use of KFAS package
    fit <- init_mod
    fit["y"] <- z_star
    theta <- drop(simulateSSM(fit, type="states", nsim=1))
    p <- ncol(unlist(fit$Q))
    
    a_y = a_j(y, y_max)
    a_yp1 = a_j(y + 1, y_max)
    
    nstot = nburn + (nskip + 1) * (nsave)
    skipcount = 0
    isave = 0
    
    #Allocate space to save
    H_save <- numeric(nsave)
    Q_save <- matrix(nrow=nsave, ncol=p)
    lambda_save <- numeric(nsave)
    zstar_save <- matrix(nrow=nsave, ncol=n)
    theta_save <- matrix(nrow=nsave, ncol=n)
    y_fc_save <- matrix(nrow=nsave, ncol=nfc)
    pred_save <- matrix(nrow=nsave, ncol=n)
    
    #MCMC Loop
    ptm <- proc.time()
    # pb <- txtProgressBar(0, nstot, style = 3)
    #progressBar=interactive()
    for (it in 1 : nstot){
      #if (progressBar) 
      #  setTxtProgressBar(pb, it)
      
      # Draw zstar
      z_star = rSTAR::rtruncnormRcpp(y_lower = g(a_y), 
                                     y_upper = g(a_yp1), 
                                     mu = drop(tcrossprod(drop(fit$Z), as.matrix(theta))), 
                                     sigma = rep(sqrt(drop(fit$H)), n), 
                                     u_rand = runif(n = n))
      
      fit["y"] <- z_star
      
      ## draw the states: FFBS
      theta <- drop(KFAS::simulateSSM(fit, type="states", nsim=1))

      ## Update model
      fit <- update_mod(fit, z_star, theta)
      
      #Draw forecast and transform to y scale
      #Right now is mostly fit for the local level model
      filt <- KFAS::KFS(fit)
      temp <- rnorm(1,mean=filt$a[n+1], sd=sqrt(filt$P[n+1]+fit$H))
      y_fc <- round_fun(g_inv(temp), y_max)
      
      # update and save
      if(it > nburn){
        skipcount = skipcount + 1
        if (skipcount > nskip) {
          isave = isave + 1
          H_save[isave] <- drop(fit$H)
          Q_save[isave,] <- drop(fit$Q)
          theta_save[isave,] <- theta
          y_fc_save[isave, ] <- y_fc
          #Draw from posterior predictive
          temp <- rnorm(n=n, mean = tcrossprod(drop(fit$Z), as.matrix(theta)), 
                        sd = rep(sqrt(drop(fit$H)), n))
          pred_save[isave,] <- round_fun(g_inv(temp), y_max) 
          skipcount=0
        }
      }
    }
    end <- (proc.time()-ptm)[1]
    print(paste("Time taken: ", round(end, digits=5), " seconds"))
    list(V_post = H_save, W_post =Q_save, fc_post = y_fc_save, 
           post_pred = pred_save, theta_samp = theta_save, g_func = g, g_inv_func=g_inv)
  }

#Update mod for KFAS
#Currently does not handle intercept
updatemod_unif_loclev_KFAS <-
  function(fit, z_star, theta){
    
    A <- 10^4
    n <- length(z_star)
    if(is.null(ncol(theta)))
      p <-1
    else
      p <- ncol(theta)
    theta <- as.matrix(theta, n, p)
    
    shape.y <- (n-1)/2
    shape.theta <- shape.y
    rate.y <- crossprod(z_star - theta[,1]) / 2
    
    eta <- truncdist::rtrunc(n = 1, 
                             'gamma',   # Family of distribution
                             a = 1/A^2, # Lower interval
                             b = Inf,   # Upper interval
                             shape = shape.y,
                             rate =  rate.y)
    fit["H"] <- 1/eta
    
    ## draw system precision W
    theta.center <- theta[-1, , drop = FALSE] - 
      tcrossprod(theta[-n, , drop = FALSE], drop(fit$T))
    rate.theta <- colSums(theta.center^2)/2
    if(rate.theta==0)
      rate.theta=1e-30
    fit["Q"] <- diag(1/truncdist::rtrunc(n = p, 
                              'gamma',   # Family of distribution
                              a = 1/A^2, # Lower interval
                              b = Inf,   # Upper interval
                              shape = shape.theta,
                              rate =  rate.theta), p)
    return(fit)
  }

#Wrapper function for importanceSSM to draw for Poisson and NegBin models
#Only works for the simple random walk type model
#Only works in univariate case
#Returns 4*nsim draws from predictive distribution
importSim <- function(object, nsim, n.ahead){
  #This means there is only one state
  j=1
  #Update model with appropriate NAs
  timespan <- attr(object, "n") + 1:n.ahead
  n <- attr(object, "n") <- attr(object, "n") + as.integer(n.ahead)
  endtime<-end(object$y) + c(0, n.ahead)
  object$y <- window(object$y, end = endtime, extend = TRUE)
  object$u <- rbind(object$u, matrix(object$u[1, ], nrow = n.ahead,
                                     ncol = ncol(object$u), byrow = TRUE))
  
  #Draw samples from signal (alpha_t)
  imp <- KFAS::importanceSSM(object, "signal",nsim = nsim, antithetics = TRUE)
  nsim <- as.integer(4 * nsim)
  w <- imp$weights/sum(imp$weights)
  imp$samples <- imp$samples[timespan, , , drop = F]
  imp$samples[,j, ] <- exp(imp$samples[, j, ])
  n <- as.integer(length(timespan))
  preds <- sapply(1:n, function(i) {
    sample_mu <- sample(imp$samples[i, j, ], size = nsim, replace = TRUE,
                        prob = w)
    if(object$distribution == "poisson")
      result= rpois(n = nsim, lambda = sample_mu)
    else
      result= rnbinom(n = nsim, size = object$u[timespan[i],j], mu = sample_mu)
    return(result)
  })
  return(preds)
}

#Pass in vector x, return plot of ergodic average
ergAverage <- function(x){
  ave <- cumsum(x)/(1:length(x))
  plot(ave, type='l')
  title(main="Ergodic Average")
}

#Direct sampler for local level model using selection normal theory
dirSampler_locLev <- function(y, V, W, g, g_inv, y_max, R0=3){
  n = length(y)
  
  #Create system matrices
  mu_z <- array(0,n)
  #Create Sigma_z
  Sigma_z <- matrix(0, n,n)
  j=W+R0
  for(i in 1:n){
    Sigma_z[i,] <- rep(j, n)
    j=j+W
  }
  Sigma_z[lower.tri(Sigma_z)] <- 0
  Sigma_z <- Sigma_z + t(Sigma_z)
  diag(Sigma_z) <- seq(W+R0, n*W+R0, by=W) +V
  
  #Create Sigma_ztheta
  Sigma_ztheta <- Sigma_z
  diag(Sigma_ztheta) <- diag(Sigma_ztheta)-V
  Sigma_theta <- Sigma_ztheta
  
  #Try Sampling
  #Draw from distribution
  V_1var <- Sigma_theta - t(Sigma_ztheta)%*%solve(Sigma_z)%*%Sigma_ztheta
  V_1samp <- mvrnorm(10^4, rep(0,nrow(V_1var)), V_1var)
  #V_1samp <- rmvn(10^4, rep(0,nrow(V_1var)), V_1var) #rmvn is faster but Cholesky decomposition is more unstable
  V_0samp <- rtmvnorm(10^4, rep(0,nrow(Sigma_z)),Sigma_z,lb=g(a_y), ub=g(a_yp1))
  
  theta_samp <- t(V_1samp) + t(Sigma_ztheta)%*%solve(Sigma_z)%*%t(V_0samp)
  theta_samp <- t(theta_samp)
  y_pred <- round_fun(g_inv(rnorm(10^4, mean=theta_samp[,ncol(theta_samp)], sd=sqrt(V+W))), y_max=y_max)
  
  return(y_pred)
}

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
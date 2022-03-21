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
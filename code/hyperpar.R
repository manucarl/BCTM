hyperpar <- function(Z, K1, K2, A, c = 0.1, alpha = 0.1,
  omegaseq, omegaprob, R = 1000, myseed = 123,
  thetaseq = NULL, type = "SD", truncate = 0.95)
{
  # Z = rows from the tensor product design matrix
  # K1, K2 = precision matrices
  # A = constraint matrix
  # c = threshold from eq. (8) in Klein & Kneib (2016)
  # alpha = probability parameter from eq. (8) in Klein & Kneib (2016)
  # omegaseq = sequence of weights for the anisotropy
  # omegaprob = prior probabilities for the weights
  # R = number of simulations
  # myseed = seed to be used internally for optimisation
  # thetaseq = sequence of hyperparameter values to consider.
  #        If NULL, numerical optimisation will be performed (time consuming!!)
  # type = type of hyperprior for tau/tau^2
  #        options: IG => IG(1,theta) for tau^2
  #                 SD => WE(0.5,theta) for tau^2
  #                 HN => HN(0,theta) for tau
  # note: any variability arising from the parameters in the null space is
  #       ignored!

  require("mvtnorm")
  require("MASS")

  I1 <- diag(ncol(K1))
  I2 <- diag(ncol(K2))
	
  Kinv <- list()
  multmat <- list()
  tr <- list()

  for(r in 1:length(omegaseq)) {
    cat("r (prep.): ", r, "\n")
    K <- (omegaseq[r]*(K1 %x% I2) + (1-omegaseq[r])*(I1 %x% K2))
    Kinv[[r]] <- ginv((t(K)+K)/2)
    e <- eigen(Kinv[[r]], symmetric=TRUE)
    ecum <- cumsum(e$values)
    tr[[r]] <- length(which(ecum/sum(e$values)<truncate))
    multmat[[r]] <- e$vectors[,1:tr[[r]]]%*%diag(sqrt(e$values[1:tr[[r]]]))
  }

  simsup <- function(theta) {
    set.seed(myseed)

    tau2sample <- rep(0,R)
    if(type == "IG")
      tau2sample <- 1/rgamma(R, shape=1, scale=1/theta)
    else if(type == "SD")
      tau2sample <- rweibull(R, shape = 0.5, scale = theta)
    else if(type == "HN")
      tau2sample <- (rnorm(R, mean=0, sd=sqrt(theta)))^2

    omegasample <- sample(1:length(omegaseq), size=R, replace=TRUE, prob=omegaprob)
    supsample <- rep(0,R)

    for(r in 1:length(omegaseq)) {
      ind <- which(omegasample==r)

      betasample <- multmat[[r]]%*%matrix(rnorm(length(ind)*tr[[r]]), nrow=tr[[r]], ncol=length(ind))

      # correct for constraint
      V <- Kinv[[r]]%*%t(A)
      W <- A%*%V
      U <- ginv(W)%*%t(V)
      betasample <- betasample - t(U)%*%A%*%betasample

      fsample <- Z%*%betasample
      supsample[ind] <- apply(abs(fsample), 2, max)
    }
    supsample <- sqrt(tau2sample)*supsample
    return(supsample)
  }

  if(is.null(thetaseq))
    {
    optfn1 <- function(theta)
      {
      quantile(simsup(theta=theta), probs=alpha) - c
      }
    bopt <- uniroot(f=optfn1, interval=c(1e-8,upper=1e8), trace=1)
    bopt <- bopt$root
    res <- bopt
    }
  else
    {
    res <- rep(0, length(thetaseq))
    for(theta in 1:length(thetaseq))
      {
      cat("theta: ", thetaseq[theta], "\n")
      res[theta] <- quantile(simsup(theta=thetaseq[theta]), probs=alpha) - c
      }
    }
  return(res)
}

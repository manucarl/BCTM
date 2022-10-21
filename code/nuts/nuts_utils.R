
test_nuts <- function(beta_plus, beta_minus, r_plus, r_minus){
  beta_temp <- beta_plus - beta_minus
  temp <- (crossprod(beta_temp,r_minus) >= 0) *
    (crossprod(beta_temp, r_plus) >= 0)
  return(temp)
}

# implements Euclidian metric (qncour (2017)) by rotating and rescaling the target space, i.e. 
rotate_rescale <- function(f, gr, xx, M,  beta_y){
  
  if(!is.vector(M)) stop("only diagonal mass matrices (with positive diagonal entries) allowed")
  M_sq <- sqrt(M)
  
  # apply opposite transformation to the parameters
  f2 <- function(beta, xx) f(M_sq * beta, xx)
  gr2 <- function(beta, xx) as.vector(gr(M_sq * beta, xx) ) * M_sq
  
  # apply transformation (rotation and scaling) to simplify kinetic energy
  beta_x <- (1/M_sq) * beta_y
  
  list(f2 = f2, gr2 = gr2, beta_x = beta_x, M_sq = M_sq)
}


# this ist the FindReasonableEpsilon function in Hoffman et al(2014)
find_reasonable_epsilon <- function(beta, f, gr, xx){
  
  step_size <- 1
  r <- rnorm(n = length(beta))
  ## Do one leapfrog step
  r_star <- r + (step_size/2)*gr(beta, xx)
  beta_star <- beta + step_size * r_star
  r_star <- r_star + (step_size/2)*gr(beta_star,  xx)
  
  H1 <- calculate_hamiltonian(beta = beta,  xx=xx, r=r, f=f)
  H2 <- calculate_hamiltonian(beta = beta_star,  xx=xx, r=r_star, f=f)
  a <- 2*(exp(H2)/exp(H1)>.5)-1
  ## If jumped into bad region, a can be NaN so setup algorithm to keep
  ## halving step_size instead of throwing error
  if(!is.finite(a)) a <- -1
  k <- 1
  ## Similarly, keep going if there are infinite values
  while (!is.finite(H1) | !is.finite(H2) | a*H2-a*H1 > -a*log(2)) {
    step_size <- (2^a)*step_size
    ## Do one leapfrog step
    r_star <- r+(step_size/2)*gr(beta, xx)
    beta_star <- beta+step_size*r_star
    r_star <- r_star+(step_size/2)*gr(beta_star, xx)
    H2 <- .calculate.H(beta=beta_star,  xx=xx, r=r_star, f=f)
    k <- k+1
    if(k>500) {
      stop("Problem.")
    }
  }
  # if(verbose) message(paste("Reasonable epsilon=", step_size, "found after", k, "steps"))
  return(invisible(step_size))
}

calculate_hamiltonian <- function(beta, xx, r, f) f(beta, xx) - (1/2)*sum(r^2)


# .update_control <- function(control){
#   default <- list(adapt_delta=0.8, metric=NULL, step_size=NULL,
#                   adapt_mass=TRUE, max_treedepth=12, w1=75, w2=50, w3=25)
# 
#   new <- default
#   if(!is.null(control))
#     for(i in names(control))  new[[i]] <- control[[i]]
#   if(is.matrix(new$metric)) new$adapt_mass <- FALSE
#   return(new)
# }
# All code in this script is obtained from:
# Hothorn, T. (2017b).MLT: Most likely transformations: The mlt package, R package vignette version 0.1-5.Available at: https://CRAN.R-project.org/package=mlt.docreg.
# library("ctm")
# library("np")
# library("gamboostLSS")
library("lattice")
# library("multicore")

Ftruth <- function() {
  
  a <- list(a1 = function(x) rep(0, nrow(x)),
            a2 = function(x) x[, "x2"])
  b <- list(b1 = function(x) x[, "x1"],
            b2 = function(x) rep(.5, nrow(x)))
  
  standd <- function(x) 
    1 / rowSums(sapply(b, function(f) f(x)))
  expect <- function(x) 
    rowSums(sapply(a, function(f) f(x))) * standd(x)
  
  ret <- list(mean = expect, sd = standd)
  class(ret) <- "truth"
  ret
}



truth <- Ftruth()

dgp <- function(n = NULL, pnon = 0) {
  if (is.null(n)) {
    gr1 <- seq(from = .1, to = .9, length = 20)	
    gr2 <- seq(from = 0, to = 1.9, length = 20)	
    xdf <- expand.grid(x1 = gr1, x2 = gr2)
    n <- nrow(xdf)
    ng <- factor(rep("test", n))
  } else {
    stopifnot(length(n) == 2)
    ng <- factor(rep(c("learn", "eval"), n))
    n <- sum(n)
    xdf <- data.frame(x1 = runif(n, min = 0, max = 1),
                      x2 = runif(n, min = 0, max = 2))
  }
  if (pnon > 0) {
    xnon <- matrix(runif(n * pnon), nrow = n)
    colnames(xnon) <- paste("x", (1:pnon) + 2, sep = "")
    xdf <- cbind(xdf, xnon)
  }
  xdf$y <- rnorm(n, mean = truth$mean(xdf),
                 sd = truth$sd(xdf))
  xdf$ng <- ng
  xdf
}




xdf <- dgp()
ys <- seq(from = min(xdf$y), to = max(xdf$y), length = 100)

pdf <- function(object, newdata, q)
  UseMethod("pdf")   

pdf.truth <- function(object, newdata, q = ys) {
  m <- truth$mean(newdata)
  s <- truth$sd(newdata)
  sapply(1:length(m), function(i) 
    pnorm(q, mean = m[i], sd = s[i])
  )
}	

pdens.truth <- function(object, newdata, q = ys) {
  m <- truth$mean(newdata)
  s <- truth$sd(newdata)
  sapply(1:length(m), function(i) 
    dnorm(q, mean = m[i], sd = s[i])
  )
}

Fboost <- function(xydf, monotone = FALSE, ylin = FALSE, ngrid = 50) {
  ldf <- subset(xydf, ng == "learn")
  edf <- subset(xydf, ng == "eval")
  varn <- colnames(xydf)
  varn <- varn[grep("^x", varn)]
  for (xv in c(varn, "y")) {
    rx <- range(ldf[[xv]])
    edf <- subset(edf, edf[[xv]] > rx[1] & edf[[xv]] < rx[2])
  }
  
  if (ylin) {
    yfm <- "bols(y, lambda = 0)"
    ### df  = 8
    xfm <- paste("bbs(", varn, ", df = 8)", collapse = "+")
  } else {
    yfm <- "bbs(y, df = 2.75)"
    xfm <- paste("bbs(", varn, ", df = 2.75)", collapse = "+")
  }
  fm <- paste(yfm, xfm, sep = "~")
  fm <- as.formula(fm)
  
  ctrl <- boost_control(nu = 0.2, mstop = 1000, trace = TRUE)
  mod <- ctm(fm, data = ldf, family = Binomial(link = "probit"),
             control = ctrl, monotone = monotone, ngrid = ngrid)
  r <- outrisk(mod, newdat = edf)
  while(which.min(r) == mstop(mod) && mstop(mod) < 1200) {
    print(mstop(mod))
    mod[mstop(mod) + 100]	
    r <- outrisk(mod, newdat = edf)
  }
  ms <- which.min(r)
  ctrl$mstop <- ms
  mod <- ctm(fm, data = xydf, 
             family = Binomial(link = "probit"),
             control = ctrl, monotone = monotone, ngrid = ngrid)
  mod
}

pdf.mboost <- function(object, newdata, q = ys) {
  sapply(1:nrow(newdata), function(i)
    predict(object, newdata = newdata[i,, drop = FALSE], y = q,
            type = "response")[,1]
  )
}

Fnp <- function(xydf) {
  varn <- colnames(xydf)
  varn <- varn[grep("^x", varn)]    
  fm <- as.formula(paste("y ~ ", paste(varn, collapse = "+")))
  args <- list(bw = fm, data = xydf) 
  do.call("npcdist", args)
}

pdf.condistribution <- function(object, newdata, q = ys) {
  sapply(1:nrow(newdata), function(i) {
    tmp <- newdata[i,,drop = FALSE]
    tmp$y <- NULL
    predict(object, newdata = data.frame(y = q, tmp))
  })
}

Fgamlss <- function(xydf) {
  
  varn <- colnames(xydf)
  varn <- varn[grep("^x", varn)]    
  varn <- varn[!(varn %in% c("x1", "x2"))]
  mu.fm <- "y ~ bbs(x = x1, by = x2)"
  if (length(varn) > 0) {
    mu.fm <- paste(mu.fm, "+", paste(paste("bbs(", varn, ")"), 
                                     collapse = "+"))
  }
  mu.fm <- as.formula(mu.fm)
  si.fm <- y ~ bbs(x1)
  w <- as.integer(xydf$ng == "learn")
  
  mod <- gamboostLSS(formula = list(mu = mu.fm,
                                    sigma = si.fm),
                     families = gamboostLSS:::GaussianLSS(), 
                     data = xydf, weights = w, control = 
                       boost_control(risk = "oobag"))
  
  mr <- risk(mod, parameter = "mu")$mu    
  sr <- risk(mod, parameter = "sigma")$sigma    
  
  while((which.min(mr) == mstop(mod, parameter = "mu")) ||
        (which.min(sr) == mstop(mod, parameter = "sigma")) && max(mstop(mod)) < 300) {
    mod[mstop(mod) + 100]   
    mr <- risk(mod, parameter = "mu")$mu    
    sr <- risk(mod, parameter = "sigma")$sigma    
  }
  mr <- risk(mod, parameter = "mu")$mu    
  sr <- risk(mod, parameter = "sigma")$sigma    
  
  mod <- gamboostLSS(formula = list(mu = mu.fm,
                                    sigma = si.fm),
                     families = gamboostLSS:::GaussianLSS(), 
                     data = xydf, control = 
                       boost_control(risk = "oobag",
                                     mstop = list(mu = which.min(mr), 
                                                  sigma = which.min(sr))))
  
  
  mod
}

pdf.gamboostLSS <- function(object, newdata, q = ys) {
  
  mu <- predict(object, newdata = newdata, parameter = "mu")
  sd <- predict(object, newdata = newdata, parameter = "sigma", 
                type = "response")
  
  sapply(1:length(mu), function(i) {
    pnorm(ys, mean = mu[i], sd = sd[i])
  })
}

if (FALSE) {
  
  xdf <- dgp()
  xdf$m <- truth$m(xdf)
  xdf$sd <- truth$sd(xdf)
  wireframe(m ~ x1 + x2, data = xdf)
  wireframe(sd ~ x1 + x2, data = xdf)
  
  xdf <- dgp(c(100, 100), pnon = 0)
  boostl <- Fboost(xdf, ylin = TRUE)
  boost <- Fboost(xdf)
  kern <- Fnp(xdf)
  glss <- Fgamlss(xdf)
  x <- dgp(pnon = 0)
  pt <- pdf(truth, x)
  pbl <- pdf(boostl, x)
  pb <- pdf(boost, x)
  pk <- pdf(kern, x)
  pg <- pdf(glss, x)
  
  summary(colMeans(abs(pt - pb)))
  summary(colMeans(abs(pt - pbl)))
  summary(colMeans(abs(pt - pk)))
  summary(colMeans(abs(pt - pg)))
  
  par(ask = TRUE)
  for (i in 1:ncol(pt)) {
    print(truth$sd(x)[i])
    plot(pt[,i])
    lines(pb[,i], col = "red", lwd = 2)
    lines(pbl[,i], col = "blue", lwd = 2)
    lines(pk[,i], col = "green", lwd = 2)
    lines(pg[,i], col = "black", lwd = 2)
  }
  
}


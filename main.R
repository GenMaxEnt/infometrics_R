##### JAYNES #####
jaynes.obj <- function(lambda, x, y) {
  p_i <- exp(-lambda*x)
  Phi <- sum(p_i)
  ll  <- lambda*y + log(Phi)
  return(ll)
}

jaynes.grad <- function(lambda, x, y) {
  p_i <- exp(-lambda*x)
  Phi <- sum(p_i)
  p_i <- p_i / Phi
  foc <- y - sum(x*p_i)
  return(foc)
}

jaynes.est <- function(x, y) {
  
  lambda0 <- 0
  temp <- optim(lambda0, jaynes.obj, jaynes.grad, x = x , y = y, method = "BFGS")
  lambda <- temp$par
  p_i <- exp(-lambda*x)
  Phi <- sum(p_i)
  p_i <- p_i / Phi
  
  out <- list(
    par = p_i,
    value = temp$value,
    convergence = temp$convergence
  )
  return(out)

}

##### INVERSE #####
inverse.obj <- function(lambda, X, y, p0) {
  
  p_j  <- p0 * exp(lambda %*% X)
  Phi  <- sum(p_j)
  
  ll <- sum(lambda*y) - log(Phi)
  return(ll)
  
}

inverse.grad <- function(lambda, X, y, p0) {
  
  p_j  <- p0 * exp(lambda %*% X)
  Phi  <- sum(p_j)
  p_j  <- p_j/Phi
  
  grad <- y - X %*% t(p_j)
  return(grad)
  
}

inverse.test <- function(formula, data, p0) {
  
  mf <- model.frame(formula, data)
  mt <- attr(mf, "terms")
  
  y <- model.response(mf)
  X <- model.matrix(mt, mf)
  
  n <- length(y)
  p <- ncol(X)
  
  if (missing(p0)) {
    p0 <- rep(1/p, p)
  }
  
  lambda0 <- rep(0, n)
  temp <- optim(lambda0, inverse.obj, inverse.grad, 
                y = y, X = X, p0 = p0, 
                method = "BFGS", control = list(fnscale = -1))
  
  lambda <- temp$par
  
  p_j   <- p0 * exp(lambda %*% X)
  Phi   <- sum(p_j)
  bhat  <- p_j/Phi
  
  yhat  <- y - X %*% t(bhat)
  ehat  <- y - yhat
  df.residual  <- n - p
  
  sigma2 <- sum(ehat^2) / df.residual
  vcov <- sigma2 * chol2inv(t(X) %*% X)

  out <- list(
    call = match.call(),
    terms = mt,
    coefficients = bhat,
    residuals = ehat,
    fitted.values = yhat,
    rank = p,
    df.residual = df.residual,
    sigma = sqrt(sigma2),
    vcov = vcov,
    model = mf,
    xlevels = .getXlevels(mt, mf),
    contrasts = attr(X, "contrasts")
  )
  
  class(out) <- "jaynes"
  out
}

summary.jaynes <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  tval <- object$coefficients / se
  pval <- 2 * pt(abs(tval), df = object$df.residual, lower.tail = FALSE)
  
  coefmat <- cbind(
    Estimate = object$coefficients,
    Std.Error = se,
    t.value = tval,
    Pr = pval
  )
  
  out <- list(
    call = object$call,
    coefficients = coefmat,
    sigma = object$sigma,
    df = c(object$rank, object$df.residual),
    r.squared = 1 - sum(object$residuals^2) /
      sum((object$fitted.values - mean(object$fitted.values))^2)
  )
  
  class(out) <- "summary.jaynes"
  out
}

print.summary.jaynes <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  
  cat("\nCoefficients:\n")
  printCoefmat(
    x$coefficients,
    digits = max(3L, getOption("digits") - 2L),
    signif.stars = TRUE
  )
  
  cat(
    "\nResidual standard error:",
    formatC(x$sigma, digits = 4),
    "on", x$df[2], "degrees of freedom\n"
  )
  
  cat(
    "Multiple R-squared:",
    formatC(x$r.squared, digits = 4), "\n"
  )
  
  invisible(x)
}

inverse.pure.obj <- function(lambda, X, y, p0) {
  
  p_j  <- p0 * exp(lambda %*% X)
  Phi  <- sum(p_j)

  ll <- sum(lambda*y) - log(Phi)
  return(ll)
  
}

inverse.pure.grad <- function(lambda, X, y, p0) {
  
  p_j  <- p0 * exp(lambda %*% X)
  Phi  <- sum(p_j)
  p_j  <- p_j/Phi

  grad <- y - X %*% t(p_j)
  return(grad)
  
}

inverse.pure <- function(X, y, p0) {
  
  tt <- dim(X)
  n <- tt[1]
  J <- tt[2]
  
  if (missing(p0)) {
    p0 <- rep(1/J, J)
  }
  
  lambda0 <- rep(0, n)
  temp <- optim(lambda0, inverse.pure.obj, inverse.pure.grad, 
                y = y, X = X, p0 = p0, 
                method = "BFGS", control = list(fnscale = -1))
  
  lambda <- temp$par
  
  p_j  <- p0 * exp(lambda %*% X)
  Phi  <- sum(p_j)
  p_j  <- p_j/Phi

  out <- list(
    lambda = lambda,
    p = p_j,
    value = temp$value,
    convergence = temp$convergence
  )
  
  return(out)
  
}

inverse.noise.obj <- function(lambda, X, y, p0, w0, V, nu) {
  
  p_j  <- p0 * exp(lambda %*% X / (1 - nu))
  Phi  <- sum(p_j)
  w_im <- w0 * exp(lambda*V/nu)
  Psi  <- apply(w_im, 1, sum)
  
  ll <- sum(lambda*y) - (1-nu)*log(Phi) - nu*sum(log(Psi))
  return(ll)
  
}

inverse.noise.grad <- function(lambda, X, y, p0, w0, V, nu) {
  
  p_j  <- p0 * exp(lambda %*% X / (1 - nu))
  Phi  <- sum(p_j)
  p_j  <- p_j/Phi
  w_im <- w0 * exp(lambda*V/nu)
  Psi  <- apply(w_im, 1, sum)
  w_im <- w_im/Psi
  e_i  <- apply(w_im * V, 1, sum)
  
  grad <- y - X %*% t(p_j) - e_i
  return(grad)

}

inverse.noise <- function(X, y, p0, w0, V, nu) {
  
  tt <- dim(X)
  n <- tt[1]
  J <- tt[2]
  
  if (missing(p0)) {
    p0 <- rep(1/J, J)
  }
  
  if (missing(w0)) {
    w0 <- matrix(1/3, n, 3)
  }
  
  if (missing(V)) {
    V <- matrix(0, n, 3)
    for (i in 1:n) {
      stdev <- sd(X[i,])
      V[i,1] <- -3*stdev
      V[i,3] <-  3*stdev
    }
  }
  
  if (missing(nu)) {
    nu <- 0.5
  }
  
  lambda0 <- rep(0, n)
  temp <- optim(lambda0, inverse.noise.obj, inverse.noise.grad, 
                y = y, X = X, p0 = p0, w0 = w0, V = V, nu = nu, 
                method = "BFGS", control = list(fnscale = -1))
  
  lambda <- temp$par
    
  p_j  <- p0 * exp(lambda %*% X / (1 - nu))
  Phi  <- sum(p_j)
  p_j  <- p_j/Phi
  w_im <- w0 * exp(lambda*V/nu)
  Psi  <- apply(w_im, 1, sum)
  w_im <- w_im/Psi
  e_i  <- apply(w_im * V, 1, sum)
  
  out <- list(
    lambda = lambda,
    p = p_j,
    value = temp$value,
    convergence = temp$convergence
  )
  
  return(out)

}

##### REG #####
gce.reg.obj <- function(lambda, X, y, p0, Z, w0, V, nu) {
  
  p_km <- p0 * exp(Z * as.vector(lambda %*% X) / (1-nu))
  Phi  <- apply(p_km, 1, sum)
  w_tm <- w0 * exp(as.matrix(lambda) %*% t(as.matrix(V)) / nu)
  Psi  <- apply(w_tm, 1, sum)

  ll <- sum(lambda*y) - (1-nu)*sum(log(Phi)) - nu*sum(log(Psi))
  return(ll)
  
}

gce.reg.grad <- function(lambda, X, y, p0, Z, w0, V, nu) {
  
  p_km <- p0 * exp(Z * as.vector(lambda %*% X) / (1-nu))
  Phi  <- apply(p_km, 1, sum)
  p_km <- p_km / as.vector(Phi)
  beta <- apply(Z * p_km, 1, sum)
  w_tm <- w0 * exp(as.matrix(lambda) %*% t(as.matrix(V)) / nu)
  Psi  <- apply(w_tm, 1, sum)
  w_tm <- w_tm / as.vector(Psi)
  e_t  <- w_tm %*% as.matrix(V)
  
  grad <- y - X %*% as.matrix(beta) - e_t
  return(grad)
  
}

gce.reg <- function(X, y, p0, Z, w0, V, nu) {
  
  # note X, p0, Z w0 are matrices, y and V are vectors
  h <- dim(X)
  n <- h[1]
  K <- h[2]
  M <- dim(Z)[2]
  
  if (missing(p0)) {
    p0 <- rep(1/M, M)
    p0 <- matrix(p0, nrow = K, ncol = M, byrow = TRUE)
  }
  
  if (missing(w0)) {
    J <- 3
    w0 <- rep(1/J, J)
    w0 <- matrix(w0, nrow = n, ncol = J, byrow = TRUE)
  }
  
  if (missing(V)) {
    sd_y <- sd(y)
    V <- c(-3*sd_y, 0, 3*sd_y)
  }

  if (missing(nu)) {
    nu <- 0.5
  }
  
  lambda0 <- rep(0, n)
  temp <- optim(lambda0, gce.reg.obj, gce.reg.grad, 
                X=X, y=y, p0=p0, Z=Z, w0=w0, V=V, nu=nu,
                method = "BFGS", control = list(fnscale = -1))
  
  lambda <- temp$par
  p_km <- p0 * exp(Z * as.vector(lambda %*% X) / (1-nu))
  Phi  <- apply(p_km, 1, sum)
  p_km <- p_km / as.vector(Phi)
  beta <- apply(Z * p_km, 1, sum) 
  w_tm <- w0 * exp(as.matrix(lambda) %*% t(as.matrix(V)) / nu)
  Psi  <- apply(w_tm, 1, sum)
  w_tm <- w_tm / as.vector(Psi)
  e_t  <- w_tm %*% as.matrix(V)
  
  return(beta)
  
}

##### TEST JAYNES #####
x <- 1:6
y <- 2

tt <- jaynes.est(x, y)
plot(tt$par, type = "h")

##### TEST INVERSE #####
rm(list=ls())
n <- 100
J <- 3
set.seed(256)
X <- matrix(runif(n*J),n,J)
#X <- matrix(1:12,n,J)
p0 <- rep(1/J, J)
w0 <- matrix(1/3, n, 3)
V <- matrix(0, n, 3)
for (i in 1:n) {
  stdev <- sd(X[i,])
  V[i,1] <- -3*stdev
  V[i,3] <-  3*stdev
}
nu <- 0.5

lambda <- rep(0,n)
lambda <- 1:n
p <- c(.8, .05, .15)
y <- X %*% p + rnorm(n, 0, .05)

mydata <- data.frame(y=y, x1=X[,1], x2=X[,2], x3=X[,3])

s3.test <- inverse.test(y~x1+x2+x3-1, mydata)
summary(s3.test)

zz <- inverse.pure(X, y)
zz$p

tt <- inverse.noise(X, y)
tt$p

ttt <- inverse.noise(X, y, p0, w0, V, .7)
ttt$p

##### TEST REGRESSION #####
rm(list=ls())
n <- 50
K <- 3
M <- 3
J <- 3
set.seed(256)
X <- matrix(runif(n*(K-1)), n, (K-1))
X <- cbind(1, X)
p0 <- rep(1/M, M)
p0 <- matrix(p0, nrow = K, ncol = M, byrow = TRUE)
w0 <- rep(1/J, J)
w0 <- matrix(w0, nrow = n, ncol = J, byrow = TRUE)
beta <- matrix(c(2, 1, 3), ncol = 1)
Z <- matrix(c(-5,0,5), nrow = K, ncol = M, byrow = TRUE)
y <- X %*% beta + rnorm(n, 0, 0.05)
y <- as.vector(y)
sd_y <- sd(y)
V <- c(-3*sd_y, 0, 3*sd_y)
nu <- 0.5

gce.reg(X,y,,Z)

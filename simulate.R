simulate <- function(beta0, sig20, lambda0, n, tau, Gr) {
  
  nbeta <- length(beta0)
  
  b <- sqrt(sig20) * rnorm(n, 0, 1)
  NN <- rep(2, n)
  AA <- runif(n, 0, 1)
  NN[AA < 0.2] <- 1
  NN[AA > 0.9] <- 3
  
  W <- list()
  X1 <- list()
  X2 <- list()
  ID <- numeric()
  L <- numeric()
  U <- numeric()
  
  for (i in 1:n) {
    ni <- NN[i]
    Z0 <- runif(1, 0, 1)
    Z1 <- runif(1, 0, 1)
    Z2 <- runif(ni, 0, 1)
    bi <- b[i]
    X1i <- cbind(Z0 < 0.5, Z2)
    X2i <- cbind(Z1 > 0.5, Z2)
    
    Xbeta1 <- X1i %*% beta0 + bi
    Xbeta2 <- X2i %*% beta0 + bi
    
    Wi <- rep(8, ni)
    
    Ti <- -log(runif(ni, 0, 1))
    G_out <- Gtransform(Ti, Gr)  # Assuming this function is available in your R environment
    temp1 <- (exp(Xbeta1) * log(1 + lambda0 * Wi))
    ind <- G_out$G < temp1
    
    Ti <- (ind * G_out$G / exp(Xbeta1)) + 
      (!ind * ((G_out$G - temp1) / exp(Xbeta2) + log(1 + lambda0 * Wi)))
    Ti <- (exp(Ti) - 1) / lambda0
    
    TLi <- runif(ni, 0, 1) * 3
    TUi <- TLi + 0.1 + runif(ni, 0, 1)
    
    Li <- TLi
    Ui <- TUi
    
    Li[Ti < TLi] <- 0
    Ui[Ti < TLi] <- TLi[Ti < TLi]
    Li[Ti > TUi] <- TUi[Ti > TUi]
    Ui[Ti > TUi] <- 99999
    
    ID <- c(ID, rep(i, ni))
    L <- c(L, Li)
    U <- c(U, Ui)
    W <- c(W, Wi)
    X1 <- rbind(X1, X1i)
    X2 <- rbind(X2, X2i)
  }
  
  up <- max(U[U < tau])
  while (up > max(L)) {
    U[U == up] <- 99999
    up <- max(U[U < tau])
  }
  
  time <- unique(c(L[L > 0], U[U <= tau]))
  m <- length(time)
  
  indW <- outer(W, time, ">")
  Xtime <- vector("list", length = nbeta)
  
  for (s in 1:nbeta) {
    Xtime[[s]] <- X1[, s, drop = FALSE] * indW + X2[, s, drop = FALSE] * (!indW)
  }
  
  return(list(ID = ID, L = L, U = U, time = time, Xtime = Xtime))
}

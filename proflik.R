proflik <- function(par, templambda, time, Xtime, TindL, TindU, TindR, 
                    U, tau, Gr, ys, mus, nbeta, bt, wt, n, indID) {
  
  N <- length(ys)
  NN <- length(U)
  m <- length(time)
  
  oldbeta <- par[1:nbeta]
  oldsig2 <- par[nbeta + 1]
  
  bb <- sqrt(2 * oldsig2) * ys
  Ebb <- exp(bb)
  
  Xbeta <- matrix(0, nrow = NN, ncol = m)
  for (p in 1:nbeta) {
    Xbeta <- Xbeta + Xtime[[p]] %*% oldbeta[p]
  }
  
  EXbeta <- exp(Xbeta)
  
  oldlambda <- templambda
  
  maxiter <- 5000
  error <- 1
  iter <- 0
  
  while (iter < maxiter & error > 0.005) {
    
    lambdaEXbeta <- oldlambda * EXbeta
    
    S1 <- (t(t(TindL * lambdaEXbeta) %*% rep(1, m)) %*% Ebb)
    S2 <- (t(t(TindU * lambdaEXbeta) %*% rep(1, m)) %*% Ebb)
    S2[U > tau, ] <- 99999 + max(max(S1))
    
    G1 <- Gtransform(S1, Gr)$G
    G2 <- Gtransform(S2, Gr)$G
    
    Em <- exp(indID * (-G1 + log(1 - exp(-G2 + G1))))
    
    EpEm <- t(indID) %*% Em
    D <- EpEm / (exp(-G1) - exp(-G2))
    
    EpEm <- EpEm * (mus / (EpEm %*% mus))
    
    v <- rep(0, NN)
    w <- matrix(0, nrow = NN, ncol = m)
    
    ind <- matrix(U <= tau, ncol = N, nrow = length(U), byrow = TRUE)
    
    temp0 <- lambdaEXbeta
    DS <- S2 - S1
    DS[DS <= 0] <- 99999
    
    if (Gr > 0) {
      r <- Gr
      A <- 1/r + S1
      C <- 1/r + S2
      f <- sum(wt) * ((A^(-1/r) - ind * C^(-1/r)) * D %*% mus)
      v <- (sum(wt * bt) * ((A^(-1/r - 1) - ind * C^(-1/r - 1)) * D) %*% (Ebb * mus)) / f
      
      temp1 <- matrix(0, nrow = NN, ncol = N)
      temp2 <- matrix(0, nrow = NN, ncol = N)
      for (q in 1:length(bt)) {
        temp1 <- temp1 + bt[q] * wt[q] * (1 / (1 - exp(-DS * bt[q] / A)))
        temp2 <- temp2 + bt[q] * wt[q] * (1 / (1 - exp(-DS * bt[q] / C)))
      }
      
      temp3 <- ((A^(-1/r - 1) * temp1 - ind * C^(-1/r - 1) * temp2) * D) %*% (Ebb * mus)
      temp <- temp0 * (temp3 / f)
    } else {
      v <- EpEm %*% Ebb
      temp <- temp0 * ((1 / (1 - exp(-DS)) * EpEm) %*% Ebb)
    }
    
    w[TindL == 0 & TindU == 1,] <- temp[TindL == 0 & TindU == 1,]
    temp4 <- temp0 * v
    w[TindU == 0,] <- temp4[TindU == 0,]
    
    tempind <- (U > tau) & (TindL == 0)
    w[tempind,] <- temp4[tempind,]
    
    newlambda <- (colSums(w * TindR) / (v %*% (EXbeta * TindR)))
    error <- sum(abs(newlambda - oldlambda))
    
    iter <- iter + 1
    
    oldlambda <- newlambda
  }
  
  f <- loglik(oldbeta, oldsig2, newlambda, time, Xtime, TindL, TindU, TindR, U, tau, Gr, ys, mus, indID)
  
  return(f)
}

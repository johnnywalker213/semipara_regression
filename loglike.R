loglik <- function(oldbeta, oldsig2, oldlambda, time, Xtime, TindL, 
                   TindU, TindR, U, tau, Gr, ys, mus, indID) {
  
  bb <- sqrt(2 * oldsig2) * ys
  Ebb <- exp(bb)
  
  m <- length(time)
  nbeta <- length(oldbeta)
  NN <- length(U)
  
  Xbeta <- matrix(0, nrow = NN, ncol = m)
  
  for (p in 1:nbeta) {
    Xbeta <- Xbeta + Xtime[[p]] %*% oldbeta[p]
  }
  
  EXbeta <- exp(Xbeta)
  lambdaEXbeta <- oldlambda * EXbeta
  
  S1 <- t(t(TindL * lambdaEXbeta) %*% rep(1, m)) %*% Ebb
  S2 <- t(t(TindU * lambdaEXbeta) %*% rep(1, m)) %*% Ebb
  
  S2[U > tau, ] <- 99999 + max(max(S1))
  
  Gtransform <- function(x, r) {
    if (r == 0) {
      G <- x
      dG <- rep(1, length(x))
      d2G <- rep(0, length(x))
      IG <- x
    } else {
      G <- log(1 + r * x) / r
      dG <- 1 / (1 + r * x)
      d2G <- -r / ((1 + r * x) ^ 2)
      IG <- (exp(r * x) - 1) / r
    }
    return(list(G = G, dG = dG, d2G = d2G, IG = IG))
  }
  
  G1 <- Gtransform(S1, Gr)$G
  G2 <- Gtransform(S2, Gr)$G
  
  Em <- exp(indID * (-G1 + log(1 - exp(-G2 + G1))))
  
  fvec <- (log((Em * mus) / sqrt(3.14159265358)))
  
  return(fvec)
}

# Example usage:
# Ensure you pass appropriate arguments when calling the function
# result <- loglik(oldbeta, oldsig2, oldlambda, time, Xtime, TindL, TindU, TindR, U, tau, Gr, ys, mus, indID)

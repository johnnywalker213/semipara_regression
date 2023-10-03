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

# Example usage:
# x <- 1:5
# r <- 2
# result <- Gtransform(x, r)
# G <- result$G
# dG <- result$dG
# d2G <- result$d2G
# IG <- result$IG

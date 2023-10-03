GammaInt <- function(aa, bb, r, bt, wt) {
  
  # f = 1./(1-exp(-aa)).*exp(-bb);
  f <- exp(-bb) / (1 - exp(-aa))
  
  if (r > 0) {
    # initialize f with a zero vector of the same length as aa
    f <- 0 * aa
    
    # for k=1:1:length(bt)
    for (k in 1:length(bt)) {
      # x=bt(k)./(bb+1/r);
      x <- bt[k] / (bb + 1/r)
      
      # f=f+wt(k).*x./(1-exp(-x.*aa)).*power(bb+1/r, -1/r)/gamma(1/r)/power(r,1/r);
      f <- f + wt[k] * x / (1 - exp(-x * aa)) * 
        ((bb + 1/r)^(-1/r)) / gamma(1/r) / (r^(1/r))
    }
  }
  
  return(f)
}

#' @export
Lamb_gmonl <- function(times, subj, X, y, d, tau, kn, degree, lambda,omega,range,mv){
  dim = length(subj)
  W = Weight_Ni(y, subj)$W
  X = matrix(X, nrow = dim)
  px = ncol(X)
  n = length(unique(subj))
  lambda = unique(lambda)
  nlam = length(lambda)
  dim = length(y)
  yhatsic = matrix(NA, dim, nlam)
  ressic = matrix(NA, dim, nlam)
  SIC = NULL
  plam = NULL
  for (i in 1:nlam) {
    qvcsic = qrvcp_gmonl(times, subj, y, X, tau, kn, degree, 
                         lambda = lambda[i], d,omega,range,mv)
    hat_btsic = qvcsic$hat_bt
    yhatsic_k = matrix(NA, dim, px)
    for (k in 1:px) {
      yhatsic_k[, k] = hat_btsic[seq((k - 1) * dim + 1, 
                                     k * dim)] * X[, k]
    }
    yhatsic[, i] = rowSums(yhatsic_k)
    ressic[, i] = y - yhatsic[, i]
    plam[i] = length(ressic[which(abs(ressic[, i]) < (10^(-2)))])
    SIC[i] = log(sum(W * ressic[, i] * (tau - 1 * (ressic[, 
                                                          i] < 0)))/n) + log(dim) * plam[i]/(dim * 2)
  }
  lambdasic = lambda[which(SIC == min(SIC))]
  Lout = list(lambdasic = lambdasic)
  return(Lout)
}

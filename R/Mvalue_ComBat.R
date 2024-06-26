# Adjust for batch effects in DNA methylation data by converting beta-values to M-values followed by ComBat
Mvalue_ComBat <- function(bv, batch, group = NULL, 
                          covar_mod = NULL, full_mod = TRUE, 
                          mean.only = FALSE, pseudo_beta = 1e-4,
                          ref.batch = NULL) {
  ## convert extreme 0 or 1 values to pseudo-beta
  if (pseudo_beta <= 0 | pseudo_beta >= 0.5) {
    stop("Invalid pseudo beta-values.")
  }
  bv[bv == 0] <- pseudo_beta
  bv[bv == 1] <- 1 - pseudo_beta
  
  ## convert beta values to M values
  mv <- log(bv / (1 - bv))
  
  ## construct the model matrix
  mod <- covar_mod
  if (full_mod) {
    mod <- cbind(group, mod)
  }
  
  ## run ComBat
  if (mean.only) {
    adj_mv <- ComBat(mv, batch, mod = mod, mean.only = TRUE, 
                     par.prior = TRUE, ref.batch = ref.batch)
  } else {
    adj_mv <- ComBat(mv, batch, mod = mod, mean.only = FALSE, 
                     par.prior = TRUE, ref.batch = ref.batch)
  }
  
  ## convert adjusted M values back to beta values
  adj_bv <- exp(adj_mv) / (1 + exp(adj_mv))
  return(adj_bv)
}

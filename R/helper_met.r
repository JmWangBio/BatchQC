# No longer necessary when integrated into sva
####  Expand a vector into matrix (columns as the original vector)
vec2mat <- function(vec, n_times){
  return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}

# To be integrated into sva
#### Monte Carlo integration functions for beta distribution
monte_carlo_int_beta <- function(dat, mu, gamma, phi, delta, feature.subset.n) {
  weights <- pos_res <- list()
  for (i in 1:nrow(dat)) {
    m <- mu[-i, !is.na(dat[i, ])]
    p <- phi[-i, !is.na(dat[i, ])]
    x <- dat[i, !is.na(dat[i, ])]
    gamma_sub <- gamma[-i]
    delta_sub <- delta[-i]
    
    # take a subset of features to do integration - save time
    if (!is.null(feature.subset.n) & is.numeric(feature.subset.n) & length(feature.subset.n) == 1) {
      if (i == 1) {
        cat(sprintf("Using %s random features for Monte Carlo integration\n", feature.subset.n))
      }
      mcint_ind <- sample(1:(nrow(dat) - 1), feature.subset.n, replace = FALSE)
      m <- m[mcint_ind, ]
      p <- p[mcint_ind, ]
      gamma_sub <- gamma_sub[mcint_ind]
      delta_sub <- delta_sub[mcint_ind]
      G_sub <- feature.subset.n
    } else {
      if (i == 1) {
        cat(
          "Using all features for Monte Carlo integration; 
        the function runs very slow for large number of features\n"
        )
      }
      G_sub <- nrow(dat) - 1
    }
    
    LH <- sapply(1:G_sub, function(j) {
      prod(stats::dbeta(x, shape1 = m[j, ] * (p[j, ] * exp(delta_sub[j])), 
                        shape2 = (1 - m[j, ]) * (p[j, ] * exp(delta_sub[j]))))
    })
    LH[is.na(LH)] <- 0
    if (sum(LH) == 0 | is.na(sum(LH))){
      pos_res[[i]] <- c(gamma.star = as.numeric(gamma[i]),
                        delta.star = as.numeric(delta[i]))
    } else {
      pos_res[[i]] <- c(gamma.star = NA, 
                        delta.star = NA)
      pos_res[[i]]["gamma.star"] <- sum(gamma_sub * LH, na.rm = TRUE) / 
        sum(LH[!is.na(gamma_sub)])
      pos_res[[i]]["delta.star"] <- sum(delta_sub * LH, na.rm = TRUE) / 
        sum(LH[!is.na(delta_sub)])
    }
    
    weights[[i]] <- as.matrix(LH / sum(LH))
  }
  pos_res <- do.call(rbind, pos_res)
  res <- list(gamma_star = pos_res[, "gamma.star"], 
              delta_star = pos_res[, "delta.star"]) 
  return(res)
} 

# To be integrated into sva
#### Match quantiles for beta distribution
match_quantiles_beta <- function(bv_sub, old_mu, old_phi, new_mu, new_phi) {
  new_bv_sub <- matrix(NA, nrow = nrow(bv_sub), ncol = ncol(bv_sub))
  for (a in 1:nrow(bv_sub)) {
    for (b in 1:ncol(bv_sub)) {
      tmp_p <- stats::pbeta(bv_sub[a, b], shape1 = old_mu[a, b] * old_phi[a, b], 
                            shape2 = (1 - old_mu[a, b]) * old_phi[a, b])
      if (is.na(tmp_p)) {
        new_bv_sub[a, b] <- bv_sub[a, b]  
      } else {
        new_bv_sub[a, b] <- stats::qbeta(tmp_p, shape1 = new_mu[a, b] * new_phi[a, b], 
                                         shape2 = (1 - new_mu[a, b]) * new_phi[a, b])
      }
    }
  }
  return(new_bv_sub)
}

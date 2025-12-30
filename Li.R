#load packages in main first

###CRE
worst_case_randomization_test <- function(Y,Z,M, delta, assumption="Li_A1") {

  n <- length(Y)
  n_t <- sum(Z)
  n_c <- n - n_t
  
  if (n != length(Z) || n != length(M)) {
    stop("Y, Z, and M must have the same length")}
  if (!all(Z %in% c(0,1))) stop("Z must be 0 or 1")
  if (!all(M %in% c(0,1))) stop("M must be 0 or 1")
  if (n_t == 0 || n_c == 0) stop("need at least one treated and one control")
  if (!assumption %in% c("Li_A1", "Li_A2", "Li_A3", "Li_A4")) stop("Invalid assumption")
  if (length(delta) == 1) {
    delta <- rep(delta, n)
  }
  # Treated and control indicators
  treated_idx <- which(Z == 1)
  control_idx <- which(Z == 0)
  
  if (assumption == "Li_A4") {
    # Li_A4: ignore missingness entirely
    obs_idx <- which(M == 0)                # only observed outcomes
    Z_obs <- Z[obs_idx]
    Y_obs <- Y[obs_idx] - Z_obs*delta[obs_idx]    # adjust for delta
    
    # Wilcoxon test
    wilcox_result <- wilcox.test(Y_obs[Z_obs == 1], Y_obs[Z_obs == 0], 
                                  alternative = "greater", exact = FALSE, continuity = TRUE)
    
    return(list(p_value = wilcox_result$p.value))
  }
  
  # Impute worst-case Y(0)
  Y_0_tilde <- numeric(n)
  Y_0_tilde <- getY_0_tilde(Y,Z,M,delta, assumption)

  # Wilcoxon test
  wilcox_result <- wilcox.test(Y_0_tilde[Z == 1], Y_0_tilde[Z == 0], 
                               alternative = "greater", exact=FALSE, continuity=TRUE)
  
  return(list(p_value = wilcox_result$p.value))
}
getY_0_tilde <- function(Y, Z, M, delta, assumption="Li_A1") {
  n <- length(Y)
  if (length(delta) == 1) delta <- rep(delta, n)
  big <- 1e7
  
  # fill NAs temporarily
  Y_filled <- Y
  Y_filled[is.na(Y_filled)] <- 0
  
  Y_tilde <- numeric(n)
  
  if (assumption == "Li_A1") {
    Y_tilde[Z==1 & M==0] <- Y_filled[Z==1 & M==0] - delta[Z==1 & M==0]
    Y_tilde[Z==1 & M==1] <- -big
    Y_tilde[Z==0 & M==0] <- Y_filled[Z==0 & M==0]
    Y_tilde[Z==0 & M==1] <- big
  } else if (assumption == "Li_A2") {
    Y_tilde[Z==1 & M==0] <- Y_filled[Z==1 & M==0] - delta[Z==1 & M==0]
    Y_tilde[Z==1 & M==1] <- big
    Y_tilde[Z==0 & M==0] <- Y_filled[Z==0 & M==0]
    Y_tilde[Z==0 & M==1] <- big
  } else { # Li_A3
    Y_tilde[Z==1 & M==0] <- Y_filled[Z==1 & M==0] - delta[Z==1 & M==0]
    Y_tilde[Z==1 & M==1] <- -big
    Y_tilde[Z==0 & M==0] <- Y_filled[Z==0 & M==0]
    Y_tilde[Z==0 & M==1] <- -big
  }
  
  return(Y_tilde)
}

###2-STEP-procedure (doesn't work)
worst_case_randomization_test_two_step <- function(Y, Z, M, delta, beta = 0.04, assumption = "Li_A2") {
  
  n <- length(Y)
  
  if (length(delta) == 1) {
    delta <- rep(delta, n)
  }
  n_t <- sum(Z)   # treated
  n_c <- n - n_t  # controls
  n_11 <- sum(Z == 1 & M == 0)  # treated observed
  n_01 <- sum(Z == 0 & M == 0)  # control observed
  
  if (n != length(Z) || n != length(M)) stop("Y, Z, and M must have the same length")
  if (!all(Z %in% c(0,1))) stop("Z must be 0 or 1")
  if (!all(M %in% c(0,1))) stop("M must be 0 or 1")
  if (assumption != "Li_A2") stop("Two-step procedure only for monotone2 assumption")
  if (sum(M) == 0) stop("no missing values...")
  
  
  # Bound from Wang
  M_hat <- wang_upper_M(n_01, n_c, n, beta)
  m <- max(n_11 + n_01 - M_hat, 0)
  
  # Ai, Bi and Di
  treated_observed_idx <- which(Z == 1 & M == 0)
  D <- numeric(length(treated_observed_idx))
  big <- .Machine$double.xmax/10
  
  for (k in seq_along(treated_observed_idx)) {
    i <- treated_observed_idx[k]
    
    # Case A_i 
    Y_A <- numeric(n)
    for (j in 1:n) {
      if (Z[j] == 1 & M[j] == 0 & j == i) Y_A[j] <- big
      else if (Z[j] == 1 & M[j] == 0) Y_A[j] <- Y[j] - delta[j]
      else if (Z[j] == 1 & M[j] == 1) Y_A[j] <- big
      else if (Z[j] == 0 & M[j] == 0) Y_A[j] <- Y[j]
      else if (Z[j] == 0 & M[j] == 1) Y_A[j] <- big
    }
    A_i <- sum(Z * rank(Y_A, ties.method = "average"))
    
    # Case B_i
    Y_B <- numeric(n)
    for (j in 1:n) {
      if (Z[j] == 1 & M[j] == 0) Y_B[j] <- Y[j] - delta[j]
      else if (Z[j] == 1 & M[j] == 1) Y_B[j] <- big
      else if (Z[j] == 0 & M[j] == 0) Y_B[j] <- Y[j]
      else if (Z[j] == 0 & M[j] == 1) Y_B[j] <- big
    }
    B_i <- sum(Z * rank(Y_B, ties.method = "average"))
    
    D[k] <- A_i - B_i
  }
  
  # Sort treated observed by Di
  sorted_indices <- treated_observed_idx[order(D)]
  
  # Impute
  Y_tilde_two_step <- numeric(n)
  for (i in 1:n) {
    if (Z[i] == 1 & M[i] == 1) {  # treated missing
      Y_tilde_two_step[i] <- big
    } else if (Z[i] == 1 & M[i] == 0) {  # treated observed
      if (i %in% sorted_indices[1:m]) {  # unfavorable 
        Y_tilde_two_step[i] <- big
      } else {                          # favorable
        Y_tilde_two_step[i] <- Y[i] - delta[i]
      }
    } else if (Z[i] == 0 & M[i] == 0) {  # control observed
      Y_tilde_two_step[i] <- Y[i]
    } else if (Z[i] == 0 & M[i] == 1) {  # control missing
      Y_tilde_two_step[i] <- big
    }
  }
  
  # Wilcoxon test
  wilcox_result <- wilcox.test(Y_tilde_two_step[Z == 1], Y_tilde_two_step[Z == 0], alternative = "greater")
  
  
  p_value_two_step <- wilcox_result$p.value + beta
  
  return(list(
    p_value = p_value_two_step,
    t_obs = wilcox_result$statistic,
    Y_tilde = Y_tilde_two_step,
    M_hat = M_hat,
    m = m
  ))
}
wang_upper_M <- function(x, n0, N, beta = 0.04) {
  
  #   x     = observed number of "successes" (here n_01, i..e number of control missing)
  #   n0    = size of control group
  #   N     = total population size (i.e small n)
  #   beta  = significance level 

  if (x <= 0) return(Inf)
  if (x >= n0) return(n0)
  
  M_candidates <- seq(x, N)
  
  # Compute cumulative probabilities until coverage >= 1 - beta
  for (M in M_candidates) {
    p_val <- phyper(x - 1, M, N - M + 1, n0)  # Eq.(10) in Wang (2015)
    if (p_val >= 1 - beta) return(M)
  }
  
  stop("couldn't bound M with given beta, try with different beta")
}

###SRE
worst_case_stratified_randomization_test <- function(
    data, delta_blocks, assumption="Li_A1",
    method="median_aligned_ranks", L=2000) {
  
  Z             <- data$Z
  block         <- data$block
  prep          <- prepare_Z_permutations(Z, block, L)
  Z_perms       <- prep$Z_perms       
  block_indices <- prep$block_indices
  
  n             <- length(Z)
  K             <- length(block_indices)
  
  if (length(delta_blocks) == 1) {
    delta <- rep(delta_blocks, n)
  } else if (length(delta_blocks) == K) {
    delta <- numeric(n)
    for (k in seq_len(K)) delta[block_indices[[k]]] <- delta_blocks[k]
  } else stop("delta_blocks has wrong length.")
  
  if (method == "combined_ranks") {
    
    Y_0_tilde <- getY_0_tilde(data$Y_obs, Z, data$M_obs, delta, assumption)
    
    # create temp data frames so coin's formula finds variables
    df_oracle <- data.frame(Y = data$Y_true, Z = as.factor(Z), block = as.factor(block))
    oracle_test <- coin::wilcox_test(Y ~ Z | block, data = df_oracle,
                                     distribution = approximate(nresample = 10000),
                                     alternative = "less")
    p_oracle <- as.numeric(pvalue(oracle_test))
    
    df_non <- data.frame(Y = Y_0_tilde, Z = as.factor(Z), block = as.factor(block))
    van_elteren_test <- coin::wilcox_test(Y ~ Z | block, data = df_non,
                                          distribution = approximate(nresample = 10000),
                                          alternative = "less")
    p_non <- as.numeric(pvalue(van_elteren_test))
    
    return(list(
      p_oracle = p_oracle,
      p_value  = p_non
    ))
  }
  else {
    #oracle ranks
    Y_true_aligned <- data$Y_true - ave(data$Y_true, block)
    Y_ranks_oracle <- rank(Y_true_aligned, ties.method = "average")
    obs_T_oracle   <- sum(Y_ranks_oracle[Z == 1])
    
    #ranks
    Y_aligned  <- getY_aligned(data$Y_obs, data$M_obs, Z, block, delta, method = method) 
    Y_0_tilde  <- getY_0_tilde(Y_aligned, Z, data$M_obs, delta, assumption)
    Y_ranks_non <- rank(Y_0_tilde, ties.method = "average")
    obs_T_non   <- sum(Y_ranks_non[Z == 1])
    
    #permutation distribution
    T_perm_oracle <- numeric(L)
    T_perm_non    <- numeric(L)
    
    for (i in seq_len(L)) {
      Zp <- Z_perms[[i]]
      T_perm_oracle[i] <- sum(Y_ranks_oracle[Zp == 1])
      T_perm_non[i]    <- sum(Y_ranks_non   [Zp == 1])
    }
    
    #pvals
    p_oracle <- mean(T_perm_oracle >= obs_T_oracle)
    p_non    <- mean(T_perm_non    >= obs_T_non)
    
    return(list(
      p_oracle = p_oracle,
      p_value  = p_non
    ))
  }
}

getY_aligned <- function(Y_obs, M_obs, Z, block, delta, method = "median_aligned_ranks") {
  K         <- length(unique(block))
  Y_aligned <- numeric(length(Y_obs))
  Y_0       <- Y_obs - Z * delta
  
  for (k in seq_len(K)) {
    idx <- which(block == unique(block)[k])
    if (method == "mean_aligned_ranks") {
      mean_Y_0_k <- mean(Y_0[M_obs == 0 & block == unique(block)[k]], na.rm = TRUE)
    } else if (method == "median_aligned_ranks") {
      mean_Y_0_k <- median(Y_0[M_obs == 0 & block == unique(block)[k]], na.rm = TRUE)
    } else {
      stop("Unknown method for stratified experiment")
    }
    Y_aligned[idx] <- Y_obs[idx] - mean_Y_0_k
  }
  return(Y_aligned)
}

prepare_Z_permutations <- function(Z, block, L=2000) {
  block_ids <- sort(unique(block))
  K <- length(block_ids)
  
  #block indices
  block_indices <- lapply(block_ids, function(b) which(block == b))
  
  #Z permutations
  Z_perms <- vector("list", L)
  for (i in seq_len(L)) {
    Zp <- Z
    for (k in seq_len(K)) {
      idx <- block_indices[[k]]
      Zp[idx] <- sample(Z[idx])
    }
    Z_perms[[i]] <- Zp
  }
  
  list(
    block_indices = block_indices,
    Z_perms = Z_perms
  )
}





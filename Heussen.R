####Heussen

conditional_randomization_test <- function(Y, Z, M, L = 10000, test_stat, 
                                           alternative = "greater") {
  
  n  <- length(Y)
  n_t <- sum(Z)
  n_c <- n - n_t
  m   <- sum(M)
  m_t <- sum(M[Z == 1])
  m_c <- m - m_t
  obs_t <- n_t - m_t
  obs_c <- n_c - m_c
  
  missing_idx  <- which(M == 1)
  observed_idx <- which(M == 0)

  Y_imp <- ifelse(is.na(Y), 0, Y)
  t_obs <- test_stat(Y_imp, Z)
  
  Z_miss_template <- c(rep(1, m_t), rep(0, m_c))
  Z_obs_template  <- c(rep(1, obs_t), rep(0, obs_c))
  
  t_sim <- numeric(L)
  
  for (i in seq_len(L)) {
  
    Z_sim <- Z
    Z_sim[missing_idx]  <- sample(Z_miss_template) #shuffle missing
    Z_sim[observed_idx] <- sample(Z_obs_template)  #suffle observed
    
    t_sim[i] <- test_stat(Y_imp, Z_sim)
  }

  if (alternative == "greater") {
    p_value <- mean(t_sim >= t_obs)
  } else if (alternative == "less") {
    p_value <- mean(t_sim <= t_obs)
  } else {
    p_value <- mean(abs(t_sim) >= abs(t_obs))
  }
  
  list(p_value = p_value)
}

test_stat <- function(Y_imp, Z) {
  mean(Y_imp[Z == 1]) - mean(Y_imp[Z == 0])
}


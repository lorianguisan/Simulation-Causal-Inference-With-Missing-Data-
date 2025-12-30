#load packages in main first

generate_data <- function(n = 100, n_t = 50, distribution ="normal", scenario = "MCAR", 
                          covariates = FALSE, delta = 0, p_miss = 0.05, mu = 0, sd = 1, unif_bound = 2) {    
  
    if (distribution == "normal") {
      Y0 <- rnorm(n, mean = mu, sd = sd)      # N(mu,sd)
      if (covariates == TRUE) X <- rnorm(n, mean = 0, sd = 1)
    } else if (distribution == "uniform") {
      Y0 <- runif(n, min = mu-abs(unif_bound), max = mu+abs(unif_bound))      # U([mu-2,mu+2])
      if (covariates == TRUE) X <- runif(n, min = -2, max = 2)
    } else {
      stop("distribution must be 'normal' or 'uniform'")
    }
  
  
    if (length(delta) == 1) {
      delta <- rep(delta, n)
    } else if (length(delta) != n) {
      stop("delta must be of length 1 or length n")
    }
  
    Y1 <- Y0 + delta # Generate Y1 under H_delta
    
    # CRE with fixed n_t and n_c
    n_c <- n - n_t
    Z <- sample(c(rep(1, n_t), rep(0, n_c)), size = n, replace = FALSE)
    
    Y_true <- Z * Y1 + (1 - Z) * Y0 # True Y_obs
    
    # (M=1 is missing, M=0 is observed !!!!!)
    M0 <- rep(NA, n)
    M1 <- rep(NA, n)
    M_obs <- rep(NA, n)
    
    if (p_miss == 0.05) {
      q <- 0.03  
    }else {
      q <- p_miss/2
    }
    
    switch(scenario,
           
           "MCAR" = { #satisfied for all assumptions
             M_obs <- rbinom(n, 1, p_miss)    #bernoulli
            },
           
           "MAR" = { #satisfies Li A4
             M0 <- rbinom(n, 1, p_miss)     
             M1 <- rbinom(n, 1, p_miss)      
             M_obs <- Z * M1 + (1 - Z) * M0    #M_obs depends on Z, not Y
           },
           
           "threshold_outcome" = { #satisfies Heng A2, easy to control p_miss
             M_obs <- as.integer(
                Y_true <= quantile(Y_true, p_miss/2) |
                Y_true >= quantile(Y_true, 1 - p_miss/2)
             )
           },
           
           "threshold_mixed" = { #satisfies Heng A2, easy to control p_miss
             M_obs <- as.integer(
                Y_true <= quantile(Y_true, p_miss/2) |
                X      >= quantile(X, 1 - p_miss/2)
             )
           },
           
           "threshold_covariates" = { #satisfies Heng A3, easy to control p_miss
             M_obs <- as.integer(
                X <= quantile(X, p_miss/2) |
                X >= quantile(X, 1 - p_miss/2)
             )
           },
           
           "threshold_monotone_pos" = { 
             # Li Ass. 2  Monotonicity ; treatment discourages missingness
             p <- 1-p_miss
             M1 <- as.integer(Y0 >= quantile(Y0, p+q))
             M0 <- as.integer(Y0 >= quantile(Y0, p-q))
             M_obs <- Z * M1 + (1 - Z) * M0 
           },
           
           "threshold_monotone_neg" = {
             # Li Ass. 3  Monotonicity ; treatment encourages missingness
             p <- p_miss
             M1 <- as.integer(Y0 <= quantile(Y0, p+q))
             M0 <- as.integer(Y0 <= quantile(Y0, p-q))
             M_obs <- Z * M1 + (1 - Z) * M0 
           },
           
           "threshold_worst_case" = { 
             M0 <- as.integer(Y0 <= quantile(Y0, p_miss))  # Miss if Y0 is low
             M1 <- as.integer(Y0 >= quantile(Y0, 1-p_miss))  # Miss if Y0 is high 
             M_obs <- Z * M1 + (1 - Z) * M0
           },
           
           "threshold_best_case" = {
             M0 <- as.integer(Y0 >= quantile(Y0, 1-p_miss))  # Miss if Y0 is high
             M1 <- as.integer(Y0 <= quantile(Y0, p_miss))  # Miss if Y0 is low (inverted)
             M_obs <- Z * M1 + (1 - Z) * M0 
           },
           
           stop("Invalid scenario.")
    )
    
    # Observed Outcome (NA if Mi=1)
    Y_obs <- Y_true
    Y_obs[M_obs == 1] <- NA
    
    # Impute potential outcome missingness
    M0[Z == 0] <- M_obs[Z == 0] #unobserved remain NA
    M1[Z == 1] <- M_obs[Z == 1] #unobserved remain NA

    if (covariates == TRUE) {
      return(data.frame(X, Z, Y0, Y1, Y_true, M0, M1, M_obs, Y_obs))
    } else {
      return(data.frame(Z, Y0, Y1, Y_true, M0, M1, M_obs, Y_obs))
    }
  }

generate_stratified_data <- function(n_blocks=c(20,20,20,20,20), n_t_blocks = c(10,10,10,10,10), mu_blocks = c(0,0,0,0,0),
                                     sd_blocks = c(1,1,1,1,1), distribution_blocks ="normal", delta_blocks = c(0),
                                     scenario_blocks = "Li_A1", p_miss_blocks = c(0.05), unif_blocks = c(2))
  {
  K = length(n_blocks)
  
  expand_to_K <- function(x) {
    if (length(x) == 1) rep(x, K)
    else if (length(x) == K) x
    else stop(paste("Parameter has wrong length:", paste(x, collapse=",")))
  } #checks length
  
  n_blocks        <- expand_to_K(n_blocks)
  n_t_blocks      <- expand_to_K(n_t_blocks)
  mu_blocks       <- expand_to_K(mu_blocks)
  sd_blocks       <- expand_to_K(sd_blocks)
  delta_blocks    <- expand_to_K(delta_blocks)
  distribution_blocks <- expand_to_K(distribution_blocks)
  scenario_blocks <- expand_to_K(scenario_blocks)
  p_miss_blocks   <- expand_to_K(p_miss_blocks)
  unif_blocks     <- expand_to_K(unif_blocks)
  
  if (any(n_blocks <= 0)) stop("All n_blocks must be > 0.")
  if (any(n_t_blocks < 0)) stop("All n_t_blocks must be >= 0.")
  if (any(n_t_blocks > n_blocks)) stop("n_t_blocks[i] must be less than n_blocks[i].")
  if (sum(n_t_blocks) >= sum(n_blocks)) stop("n_t must be less than n.")
  
  
  strata_list <- vector("list", K)
  
  for (i in 1:K) {
    strata_list[[i]] <- generate_data(
      n = n_blocks[i], 
      n_t = n_t_blocks[i],
      distribution = distribution_blocks[i],
      scenario = scenario_blocks[i],
      covariates = FALSE,
      delta = delta_blocks[i],
      p_miss = p_miss_blocks[i],
      mu = mu_blocks[i],
      sd = sd_blocks[i],
      unif_bound = unif_blocks[i]
      
    ) %>% mutate(block = i)  # block number
  }
  
  data_stratified <- bind_rows(strata_list)
  return(data_stratified)
}
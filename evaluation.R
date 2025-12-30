#load packages in main first

compute_pvalue <- function(Y, Z, M, X, delta, assumption) {

  if (assumption %in% c("Li_A1", "Li_A2", "Li_A3", "Li_A4")) {
    if (length(delta) == 1) delta <- rep(delta, length(Y))
    li_test <- worst_case_randomization_test(Y, Z, M, delta, assumption)
    return(li_test$p_value)
    
  } else
  
  if (assumption %in% c("Heng")) {
    heng_test <- iArt.test(Z = Z, X = X, Y = matrix(Y, ncol=1),
                         L = 200, covariate_adjustment = FALSE,
                         verbose = FALSE)
    return(heng_test$p_value[1])
    
  } else 
    
  if (assumption %in% c("Heussen")) {
    heussen_test <- conditional_randomization_test(Y, Z, M, L = 1000, test_stat)
    return(heussen_test$p_value)
  }
  
  else stop("Invalid assumption")
}

run_sim <- function(
    data_scenarios = c(
      "MCAR", "MAR",
      "threshold_outcome", "threshold_mixed", "threshold_covariates",
      "threshold_monotone_pos", "threshold_monotone_neg",
      "threshold_worst_case", "threshold_best_case"
    ),
    assumptions = c("Li_A1", "Li_A2", "Li_A3", "Li_A4", "Heng", "Heussen"),
    n_configs = list(list(n = 1000, n_t = 500)),
    p_miss_vals = seq(0.05, 0.65, by = 0.1),
    delta_vals = seq(0, 3, by = 0.1),
    reps = 500,
    distribution = "normal",
    compute_H_delta = FALSE
) {

  base_grid <- expand_grid(
    data_scenario = data_scenarios,
    p_miss = p_miss_vals,
    delta = delta_vals
  )

  n_configs_df <- do.call(rbind, lapply(n_configs, as.data.frame))

  sim_grid <- base_grid %>%
    tidyr::crossing(n_configs_df) %>%
    tidyr::crossing(assumption = assumptions)
  
  cat("Running simulation...\n")
  cat("Total unique scenarios:", nrow(sim_grid), "\n")
  cat("Replications per scenario:", reps, "\n")
  cat("Total runs:", nrow(sim_grid) * reps, "\n\n")
  
  results_list <- vector("list", nrow(sim_grid))
  
  for (i in seq_len(nrow(sim_grid))) {
    sc_data <- sim_grid$data_scenario[i]
    sc_assum <- sim_grid$assumption[i]
    n <- sim_grid$n[i]
    n_t <- sim_grid$n_t[i]
    pm <- sim_grid$p_miss[i]
    del <- sim_grid$delta[i]
    
    if (sc_assum %in% c("Heng", "Heussen")) cat("...running scenario", i, "of", nrow(sim_grid), "\n")
    else if (i %% 10 == 0) cat("...running scenario", i, "of", nrow(sim_grid), "\n")
    
    # Decide if covariates are needed
    cov_flag <- sc_assum %in% c("Heng")
    
    # Storage for p-values
    li_flag <- sc_assum %in% c("Li_A1", "Li_A2", "Li_A3", "Li_A4")
    
    if (compute_H_delta && li_flag) {
      p_value_results <- matrix(NA_real_, nrow = reps, ncol = 3)
      colnames(p_value_results) <- c("pval_H_zero", "pval_H_delta", "pval_true_data")
      
    } else if (!compute_H_delta && li_flag) {
      p_value_results <- matrix(NA_real_, nrow = reps, ncol = 2)
      colnames(p_value_results) <- c("pval_H_zero", "pval_true_data")
      
    } else {
      p_value_results <- matrix(NA_real_, nrow = reps, ncol = 1)
      colnames(p_value_results) <- "pval_H_zero"
    }
    
    for (r in 1:reps) {
      set.seed((i * 1000) + r) #for reproducibility
      
      # Generate data
      dat <- generate_data(
        n          = n,
        n_t        = n_t,
        scenario   = sc_data,
        p_miss     = pm,
        delta      = del,
        covariates = cov_flag,
        distribution = distribution
      )
      
      # Construct X as a matrix if covariates exist
      X_mat <- if (cov_flag && "X" %in% names(dat)) cbind(dat$X) else NULL
      
      # Compute p-value for H0
      pval_h0 <- compute_pvalue(
        Y          = dat$Y_obs,
        Z          = dat$Z,
        M          = dat$M_obs,
        X          = X_mat,
        delta      = 0,
        assumption = sc_assum
      )
      
      # Compute p-value for H_delta only if we want
      if (compute_H_delta) {
        if (cov_flag) { # Heng's test can't handle delta != 0
          pval_hdelta <- pval_h0
        } else {
          pval_hdelta <- compute_pvalue(
            Y          = dat$Y_obs,
            Z          = dat$Z,
            M          = dat$M_obs,
            X          = X_mat,
            delta      = del,
            assumption = sc_assum
          )
        }
      }
      
      # Oracle p-value 
      if (li_flag) {
        true_test <- wilcox.test(
          dat$Y_true[dat$Z == 1],
          dat$Y_true[dat$Z == 0],
          alternative = "greater"
        )
      }
      
      # Store p-values
      if (li_flag) {
        if (compute_H_delta) {
          p_value_results[r, ] <- c(pval_h0, pval_hdelta, true_test$p.value)
        } else {
          p_value_results[r, ] <- c(pval_h0, true_test$p.value)
        }
      } else {
          p_value_results[r, ] <- pval_h0
      }
    
    results_list[[i]] <- as.data.frame(p_value_results)
    }
  }
  
  cat("...simulation complete.\n")

  all_p_values_df <- bind_rows(results_list)
  
  sim_grid_expanded <- sim_grid[rep(seq_len(nrow(sim_grid)), each = reps), ]
  
  final_results_df <- bind_cols(sim_grid_expanded, all_p_values_df)
  
  return(final_results_df)
}

compute_rejection_rates <- function(sim_results, alpha = 0.05) {
  
  cat("Calculating rejection rates...\n")
  
  has_H_delta <- "pval_H_delta" %in% colnames(sim_results)
  has_oracle <- "pval_true_data" %in% colnames(sim_results)
  
  summary_expr <- list(
    rej_rate_H0        = ~ mean(pval_H_zero    < alpha, na.rm = TRUE),
    rej_rate_true_data = ~ mean(pval_true_data < alpha, na.rm = TRUE)
  )
  if (has_H_delta) {
    summary_expr$rej_rate_H_delta <- ~ mean(pval_H_delta < alpha, na.rm = TRUE)
  }
  rej_rates_wide <- sim_results %>%
    group_by(data_scenario, assumption, n, n_t, p_miss, delta) %>%
    summarise(
      across(any_of(c(
        "pval_H_zero",
        if (has_H_delta) "pval_H_delta",
        if (has_oracle)  "pval_true_data"
      )),
             ~ mean(.x < alpha, na.rm = TRUE),
             .names = "rej_rate_{.col}"),
      n_reps = n(),
      .groups = "drop"
    )
  
  rej_rates_long <- rej_rates_wide %>%
    tidyr::pivot_longer(
      cols = starts_with("rej_rate"),
      names_to = "curve_type",
      values_to = "rej_rate"
    ) %>%
    mutate(
      curve_type = dplyr::recode(
        curve_type,
        rej_rate_pval_H_zero  = "H_0 test",
        !!!if (has_H_delta) list(rej_rate_pval_H_delta   = "H_delta test"),
        !!!if (has_oracle)  list(rej_rate_pval_true_data = "Oracle (true data)")
      )
      
    )
  
  return(rej_rates_long)
}


run_sim_SRE <- function(
    data_scenarios   = c(
      "MCAR", "MAR", "threshold_outcome",
      "threshold_monotone_pos", "threshold_monotone_neg",
      "threshold_worst_case", "threshold_best_case"
    ),
    assumptions      = c("Li_A1", "Li_A2", "Li_A3"),
    n_configs = list(list(
      n_blocks   = c(20, 20, 20, 20, 20),
      n_t_blocks = c(10, 10, 10, 10, 10),
      mu_blocks  = c(0, 1, -1, 2, -2),
      sd_blocks  = c(1, 1, 1, 1, 1)
    )),
    p_miss_vals      = seq(0.05, 0.65, by = 0.1),
    delta_vals       = seq(0, 3, by = 0.1),
    reps             = 500,
    distribution     = "normal",
    method           = "median_aligned_ranks",
    compute_H_delta  = FALSE
) {
  
  # base grid over scenario-level parameters
  base_grid <- expand_grid(
    data_scenario = data_scenarios,
    p_miss        = p_miss_vals,
    delta         = delta_vals
  )
  
  # convert n_configs list into a dataframe
  n_configs_df <- bind_rows(
    lapply(n_configs, function(cfg) {
      tibble(
        n_blocks   = list(cfg$n_blocks),
        n_t_blocks = list(cfg$n_t_blocks),
        mu_blocks  = list(cfg$mu_blocks),
        sd_blocks  = list(cfg$sd_blocks)
      )
    })
  )
  
  # full cross-product grid (scenario × configs × assumptions)
  sim_grid <- base_grid %>%
    tidyr::crossing(n_configs_df) %>%
    tidyr::crossing(tibble(assumption = assumptions))
  
  cat("Running SRE simulation…\n")
  cat("Total unique scenarios:", nrow(sim_grid), "\n")
  cat("Replications per scenario:", reps, "\n")
  cat("Total runs:", nrow(sim_grid) * reps, "\n\n")
  
  results_list <- vector("list", nrow(sim_grid))
  
  
  #loop over SRE scenarios
  for (i in seq_len(nrow(sim_grid))) {

    cat("...running SRE scenario", i, "of", nrow(sim_grid), "\n")
    
    # scenario-level values
    sc_data  <- sim_grid$data_scenario[i]
    sc_assum <- sim_grid$assumption[i]
    pm       <- sim_grid$p_miss[i]
    del      <- sim_grid$delta[i]
    
    # block-wise parameters (LIST COLS)
    n_blocks   <- sim_grid$n_blocks[[i]]
    n_t_blocks <- sim_grid$n_t_blocks[[i]]
    mu_blocks  <- sim_grid$mu_blocks[[i]]
    sd_blocks  <- sim_grid$sd_blocks[[i]]
    
    # prepare storage
    if (compute_H_delta) {
      pvals <- matrix(NA_real_, reps, 3)
      colnames(pvals) <- c("pval_H_zero", "pval_H_delta", "pval_true_data")
    } else {
      pvals <- matrix(NA_real_, reps, 2)
      colnames(pvals) <- c("pval_H_zero", "pval_true_data")
    }
    
    #reps
    
    for (r in 1:reps) {
      
      # delta vector per block
      delta_blocks <- rep(del, length(n_blocks))
      
      # generate stratified data 
      dat <- generate_stratified_data(
        n_blocks            = n_blocks,
        n_t_blocks          = n_t_blocks,
        mu_blocks           = mu_blocks,
        sd_blocks           = sd_blocks,
        distribution_blocks = distribution,
        scenario_blocks     = sc_data,
        delta_blocks        = delta_blocks,
        p_miss_blocks       = rep(pm, length(n_blocks))
      )
      
      # H0 test (delta = 0)
      h0_test <- worst_case_stratified_randomization_test(
        data         = dat,
        delta_blocks = rep(0, length(n_blocks)),
        assumption   = sc_assum,
        method       = method
      )
      p_h0 <- h0_test$p_value
      
      # oracle
      p_true <- h0_test$p_oracle
      
      
      if (compute_H_delta) {  # H_delta test
      p_hdelta <- worst_case_stratified_randomization_test(
        data         = dat,
        delta_blocks = delta_blocks,
        assumption   = sc_assum,
        method       = method
      )$p_value
      }
      
      if (compute_H_delta) {
        pvals[r, ] <- c(p_h0, p_hdelta, p_true)
      } else {
        pvals[r, ] <- c(p_h0, p_true)
      }
    }
    
    results_list[[i]] <- as.data.frame(pvals)
  }
  
  
  cat("\nSRE simulation complete.\n")
  
  #dataframe
  
  sim_grid_expanded <- sim_grid[rep(1:nrow(sim_grid), each = reps), ]
  final_df <- bind_cols(sim_grid_expanded, bind_rows(results_list))
  
  return(final_df)
}

compute_rejection_rates_SRE <- function(sim_results, alpha = 0.05) {
  
  cat("Calculating rejection rates...\n")
  
  has_H_delta <- "pval_H_delta" %in% colnames(sim_results)
  
  summary_expr <- list(
    rej_rate_H0        = ~ mean(pval_H_zero    < alpha, na.rm = TRUE),
    rej_rate_true_data = ~ mean(pval_true_data < alpha, na.rm = TRUE)
  )
  
  if (has_H_delta) {
    summary_expr$rej_rate_H_delta <- ~ mean(pval_H_delta < alpha, na.rm = TRUE)
  }
  
  rej_rates_wide <- sim_results %>%
    group_by(
      data_scenario, 
      assumption, 
      p_miss, 
      delta,
      n_blocks,
      n_t_blocks
    ) %>%
    summarise(
      across(
        any_of(c("pval_H_zero", "pval_H_delta", "pval_true_data")),
        ~ mean(.x < alpha, na.rm = TRUE),
        .names = "rej_rate_{.col}"
      ),
      n_reps = n(),
      .groups = "drop"
    )
  
  
  # wide to long format for plotting
  rej_rates_long <- rej_rates_wide %>%
    tidyr::pivot_longer(
      cols = starts_with("rej_rate"),
      names_to = "curve_type",
      values_to = "rej_rate"
    ) %>%
    mutate(
      curve_type = dplyr::recode(
        curve_type,
        rej_rate_pval_H_zero        = "H_0 test",
        rej_rate_pval_H_delta   = "H_delta test",
        rej_rate_pval_true_data = "Oracle (true data)"
      )
    )
  
  return(rej_rates_long)
}

compute_pvalue_CRE_on_SRE <- function(dat, delta, assumption) {
  
  # ignors block effect and treats it as a CRE
  worst_case_randomization_test(
    Y     = dat$Y_obs,
    Z     = dat$Z,
    M     = dat$M_obs,
    delta = delta,
    assumption = assumption
  )$p_value
}

run_sim_SRE_as_CRE <- function(
    data_scenarios = c("threshold_monotone_pos"),
    assumption     = "Li_A2",
    n_blocks       = c(20, 20, 20, 20, 20),
    n_t_blocks     = c(10, 10, 10, 10, 10),
    mu_blocks      = c(0, 1, -1, 2, -2),
    sd_blocks      = c(1, 1, 1, 1, 1),
    p_miss_vals    = seq(0.05, 0.25, by = 0.1),
    delta_vals     = seq(0, 2.5, by = 0.1),
    reps           = 1000,
    distribution   = "normal"
) {
  
  grid <- expand_grid(
    data_scenario = data_scenarios,
    p_miss        = p_miss_vals,
    delta         = delta_vals
  )
  
  results <- vector("list", nrow(grid))
  
  for (i in seq_len(nrow(grid))) {
    
    cat("Running scenario", i, "of", nrow(grid), "\n")
    
    sc  <- grid$data_scenario[i]
    pm  <- grid$p_miss[i]
    del <- grid$delta[i]
    
    pvals <- numeric(reps)
    
    for (r in 1:reps) {
      
      dat <- generate_stratified_data(
        n_blocks            = n_blocks,
        n_t_blocks          = n_t_blocks,
        mu_blocks           = mu_blocks,
        sd_blocks           = sd_blocks,
        distribution_blocks = distribution,
        scenario_blocks     = sc,
        delta_blocks        = rep(del, length(n_blocks)),
        p_miss_blocks       = rep(pm, length(n_blocks))
      )
      
      pvals[r] <- compute_pvalue_CRE_on_SRE(
        dat,
        delta      = 0,        # H0
        assumption = assumption
      )
    }
    
    results[[i]] <- tibble(
      data_scenario = sc,
      p_miss        = pm,
      delta         = del,
      p_value       = pvals
    )
  }
  
  bind_rows(results)
}




  


  

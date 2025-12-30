# Load packages

library(dplyr)
library(ggplot2)
library(readr)
library(iArt)
library(mice) #for heng ? not sure if mandatory (used in iArt package)
library(tidyr)
library(tibble)
library(purrr)
library(coin) 

#######################################################

# Example 1 : Li A1 vs. Li A2

sim_results_li_v_li <- run_sim(
  data_scenarios   = c(
    "MCAR", "threshold_monotone_pos",
    "threshold_worst_case", "threshold_best_case"
  ),
  assumptions      = c("Li_A1", "Li_A2"),
  n_configs      = list(
    list(n = 100, n_t = 50)
  ),
  p_miss_vals = seq(0.05, 0.35, by = 0.1), 
  delta_vals  = seq(0, 2.5, by = 0.05), 
  reps        = 1000,
  distribution     = "normal"
)
rej_rates_li_v_li <- compute_rejection_rates(sim_results_li_v_li,alpha = 0.05)

rej_rates_li_v_li <- rej_rates_li_v_li %>%
  filter(curve_type %in% c("H_0 test"))

rej_rates_li_v_li <- rej_rates_li_v_li %>%
  filter(data_scenario %in% c("threshold_monotone_pos"))

rej_rates_li_v_li <- rej_rates_li_v_li %>%
  filter(p_miss %in% c(0.25, 0.35))

plot_rejection_rates(
  data_long = rej_rates_li_v_li,
  x_var = "delta",
  y_var = "rej_rate",
  color_var = "assumption",        # different assumptions are different color curves
  facet_var = "p_miss",     # one plot per data generation
  fixed_params = list(n = 100, data_scenario = "threshold_monotone_pos")
)

write.csv(sim_results_li_v_li, "sim_results_li1_v_li2.csv", row.names = FALSE)

#######################################################

# Example 2 :comparing different % in missingness under Li A1

sim_results_Li1 <- run_sim(
  data_scenarios   = c(
    "threshold_worst_case", "threshold_best_case"
  ),
  assumptions      = c("Li_A1"),
  n_configs      = list(
    list(n = 100, n_t = 50)
  ),
  p_miss_vals      = seq(0.05, 0.35, by = 0.1), 
  delta_vals       = seq(0, 2.5, by = 0.05), 
  reps             = 1000,
  distribution     = "normal"
)
rej_rates_li1 <- compute_rejection_rates(sim_results_Li1, alpha = 0.05)


plot_rejection_rates(
  data_long = rej_rates_li1,
  x_var = "delta",
  y_var = "rej_rate",
  color_var = "p_miss",        
  facet_var = "data_scenario",     
  fixed_params = list(n = 100, assumption="Li_A1")
)

write.csv(sim_results_Li1, "sim_results_Li1.csv", row.names = FALSE)

#######################################################

# Example 3 : Heng (took ~~14 hours)

sim_results_heng <- run_sim(
  data_scenarios = c("MCAR", "threshold_outcome"),
  assumptions    = c("Heng"),
  n_configs      = list(list(n = 100, n_t = 50)),
  p_miss_vals    = seq(0.05, 0.35, by = 0.1),
  delta_vals     = seq(0, 1, by = 0.1),
  reps           = 100,   # small since Heng's method already needs a lot of reps to compute p_vals
  distribution     = "normal"
)

rej_rates_heng <- compute_rejection_rates(sim_results_heng,alpha = 0.05)

rej_H0_heng <- rej_rates_heng %>%
  filter(curve_type %in% c("H_0 test"))

plot_rejection_rates(
  data_long = rej_H0_heng,
  x_var = "delta",
  y_var = "rej_rate",
  color_var = "p_miss",        # different assumptions are different color curves
  facet_var = "data_scenario",     # one plot per data generation
  fixed_params = list(n = 100, assumption="Heng")
)

write.csv(sim_results_heng, "sim_results_heng.csv", row.names = FALSE)


  
  #######################################################
  
  # Example 4 : Heussen ~ took < 1hour
  sim_results_heussen <- run_sim(
    data_scenarios   = c("MCAR", "threshold_outcome"),
    assumptions      = c("Heussen"),
    n_configs      = list(
      list(n = 100, n_t = 50)
    ),
    p_miss_vals = seq(0.05, 0.35, by = 0.1), 
    delta_vals  = seq(0, 2.5, by = 0.1), 
    reps        = 200,
    distribution     = "normal"
  )
  rej_rates_heussen <- compute_rejection_rates(sim_results_heussen, alpha = 0.05)
  
  #rej_H0 <- rej_rates_example %>% filter(curve_type == "H_0 test", n == 1000)
  
  plot_rejection_rates(
    data_long = rej_rates_heussen,
    x_var = "delta",
    y_var = "rej_rate",
    color_var = "p_miss",        # different assumptions are different color curves
    facet_var = "data_scenario",     # one plot per data generation
    fixed_params = list(n = 100, assumption="Heussen")
  )
  
  write.csv(sim_results_heussen, "sim_results_heussen.csv", row.names = FALSE)
  
  #######################################################
  
  # Example 5 : Comparing Heng and Heussen :

  p_miss_vals <- seq(0.05, 0.35, by = 0.1)
  for (p_miss_val in p_miss_vals) {
  scenarios_to_plot <- c("MCAR", "threshold_outcome")
  delta_max <- 1  # maximum delta to include
  
  rej_rates_heng <- rej_rates_heng %>%
    filter(curve_type %in% c("H_0 test"))
  rej_all <- bind_rows(rej_rates_heussen, rej_rates_heng)
  
  rej_all <- rej_all %>%
    mutate(method = case_when(
      assumption %in% c("Heussen") ~ "Heussen",
      assumption %in% c("Heng_A2") ~ "Heng",
      TRUE ~ assumption
    ))
  
  rej_plot <- rej_all %>%
    filter(
      p_miss == p_miss_val,
      data_scenario %in% scenarios_to_plot,
      delta <= delta_max
    )
  
  plot_rejection_rates(
    data_long = rej_plot,
    x_var = "delta",
    y_var = "rej_rate",
    color_var = "method",   # color by method (Heussen vs Heng)
    facet_var = "data_scenario",  # optional: separate MAR vs threshold_outcome
    fixed_params = list(p_miss = p_miss_val)
  )
  }
  
  
  #######################################################
  
  # Example 6 : SRE_aligned
  
  sim_results_SRE_example <- run_sim_SRE(
    
    data_scenarios = c(
      "MAR", "threshold_monotone_pos"),
    assumptions = c("Li_A2"),
    n_configs = list(list(
      n_blocks   = c(20, 20, 20, 20, 20),
      n_t_blocks = c(10, 10, 10, 10, 10),
      mu_blocks  = c(0, 1, -1, 2, -2),
      sd_blocks  = c(1, 1, 1, 1, 1)
    )),
    p_miss_vals = seq(0.05, 0.25, by = 0.1),
    delta_vals = seq(0, 2.5, by = 0.1),
    reps = 1000,
    distribution = "normal",
    compute_H_delta  = FALSE
  )
  
  
  rej_rates_SRE_example <- compute_rejection_rates_SRE(
    sim_results_SRE_example,
    alpha = 0.05
  )
  
  plot_rejection_rates(
    data_long   = rej_rates_SRE_example,
    x_var       = "delta",
    y_var       = "rej_rate",
    color_var   = "p_miss",            
    facet_var   = "data_scenario",  
    fixed_params = list(
      n_blocks = c(20,20,20,20,20),
      assumption = "Li_A2"
    )
  )
  sim_results_SRE_example_clean <- tidyr::unnest(sim_results_SRE_example)
  write.csv(sim_results_SRE_example_clean, "sim_results_SRE.csv", row.names = FALSE)
  
  #######################################################
  
  
  # Example 7 : Assumptions not satisfied ?
  sim_results_example7 <- run_sim(
    data_scenarios   = c(
      "MAR", "threshold_worst_case",
      "threshold_best_case", "threshold_monotone_neg"
    ),
    assumptions      = c("Li_A2"),
    n_configs        = list(
      list(n = 100, n_t = 50)
    ),
    p_miss_vals      = seq(0.05, 0.35, by = 0.1), 
    delta_vals       = seq(0, 2.5, by = 0.05), 
    reps             = 1000,
    distribution     = "normal",
    compute_H_delta  = FALSE
  )
  rej_rates_example7 <- compute_rejection_rates(sim_results_example7,alpha = 0.05)
  
  rej_rates_example7 <- rej_rates_example7 %>%
    filter(data_scenario %in% c("MAR", "threshold_monotone_neg"))
  
  plot_rejection_rates(
    data_long = rej_rates_example7,
    x_var = "delta",
    y_var = "rej_rate",
    color_var = "p_miss",       
    facet_var = "data_scenario",    
    fixed_params = list(n = 100, assumption="Li_A2")
  )
  
  rej_rates_example7 <- compute_rejection_rates(sim_results_example7,alpha = 0.05)
  
  rej_rates_example7 <- rej_rates_example7 %>%
    filter(data_scenario %in% c("threshold_best_case", "threshold_worst_case"))
  
  plot_rejection_rates(
    data_long = rej_rates_example7,
    x_var = "delta",
    y_var = "rej_rate",
    color_var = "p_miss",       
    facet_var = "data_scenario",    
    fixed_params = list(n = 100, assumption="Li_A2")
  )
  
  #######################################################
  
  # Example 8 : Van elteren method for SRE 
  
  sim_results_SRE_van_elteren <- run_sim_SRE(
    
    data_scenarios = c(
      "MAR", "threshold_monotone_pos"),
    assumptions = c("Li_A2"),
    n_configs = list(list(
      n_blocks   = c(20, 20, 20, 20, 20),
      n_t_blocks = c(10, 10, 10, 10, 10),
      mu_blocks  = c(0, 1, -1, 2, -2),
      sd_blocks  = c(1, 1, 1, 1, 1)
    )),
    p_miss_vals = seq(0.05, 0.25, by = 0.1),
    delta_vals = seq(0, 2.5, by = 0.1),
    reps = 1000,
    distribution = "normal",
    method = "combined_ranks",
    compute_H_delta  = FALSE
  )
  
  
  rej_rates_SRE_van_elteren <- compute_rejection_rates_SRE(
    sim_results_SRE_van_elteren,
    alpha = 0.05
  )
  
  plot_rejection_rates(
    data_long   = rej_rates_SRE_van_elteren,
    x_var       = "delta",
    y_var       = "rej_rate",
    color_var   = "p_miss",            
    facet_var   = "data_scenario",  
    fixed_params = list(
      n_blocks = c(20, 20, 20, 20, 20),
      assumption = "Li_A2"
    )
  )

  
  #######################################################
  
  # Example 9 : Visualization of the randomization distribution of the 
  #             Wilcoxon rank-sum statistic with n=100, n_t=50

  
  n_t <- 50
  n_c <- 50
  n_sims <- 100000
  mu_W <- n_t * (n_t + n_c + 1) / 2        # theoretical mean (=expectation) (=2525)
   
  null_dist_U <- rwilcox(n_sims, m = n_t, n = n_c)   # rwilcox actually doesn't use the usual wilcoxon statistic
  null_dist_W <- null_dist_U + (n_t * (n_t + 1) / 2) ## to go back to wilcoxon stat.
  plot_data <- data.frame(W_stat = null_dist_W)
  
  p <- ggplot(plot_data, aes(x = W_stat)) +
    
    geom_histogram(
      fill = "gray85", 
      color = "black", 
      binwidth = 10,  
      size = 0.3
    ) +

    geom_vline(xintercept = mu_W, linetype = "dashed", color = "red", linewidth = 1) +
    
    scale_x_continuous(breaks = seq(min(plot_data$W_stat), max(plot_data$W_stat), by = 100)) +
    
    labs(
      title = "Randomization Distribution of Wilcoxon Statistic",
      subtitle = paste0("CRE with N=100 (nt=", n_t, ", nc=", n_c, ") | 100,000 Permutations"),
      x = "Wilcoxon Rank Sum Statistic (W)",
      y = "Frequency (Count)"
    ) +

    theme_bw() +
    theme(
      plot.subtitle   = element_text(hjust = 0.5, face = "italic", size = 14),
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title      = element_text(size = 14),
      axis.text       = element_text(size = 12),
      panel.grid.minor = element_blank() 
    )
  
  print(p)
  
  
  #######################################################
  
  #Example 10 : Comparison between SRE methods vs SRE generated data under CRE method
  
  
  rej_rates_SRE_aligned <- rej_rates_SRE_example %>%
    filter(data_scenario %in% c("threshold_monotone_pos"))
  
  rej_rates_SRE_vanelteren <- rej_rates_SRE_van_elteren %>%
    filter(data_scenario %in% c("threshold_monotone_pos"))
  
  sim_results_SRE_as_CRE <- run_sim_SRE_as_CRE(
    p_miss_vals = seq(0.05, 0.25, by = 0.1),
    delta_vals  = seq(0, 2.5, by = 0.1),
    reps        = 1000
  )
  
  rej_rates_SRE_incorrect <- sim_results_SRE_as_CRE %>%
    group_by(p_miss, delta) %>%
    summarise(
      rej_rate = mean(p_value < 0.05),
      .groups = "drop"
    ) %>%
    mutate(method = "CRE on SRE data")
  
  rej_SRE_aligned <- rej_rates_SRE_correct %>%
    filter(curve_type == "H_0 test") %>%   # no oracle
    mutate(
      curve_type = "SRE_aligned",
      p_miss = round(p_miss, 3)           
    ) %>%
    select(p_miss, delta, rej_rate, curve_type)
  
  rej_SRE_combined <- rej_rates_SRE_vanelteren %>%
    filter(curve_type == "H_0 test") %>%   # no oracle
    mutate(
      curve_type = "SRE_combined",
      p_miss = round(p_miss, 3)           
    ) %>%
    select(p_miss, delta, rej_rate, curve_type)
  
  rej_CRE <- rej_rates_SRE_incorrect %>%
    mutate(
      curve_type = "CRE",
      p_miss = round(p_miss, 3)            # float issues
    ) %>%
    select(p_miss, delta, rej_rate, curve_type)
  
  rej_all <- bind_rows(rej_SRE_aligned, rej_SRE_combined, rej_CRE)
  rej_all <- rej_all %>% filter(p_miss %in% c(0.05, 0.15))
  
  plot_rejection_rates(
    data_long   = rej_all,
    x_var       = "delta",
    y_var       = "rej_rate",
    color_var   = "curve_type",
    facet_var   = "p_miss",
    legend_title = "Design",
    fixed_params = list(n = 100, assumption="Li_A2")
  )
  
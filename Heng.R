if (!requireNamespace("mice", quietly = TRUE)) {
  install.packages("mice")
}
library(mice)

if (!requireNamespace("missForest", quietly = TRUE)) {
  install.packages("missForest")
}
library(missForest)


T <- function(z, y) {
  #Calculates the test statistic T Wilcoxon rank-sum test.
  
  n <- length(z)
  t <- 0
  my_list <- data.frame(z, y)
  sorted_list <- my_list[order(my_list$y),]
  
  for (i in 1:n) {
    t <- t + sorted_list$z[i] * (i + 1)
  }
  
  return(t)
}

getY <- function(G, Z, X, Y, covariate_adjustment = FALSE) {
  #Imputes the missing values in Y using the imputation model G.
  
  # Combine Z, X, Y into a single data frame
  df_Z <- data.frame(cbind(Z, X, Y))
  
  # Impute the missing values in the combined data frame
  imputed_data <- G(df_Z)
  
  # Extract imputed Y values
  lenY <- ncol(Y)
  indexY <- ncol(Z) + ncol(X) + 1 # Assuming Z and X are not NULL
  Y_head <- imputed_data[, indexY:(indexY + lenY - 1)]
  Y_head <- data.frame(Y_head)
  X <- data.frame(imputed_data[, (ncol(Z) + 1):(ncol(Z) + ncol(X))])
  
  if (covariate_adjustment) {
    # Perform linear regression
    column_names <- names(Y_head)

    # Construct the formula
    formula_str <- paste("cbind(", paste(column_names, collapse = ", "), ") ~ .")
    model_formula <- as.formula(formula_str)

    # Fit the model
    model <- lm(model_formula, data = cbind(Y_head, X))

    # Calculate residuals (Y_head - Y)
    Y_head_adj <- predict(model,newdata = X)
    Y_head  = Y_head - Y_head_adj
  }
  
  Y_head <- data.frame(Y_head)
  return(Y_head)
}

split_data <- function(y, z, M) {
  #Splits the data into missing and non-missing parts based on the missingness indicator matrix M.

  # Convert M to a logical vector to identify missing values
  missing_indices <- as.logical(M)

  non_missing_indices <- !missing_indices

  # Split y based on missing and non-missing indices
  y_missing <- y[missing_indices, ]
  y_non_missing <- y[non_missing_indices, ]

  # Split z in the same way
  z_missing <- z[missing_indices, ]
  z_non_missing <- z[non_missing_indices, ]

  return(list(y_missing = y_missing, y_non_missing = y_non_missing, 
              z_missing = z_missing, z_non_missing = z_non_missing))
}


getT <- function(y, z, lenY, M) {
  #Calculates the test statistic T in missing part and non missing part for each column based on the z,imputed y

  t <- numeric(lenY)
  
  for (i in 1:lenY) {
    # Split the data into missing and non-missing parts
    split_result <- split_data(y[, i, drop = FALSE], z, M[, i, drop = FALSE])
    
    # Calculate T for missing and non-missing parts
    t_missing <- if (length(split_result$y_missing) > 0 && length(split_result$z_missing) > 0) {
      T(split_result$z_missing, split_result$y_missing)
    } else {
      0
    }
    
    t_non_missing <- if (length(split_result$y_non_missing) > 0 && length(split_result$z_non_missing) > 0) {
      T(split_result$z_non_missing, split_result$y_non_missing)
    } else {
      0
    }
    
    # Sum the T values for both parts
    t[i] <- t_missing + t_non_missing
  }
  
  return(t)
}

getZsimTemplates <- function(Z_sorted, S) {
  #Creates a list of templates for each stratum in S, the template is then used to simulate Z generated from the same distribution as Z.

  Z_sim_templates <- list()
  
  unique_strata <- unique(S)
  
  for (stratum in unique_strata) {
    strata_indices <- which(S == stratum)
    strata_Z <- Z_sorted[strata_indices]
    p <- mean(strata_Z)
    strata_size <- length(strata_indices)
    
    # Create a template for the stratum
    Z_sim_template <- c(rep(0.0, strata_size * (1 - p)), rep(1.0, strata_size * p))
    
    Z_sim_templates[[as.character(stratum)]] <- Z_sim_template
  }
  return(Z_sim_templates)
}

getZsim <- function(Z_sim_templates) {
  Z_sim <- numeric()
  
  for (Z_sim_template in Z_sim_templates) {
    # Shuffle the current template
    strata_Z_sim <- sample(Z_sim_template)
    
    # Append the shuffled template to Z_sim
    Z_sim <- c(Z_sim, strata_Z_sim)
  }
  
  return(Z_sim)
}

check_param <- function(Z, X, Y, S, G, L, verbose, covariate_adjustment, alpha, alternative, random_state) {
  #Checks the parameters passed to the iartest function.
  
  # Check if Z, X, Y are matrices or data frames
  if (!is.matrix(X) && !is.data.frame(X)) stop("X must be a matrix or data frame")
  if (!is.matrix(Y) && !is.data.frame(Y)) stop("Y must be a matrix or data frame")
  
  # Check Z: must contain only 0, 1
  if (!all(Z %in% c(0, 1))) stop("Z must contain only 0, 1")
  
  # Check X and Y: must be 2D structures
  if (ncol(X) < 1) stop("X must have at least one column")
  if (ncol(Y) < 1) stop("Y must have at least one column")
  
  # Check S: if provided, must be a vector or 1D matrix
  if (!is.null(S) && !is.vector(S) && !(is.matrix(S) && ncol(S) == 1)) stop("S must be a vector or a single column matrix")
  
  # Check L: must be a positive integer
  if (!is.numeric(L) || L <= 0 || L != as.integer(L)) stop("L must be a positive integer")
  
  # Check verbose: must be TRUE or FALSE
  if (!is.logical(verbose)) stop("verbose must be TRUE or FALSE")
  
  # Check alpha: must be between 0 and 1
  if (!is.numeric(alpha) || alpha <= 0 || alpha > 1) stop("alpha must be between 0 and 1")
  
  # Check G: must not be NULL
  if (is.null(G)) stop("G cannot be NULL")
  
  # Check covariate_adjustment: must be TRUE or FALSE
  if (!is.logical(covariate_adjustment)) stop("covariate_adjustment must be TRUE or FALSE")
  
  # Check alternative: must be "greater","less" or "two-sided"
  if (!alternative %in% c("greater", "less", "two-sided")) stop("alternative must be 'greater', 'less' or 'two-sided'")
  
  # Check random_state: if provided, must be a positive integer or NULL
  if (!is.null(random_state) && (!is.numeric(random_state) || random_state <= 0 || random_state != as.integer(random_state))) {
    stop("random_state must be a positive integer or NULL")
  }
}

choosemodel <- function(G) {
  #Chooses the imputation model based on the value of G.

  G <- tolower(G) # Convert G to lowercase
  
  if (G == "missforest") {
    # Return missForest imputer function
    return(function(data) missForest::missForest(data)$ximp)
  } else if (G == "mice") {
    # Return MICE imputer function with default method
    return(function(data) complete(mice::mice(data, printFlag = FALSE), 1))
  } else {
    stop("Unsupported imputation method specified")
  }
}

iArt.test <- function(Z, X, Y, G = 'missforest', S = NULL, L = 10000, 
                    verbose = FALSE, covariate_adjustment = FALSE, alpha = 0.05, 
                    alternative = "greater", random_state = NULL) {
  
  # Parameter checks
  check_param(Z, X, Y, S, G, L, verbose, covariate_adjustment, alpha, alternative, random_state)
  if (verbose) print("Parameters checked successfully.")
  
  # Set random seed if provided
  if (!is.null(random_state)) {
    set.seed(random_state)
    if (verbose) print(paste("Random seed set to", random_state))
  }
  
  # Preprocess the variables
  M <- is.na(Y)
  if (is.null(S)) {
    S <- rep(1, nrow(X))
  }
  Z <- matrix(Z, ncol = 1)
  S <- matrix(S, ncol = 1)

  # Choose the imputation model
  G_model <- choosemodel(G)
  if (verbose) print(paste("Imputation model chosen:", G))

  # Impute the missing values to get observed test statistics
  Y_pred <- getY(G_model, Z, X, Y, covariate_adjustment)
  t_obs <- getT(Y_pred, Z, ncol(Y), M)
  if (verbose) print(paste("Observed test statistics =", paste(t_obs, collapse = ", ")))

  # Initialize variable for simulations
  p_values <- numeric(ncol(Y))
  
  # Perform Monte Carlo 
  # Initialize an empty list to store t_sim values
  t_sim_values <- list()

  if (verbose) print("Starting Monte Carlo loop.")
  for (i in 1:L) {
      # Simulate treatment indicators
      Z_sim <- getZsim(getZsimTemplates(Z, S))
      Z_sim <- matrix(Z_sim, ncol = 1)
      
      # Re-impute and calculate test statistics
      Y_pred_sim <- getY(G_model, Z_sim, X, Y, covariate_adjustment)
      t_sim <- getT(Y_pred_sim, Z_sim, ncol(Y), M)

      # Store t_sim values in the list
      t_sim_values[[i]] <- t_sim

      if (verbose) {
          print(paste("Iteration:", i, "- Test statistics =", paste(t_sim, collapse = ", ")))
      }
  }

  # Combine t_sim values into a matrix
  t_sim_matrix <- do.call(rbind, t_sim_values)

  # Calculate mean of t_sim across all iterations
  mean_t_sim <- colMeans(t_sim_matrix)

  # Update p-values
  for (i in 1:nrow(t_sim_matrix)) {
      t_sim <- t_sim_matrix[i,]
      if (alternative == "greater") {
          p_values <- p_values + (t_sim >= t_obs)
      } else if (alternative == "less") {
          p_values <- p_values + (t_sim <= t_obs)
      } else {
          # Here, use mean_t_sim which is the mean across all iterations
          p_values <- p_values + (abs(t_sim - mean_t_sim) >= abs(t_obs - mean_t_sim))
      }
  }

  p_values <- p_values / L
  
  # Holm-Bonferroni correction
  corrected_p_values <- p.adjust(p_values, method = "holm")
  any_rejected <- any(corrected_p_values < alpha)
  
  return(list(reject = any_rejected, p_values = corrected_p_values))
}


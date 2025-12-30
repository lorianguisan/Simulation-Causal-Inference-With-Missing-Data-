#load packages in main first

plot_rejection_rates <- function(
    data_long,
    x_var = "delta",
    y_var = "rej_rate",
    color_var = "curve_type",
    facet_var = NULL,
    fixed_params = list(),
    legend_title = NULL
) {
  
  library(dplyr)
  library(ggplot2)
  library(purrr)
  
  data_plot <- data_long
  
  # to fix fixed params which can be list now
  if (length(fixed_params) > 0) {
    for (param in names(fixed_params)) {
      if (!param %in% names(data_plot)) next
      
      value <- fixed_params[[param]]
      
      if (is.list(data_plot[[param]]) || is.list(value)) {
        # list-column case (e.g., n_blocks)
        data_plot <- data_plot %>%
          filter(map_lgl(.data[[param]], ~ identical(.x, value)))
      } else {
        # normal vector filtering
        data_plot <- data_plot %>% filter(.data[[param]] == value)
      }
    }
  }
  
  if (nrow(data_plot) == 0) stop("Filtering removed all rows — check fixed_params.")
  
  # Force discrete colors if numeric
  if (is.numeric(data_plot[[color_var]])) {
    data_plot[[color_var]] <- factor(data_plot[[color_var]])
  }
  
  # Build subtitle describing fixed params
  subtitle_text <- if (length(fixed_params) > 0) {
    paste0("Fixed: ", paste(paste(names(fixed_params), fixed_params, sep = " = "),
                            collapse = ", "))
  } else NULL
  
  # Split Oracle vs main curves
  oracle_df <- data_plot %>%
    filter(curve_type == "Oracle (true data)") %>%
    distinct(
      !!sym(facet_var) %||% "facet",  # deduplicate per facet
      !!sym(x_var),
      .keep_all = TRUE
    )
  
  main_df <- data_plot %>%
    filter(curve_type != "Oracle (true data)")
  
  # Proper grouping for main curves
  main_group <- interaction(main_df[[color_var]], main_df[[facet_var]] %||% 1)
  
  # Base plot
  p <- ggplot() +
    geom_line(
      data = main_df,
      aes(
        x = !!sym(x_var),
        y = !!sym(y_var),
        color = !!sym(color_var),
        group = main_group
      ),
      linewidth = 1.5   # increased line width
    ) +
    geom_point(
      data = main_df,
      aes(
        x = !!sym(x_var),
        y = !!sym(y_var),
        color = !!sym(color_var),
        group = main_group
      ),
      size = 3          # slightly bigger points
    )
  
  # Add Oracle curve if present
  if (nrow(oracle_df) > 0) {
    oracle_group <- interaction(oracle_df[[facet_var]] %||% 1)
    
    p <- p +
      geom_line(
        data = oracle_df,
        aes(
          x = !!sym(x_var),
          y = !!sym(y_var),
          group = oracle_group
        ),
        color = "black",
        linewidth = 2     # increased line width
      ) +
      geom_point(
        data = oracle_df,
        aes(
          x = !!sym(x_var),
          y = !!sym(y_var),
          group = oracle_group
        ),
        color = "black",
        size = 3          # bigger points
      )
  }
  
  # Formatting
  p <- p +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40") +
    labs(
      x = x_var,
      y = "Rejection Rate",
      color = if (!is.null(legend_title)) legend_title else color_var,
      subtitle = subtitle_text
    ) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
    theme_bw() +
    theme(
      plot.subtitle  = element_text(hjust = 0.5, face = "italic", size = 14),
      plot.title     = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title     = element_text(size = 14),
      axis.text      = element_text(size = 12),
      legend.title   = element_text(size = 13, face = "bold"),
      legend.text    = element_text(size = 12),
      legend.position = "bottom",
      strip.background = element_rect(fill = "gray90")
    )
  
  # Faceting
  if (!is.null(facet_var) && facet_var %in% names(data_plot)) {
    p <- p + facet_wrap(vars(!!sym(facet_var)))
  }
  
  print(p)
} ##got help from chatgpt for ggplot2 plotting

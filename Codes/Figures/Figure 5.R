# Load required libraries
library(ncdf4)
library(ggplot2)
library(dplyr)

# Define the names of the betas as per your model
beta_names <- c("Intercept",
                "TGI",
                "height",
                "real_solar",
                "skyview",
                "shade",
                "shade:TGI",
                "shade:height",
                "shade:real_solar",
                "shade:skyview")

# Function to reshape betas correctly
reshape_betas <- function(betas_array, beta_names) {
  # betas_array dimensions: [betas x samples x chains] = [10 x 2000 x 4]
  
  # Initialize an empty list to store reshaped betas from each chain
  betas_list <- list()
  
  # Loop through each chain
  for (chain in 1:dim(betas_array)[3]) {
    # Extract betas for the current chain: [betas x samples]
    betas_chain <- betas_array[,,chain]
    
    # Transpose to [samples x betas]
    betas_chain_t <- t(betas_chain)
    
    # Convert to data frame and assign column names
    betas_df <- as.data.frame(betas_chain_t)
    colnames(betas_df) <- beta_names
    
    # Append to the list
    betas_list[[chain]] <- betas_df
  }
  
  # Combine all chains into one data frame: [samples * chains] x [betas]
  betas_combined <- bind_rows(betas_list)
  
  return(betas_combined)
}

# Paths to your NetCDF files
path_before_ml <- '/home/ofir/Dropbox/pycharm_projects/pymc_models/Alon/alon_trace_phys.nc'
path_after_ml  <- '/home/ofir/Dropbox/pycharm_projects/pymc_models/Alon/alon_trace_ml.nc'

# Open the NetCDF files
nc_before_ml <- nc_open(path_before_ml)
nc_after_ml  <- nc_open(path_after_ml)

# Extract 'betas' posterior samples
# Replace "posterior/betas" with the actual variable name if different
betas_before_ml_samples <- ncvar_get(nc_before_ml, "posterior/betas")  # [10 x 2000 x 4]
betas_after_ml_samples  <- ncvar_get(nc_after_ml, "posterior/betas")   # [10 x 2000 x 4]

# Close the NetCDF files
nc_close(nc_before_ml)
nc_close(nc_after_ml)

# Reshape the betas into data frames
betas_before_ml <- reshape_betas(betas_before_ml_samples, beta_names)  # [8000 x 10]
betas_after_ml  <- reshape_betas(betas_after_ml_samples, beta_names)   # [8000 x 10]

# Inspect the first few rows to verify
head(betas_before_ml)
head(betas_after_ml)


#create dataset for predictions
data  = read.csv("regression_table.csv")
data$shade = as.factor(data$shade)
# please create a table with all the options for the variables using expand.grid: TGI, height, shade, real_solar, skyview

# Define the number of points for each variable in the prediction grid
n_points <- 50  # Adjust based on your needs and computational resources

generate_predictions_for_variable <- function(betas_before, betas_after, data, variable, n_points = 50) {
  # Ensure 'shade' is a factor
  data$shade <- as.factor(data$shade)
  
  # Identify other continuous variables
  vars_names = c("TGI", "real_solar", "height", "skyview", "shade")
  continuous_vars <- setdiff(c(vars_names), c(variable, "shade"))
  
  # Calculate median values for variables to fix
  fixed_values <- data %>%
    summarize(across(all_of(continuous_vars), mean, na.rm = TRUE))
  
  # Convert fixed_values to a named list for easier access
  fixed_values_list <- as.list(fixed_values)
  # Create a sequence for the variable to vary
  var_seq <- seq(
    from = min(data[[variable]], na.rm = TRUE),
    to = max(data[[variable]], na.rm = TRUE),
    length.out = n_points
  )
  
  # Get the levels of 'shade'
  shade_levels <- levels(data$shade)
  
  # Create the prediction grid using expand.grid
  prediction_grid <- expand.grid(
    shade = shade_levels,
    variable = var_seq
  )
  
  # Rename the 'variable' column to the actual variable name
  names(prediction_grid)[names(prediction_grid) == "variable"] <- variable
  
  #browser()
  
  # Add the fixed variables by directly mutating them using the named list
  # This avoids trying to mutate non-existent columns
  prediction_grid <- prediction_grid %>%
    mutate(!!!fixed_values_list)  
  # Reorder columns for consistency
  prediction_grid <- prediction_grid %>%
    select(vars_names)
  
  # Initialize a list to store summaries
  summary_list_before <- list()
  summary_list_after  <- list()
  
  # Loop through each row of the prediction grid to generate predictions
  for (i in 1:nrow(prediction_grid)) {
    # Extract the current scenario
    scenario <- prediction_grid[i, ]
    
    # Handle interaction terms based on 'shade' level
    # Assuming 'Low' is the reference level
    is_shade_high <- ifelse(scenario$shade == shade_levels[1], FALSE, TRUE)
    
    # Create interaction terms
    interaction_terms <- data.frame(
      `shade:TGI` = ifelse(is_shade_high, scenario$TGI, 0),
      `shade:height` = ifelse(is_shade_high, scenario$height, 0),
      `shade:real_solar` = ifelse(is_shade_high, scenario$real_solar, 0),
      `shade:skyview` = ifelse(is_shade_high, scenario$skyview, 0)
    )
    
    # Combine scenario with interaction terms
    prediction_row <- cbind(
      scenario %>% select(-shade),
      interaction_terms
    )
    
    # Ensure all required betas are present
    required_betas <- c("Intercept", variable, continuous_vars, "shade",
                        "shade:TGI", "shade:height", "shade:real_solar", "shade:skyview")
    
    # Calculate the mean betas for shade and the main variable
    mean_beta_variable_before <- mean(betas_before[[variable]])
    mean_beta_shade_before <- mean(betas_before$shade)
    
    mean_beta_variable_after <- mean(betas_after[[variable]])
    mean_beta_shade_after  <- mean(betas_after$shade)
    
    # **Compute predictions before ML correction**
    # If shade is High, use mean beta for both variable and shade; otherwise, use normal betas
    pred_before <- if (is_shade_high) {
      mean_beta_variable_before * prediction_row[1,variable] + 
        mean_beta_shade_before +
        betas_before[,paste0("shade:", variable)] * prediction_row[1, paste0("shade.", variable)]
    } else {
      betas_before[,variable] * prediction_row[1,variable]
    }
    
    # **Compute predictions after ML correction**
    pred_after <- if (is_shade_high) {
      mean_beta_variable_after * prediction_row[1,variable]+ 
        mean_beta_shade_after + 
        betas_after[,paste0("shade:", variable)] * prediction_row[1,variable]
    } else {
      betas_after[, variable] * prediction_row[1,variable]
    }
    # Summarize predictions before ML correction
    #browser()
        
    summary_before <- data.frame(
      Variable = variable,
      shade = scenario$shade,
      value = scenario[[variable]],
      mean = mean(pred_before),
      lower = quantile(pred_before, 0.025),
      upper = quantile(pred_before, 0.975),
      Correction = "Before ML"
    )
    
    # Summarize predictions after ML correction
    summary_after <- data.frame(
      Variable = variable,
      shade = scenario$shade,
      value = scenario[[variable]],
      mean = mean(pred_after),
      lower = quantile(pred_after, 0.025),
      upper = quantile(pred_after, 0.975),
      Correction = "After ML"
    )
    
    # Append to the lists
    summary_list_before[[i]] <- summary_before
    summary_list_after[[i]]  <- summary_after
  }
  
  # Combine all summaries into data frames
  summarized_before_ml <- bind_rows(summary_list_before)
  summarized_after_ml  <- bind_rows(summary_list_after)
  
  #Remove rows where the variable is not within the observed range for both shade and open**
  # And combine before and after into one data frame
  observed_range_shade <- range(data[data$shade == 1, variable], na.rm = TRUE)
  observed_range_open <- range(data[data$shade == 0, variable], na.rm = TRUE)
  
  summarized_predictions <- bind_rows(summarized_before_ml, summarized_after_ml) %>%
    filter((shade == 1 & value >= observed_range_shade[1] & value <= observed_range_shade[2]) |
             (shade == 0 & value >= observed_range_open[1] & value <= observed_range_open[2]))
  
  # Return the summarized predictions
  return(summarized_predictions)
}



# Define the continuous variables to vary
variables_to_vary <- c("TGI", "height", "real_solar", "skyview")

# Initialize a list to store all summarized predictions
all_summarized_predictions <- list()

# Loop through each variable and generate predictions
for (var in variables_to_vary) {
  cat("Generating predictions for:", var, "\n")
  
  # Generate summarized predictions
  summarized_pred <- generate_predictions_for_variable(
    betas_before = betas_before_ml,
    betas_after = betas_after_ml,
    data = data,
    variable = var,
    n_points = 50
  )
  
  # Append to the list
  all_summarized_predictions[[var]] <- summarized_pred
}

# Combine all summarized predictions into one data frame
all_summarized_predictions_df <- bind_rows(all_summarized_predictions)

# Inspect the summarized predictions
head(all_summarized_predictions_df)


# Define the plotting function to use with facets
plot_predictions_facet <- function(summarized_data) {
  # Custom color palette for shade levels
  custom_colors <- c("shade" = "black", "open" = "orange")
  
  # Ensure the Correction factor and Shade levels are correctly formatted
  summarized_data <- summarized_data %>%
    mutate(shade_label = ifelse(shade == 0, "open", "shade"))
  
  # Create custom labels for facets
  facet_labels <- c("(A)", "(B)", "(C)", "(D)")
  
  
  p = ggplot() +
    # Dashed line for "Before ML" and shade == "open"
    geom_line(data = summarized_data %>% filter(Correction == "Before ML", shade_label == "open"),
              aes(x = value, y = mean, linetype = "Before ML"), color = "orange") +
    # Solid line for "After ML" and shade == "open"
    geom_line(data = summarized_data %>% filter(Correction == "After ML", shade_label == "open"),
              aes(x = value, y = mean, linetype = "After ML"), color = "orange") +
    
    # Dashed line for "Before ML" and shade == "shade"
    geom_line(data = summarized_data %>% filter(Correction == "Before ML", shade_label == "shade"),
              aes(x = value, y = mean, linetype = "Before ML"), color = "black") +
    # Solid line for "After ML" and shade == "shade"
    geom_line(data = summarized_data %>% filter(Correction == "After ML", shade_label == "shade"),
              aes(x = value, y = mean, linetype = "After ML"), color = "black") +
    
    # Add ribbon for confidence intervals - Before ML
    geom_ribbon(data = summarized_data %>% filter(Correction == "Before ML"),
                aes(x = value, ymin = lower, ymax = upper, fill = shade_label),
                alpha = 0.2, linetype = "dashed", color = NA) +
    # Add ribbon for confidence intervals - After ML
    geom_ribbon(data = summarized_data %>% filter(Correction == "After ML"),
                aes(x = value, ymin = lower, ymax = upper, fill = shade_label),
                alpha = 0.2, linetype = "solid", color = NA) +
    
    # Custom colors
    scale_fill_manual(values = custom_colors) +
    
    # Ensure linetype is included in the legend
    scale_linetype_manual(values = c("Before ML" = "dashed", "After ML" = "solid")) +
    
    # Labels and theme
    labs(
      title = "Predictions by Variable",
      x = "Value",
      y = "Predicted Response",
      fill = "Shade",  # Legend for fill
      linetype = "Correction"  # Legend for linetype
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add box around the plot
      axis.line = element_line(color = "black"),  # Add axis lines
      axis.ticks = element_line(color = "black"),  # Add axis ticks
      axis.text.x = element_text(color = "black"),  # Ensure x-axis text is displayed
      axis.text.y = element_text(color = "black")   # Ensure y-axis text is displayed
    ) +
    # Facet the plot by Variable, arranged in a 2x2 grid, with free x-axis scales
    facet_wrap(~Variable, scales = c("free"), ncol = 2)
    
    # Add annotations for the panel labels inside the panels at top left
      p + annotate("text", x = -Inf, y = Inf, label = facet_labels[1], hjust = -0.1, vjust = 1.5, size = 5, data = summarized_data %>% filter(Variable == unique(summarized_data$Variable)[1])) +
      annotate("text", x = -Inf, y = Inf, label = facet_labels[2], hjust = -0.1, vjust = 1.5, size = 5, data = summarized_data %>% filter(Variable == unique(summarized_data$Variable)[2])) +
      annotate("text", x = -Inf, y = Inf, label = facet_labels[3], hjust = -0.1, vjust = 1.5, size = 5, data = summarized_data %>% filter(Variable == unique(summarized_data$Variable)[3])) +
      annotate("text", x = -Inf, y = Inf, label = facet_labels[4], hjust = -0.1, vjust = 1.5, size = 5, data = summarized_data %>% filter(Variable == unique(summarized_data$Variable)[4]))
    
  
}

# Call the plotting function with the full summarized dataset
p <- plot_predictions_facet(all_summarized_predictions_df)
print(p)


# Define the plotting function for a single variable
plot_predictions_single <- function(summarized_data, variable, variable_name) {
  custom_colors <- c("shade" = "black", "open" = "orange")
  
  # Ensure the Correction factor and Shade levels are correctly formatted
  summarized_data <- summarized_data %>%
    mutate(shade_label = ifelse(shade == 0, "open", "shade"))
  
  # Filter data for the specific variable
  data_to_plot <- summarized_data %>% filter(Variable == variable)
  
  # Create plot
  p <- ggplot(data_to_plot) +
    geom_line(data = data_to_plot %>% filter(Correction == "Before ML", shade_label == "open"),
              aes(x = value, y = mean, linetype = "Before ML"), color = "orange") +
    geom_line(data = data_to_plot %>% filter(Correction == "After ML", shade_label == "open"),
              aes(x = value, y = mean, linetype = "After ML"), color = "orange") +
    geom_line(data = data_to_plot %>% filter(Correction == "Before ML", shade_label == "shade"),
              aes(x = value, y = mean, linetype = "Before ML"), color = "black") +
    geom_line(data = data_to_plot %>% filter(Correction == "After ML", shade_label == "shade"),
              aes(x = value, y = mean, linetype = "After ML"), color = "black") +
    
    geom_ribbon(data = data_to_plot %>% filter(Correction == "Before ML"),
                aes(x = value, ymin = lower, ymax = upper, fill = shade_label), 
                alpha = 0.2, linetype = "dashed", color = NA) +
    geom_ribbon(data = data_to_plot %>% filter(Correction == "After ML"),
                aes(x = value, ymin = lower, ymax = upper, fill = shade_label), 
                alpha = 0.2, linetype = "solid", color = NA) +
    
    scale_fill_manual(values = custom_colors) +
    scale_linetype_manual(values = c("Before ML" = "dashed", "After ML" = "solid")) +
    
    labs(
      title = "",
      x = variable_name, 
      y = "Microclimate Bias (°C)",
      fill = "",
      linetype = ""
    ) +
    theme_minimal() +
    theme( axis.text = element_text(size = 12), legend.text = element_text(size = 14), axis.title = element_text(size = 14),
      legend.position = "bottom",
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add box around the plot
      axis.line = element_line(color = "black"),  # Add axis lines
      axis.ticks = element_line(color = "black"),  # Add axis ticks
      axis.text.x = element_text(color = "black"),  # Ensure x-axis text is displayed
      axis.text.y = element_text(color = "black")   # Ensure y-axis text is displayed
    ) +
    theme(legend.position = "bottom")
  
  return(p)
}

# Create separate plots for each variable
p1 <- plot_predictions_single(summarized_data = all_summarized_predictions_df, variable = "TGI", variable_name = "TGI")
p2 <- plot_predictions_single(summarized_data = all_summarized_predictions_df, variable = "real_solar", variable_name = "Solar Radiation (W/m²)")
p3 <- plot_predictions_single(summarized_data = all_summarized_predictions_df, variable = "height", variable_name = "Height (m)")
p4 <- plot_predictions_single(summarized_data = all_summarized_predictions_df, variable = "skyview", variable_name = "Skyview (%)")

library(ggpubr)

jpeg(filename = "Figure 5.jpeg", width=2800, height = 2500, res=300)
# Arrange the plots in a 2x2 grid using ggarrange
ggarrange(p1, p2, p3, p4,
          labels = c("     (a)", "     (b)", "     (c)", "     (d)"),   # Add panel labels
          ncol = 2, nrow = 2,                      # 2x2 grid
          common.legend = TRUE,                    # Use a common legend
          legend = "bottom")                       # Place legend at the bottom
dev.off()

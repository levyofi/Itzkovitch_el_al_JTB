# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
# Example data
data <- read.csv("../Tables/RF_CorrectionModel_summary.csv", header = TRUE)

# Create separate dataframes for ME, MAE, and MSE
df_ME <- data %>%
  select(Map, M1_ME, M2_ME, M3_ME)

df_MAE <- data %>%
  select(Map, M1_MAE, M2_MAE, M3_MAE)

df_MSE <- data %>%
  select(Map, M1_MSE, M2_MSE, M3_MSE)

# Add a common function to generate the plot for a metric
# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)

# Custom function to reshape data and calculate improvement status
plot_metric <- function(plotdata, y_label) {
  
  # Get the column names for stages M1, M2, and M3
  col_names <- names(plotdata)
  
  # Calculate if the error improved from M1 to M2, and from M1 to M3
  plotdata$m2improved <- abs(plotdata[, col_names[2]]) > abs(plotdata[, col_names[3]])  # M1 > M2
  plotdata$m3improved <- abs(plotdata[, col_names[2]]) > abs(plotdata[, col_names[4]])  # M1 > M3
  
  # Reshape data from wide to long format
  long_data <- plotdata %>%
    pivot_longer(cols = col_names[2:4], names_to = "Stage", values_to = "Value")
  
  write.csv(long_data, file = paste0(col_names[2],".csv"), row.names = FALSE)
  #browser()
  # Combine the improvement status into the long data
  long_data$Improvement <- ifelse( (long_data$Stage == col_names[3] | long_data$Stage==col_names[2]), long_data$m2improved, long_data$m3improved)
  
  # Generate jitter for each point in the dataset to match dots and lines
  long_data <- long_data %>%
    mutate(jittered_x = as.numeric(factor(Stage)) + runif(n(), min = -0.05, max = 0.05))
  
  # Calculate means for each stage
  means <- long_data %>%
    group_by(Stage) %>%
    summarise(mean_value = mean(Value))
  
  ylims = range(long_data$Value)
  # Create the plot
  #long_data[long_data$Stage=="M1_ME" | long_data$Stage=="M2_ME",]
  ggplot() +
    geom_boxplot(data = long_data, aes(group = Stage, x = as.numeric(factor(Stage)), y = Value), outlier.shape = NA, alpha = 0.7, width =0.5) +
    geom_line(data = long_data[long_data$Stage==col_names[2] | long_data$Stage==col_names[3],], aes(x = jittered_x, y = Value, color = Improvement, group = Map), alpha = 0.5) +
    geom_line(data = long_data[long_data$Stage==col_names[3] | long_data$Stage==col_names[4],], aes(x = jittered_x, y = Value, color = Improvement, group = Map), alpha = 0.5) +
    geom_point(data = long_data, aes(x = jittered_x, y = Value), color = "black", size = 2) +
    scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "red"), labels = c("Yes", "No"),  name = "Improvement") +  # Color: green for improved, red for worsened
    labs(y = y_label, x = NULL) + 
    scale_x_continuous(breaks = c(1, 2, 3), labels = c("Physical model", "After ML", "After ML + PE")) +
    theme_minimal() + ylim(ylims[1], ylims[2]*1.2) + 
    theme( axis.text = element_text(size = 12), legend.text = element_text(size = 14), axis.title = element_text(size = 14),
           legend.position = "none",
           panel.grid.major = element_blank(),  # Remove major grid lines
           panel.grid.minor = element_blank(),  # Remove minor grid lines
           panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add box around the plot
           axis.line = element_line(color = "black"),  # Add axis lines
           axis.ticks = element_line(color = "black"),  # Add axis ticks
           axis.text.x = element_text(color = "black"),  # Ensure x-axis text is displayed
           axis.text.y = element_text(color = "black")   # Ensure y-axis text is displayed
    )
}

# Apply the function for each metric (for example, ME)
p_ME <- plot_metric(df_ME, "Mean Error (°C)")
p_MAE <- plot_metric(df_MAE, "Mean Absolute Error (°C)")
p_MSE <- plot_metric(df_MSE, "Mean Square Error (°C)")

# To combine the plots
library(ggpubr)

jpeg(filename = "Figure 3.jpg", width=1200, height = 3000, res=300)
# Arrange the plots in a 2x2 grid using ggarrange
ggarrange(p_ME, p_MAE, p_MSE,
          labels = c("(a)", "(b)", "(c)"),  label.y = 0.97 ,label.x = 0.13, # Add panel labels
          nrow = 3,                      # 2x2 grid
          common.legend = TRUE,                    # Use a common legend
          legend = "bottom", align = "v")   +                    # Place legend at the bottom
  theme(plot.margin = margin(0.1,0.3,1,0.1, "cm"))
dev.off()

# run statistical tests # we ended up using Bayesian modeling
library(nlme)
ME = read.csv("M1_ME.csv")
hist(ME$Value)
model = lme(Value~Stage, data = ME, random=~1|Map)
summary(model)
model = lme(Value~Stage, data = ME, random=~1|Map, weights = varIdent(form=~1|Stage))
summary(model)

MAE = read.csv("M1_MAE.csv")
hist(MAE$Value)
library(lme4)
model = glmer(Value~Stage + (1|Map), data = MAE, family = Gamma(link="log"))
summary(model)

MSE = read.csv("M1_MSE.csv")
hist(MSE$Value)
library(lme4)
model = glmer(Value~Stage + (1|Map), data = MSE, family = Gamma(link="log"))
summary(model)

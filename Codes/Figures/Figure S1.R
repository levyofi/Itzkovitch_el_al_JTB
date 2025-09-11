library(ggplot2)
library(reshape2)
library(ggpubr)

# Load the data
df <- read.csv('../Example data/Tables/online_met_dbs_comparison.csv')

# Reshape the data for ggplot2
df_long <- melt(df, measure.vars = c("Station_ME", "ERA5_ME", "Station_MAE", "ERA5_MAE"))

# Rename the columns for easier plotting
df_long$Model <- ifelse(grepl("ERA5", df_long$variable), "ERA5", "GLDAS")
df_long$Metric <- ifelse(grepl("ME", df_long$variable), "ME", "MAE")


# Base plot for both ME and MAE with no legend and no X-axis title
p_me <- ggplot(df_long[df_long$Metric=="ME",], aes(x = Model, y = value, color = Model)) +
  geom_boxplot(outlier.shape = NA) +   # Add boxplot
  geom_jitter(width = 0.2, alpha = 0.3) +  # Add jitter
  scale_color_manual(values = c("GLDAS" = "blue", "ERA5" = "orange")) +  # Set custom colors
  geom_hline(yintercept = 0, linetype = "dashed") +  # Add horizontal line at 0
  theme_minimal() + labs(y = "Mean Error (°C)")  +
  theme(panel.border = element_rect(colour = "black", fill = NA),  # Add box around each plot
        axis.text = element_text(size = 12),  # Adjust axis text size
        axis.ticks = element_line(size = 1),  # Add thicker axis ticks
        panel.grid.major = element_blank(),   # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines,
        legend.position = "none",  # Remove legends from individual plots
        axis.title.x = element_blank(),  # Remove X-axis title
        plot.title = element_blank())  # Remove panel titles

p_mae <- ggplot(df_long[df_long$Metric=="MAE",], aes(x = Model, y = value, color = Model)) +
  geom_boxplot(outlier.shape = NA) +   # Add boxplot
  geom_jitter(width = 0.2, alpha = 0.3) +  # Add jitter
  scale_color_manual(values = c("GLDAS" = "blue", "ERA5" = "orange")) +  # Set custom colors
  geom_hline(yintercept = 0, linetype = "dashed") +  # Add horizontal line at 0
  theme_minimal() + labs(y = "Mean Absolute Error (°C)") +
  theme(panel.border = element_rect(colour = "black", fill = NA),  # Add box around each plot
        axis.text = element_text(size = 12),  # Adjust axis text size
        axis.ticks = element_line(size = 1),  # Add thicker axis ticks
        panel.grid.major = element_blank(),   # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines,
        legend.position = "none",  # Remove legends from individual plots
        axis.title.x = element_blank(),  # Remove X-axis title
        plot.title = element_blank())  # Remove panel titles


# Now add a shared legend below using ggarrange
legend <- get_legend(ggplot(df_long, aes(x = Model, y = value, color = Model)) +
                       geom_boxplot() +
                       scale_color_manual(values = c("GLDAS" = "blue", "ERA5" = "orange")) +
                       theme(legend.position = "bottom"))

# Arrange the plots with the shared legend and no panel titles or x-axis title
final_plot <- ggarrange(p_me, p_mae, ncol = 2, common.legend = TRUE, legend = "bottom")

# Display the final plot
jpeg(filename = "Figure S3.jpg", width=3000, height = 1400, res=300)
print(final_plot)
dev.off()


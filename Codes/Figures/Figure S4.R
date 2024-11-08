library(ggplot2)
library(reshape2)
library(ggpubr)

# Load the data from an Excel file
df <- readxl::read_excel('../Example data/Tables/ML_NN_comparison.xlsx')
df = df[,c("NNM1_ME", "Phy_ME", "RF_ME", "NNM2_ME", "NNM1_MAE", "Phy_MAE", "RF_MAE", "NNM2_MAE")]
# Parameters for plotting
params <- c('ME', 'MAE')

# Reshape the data for ggplot2 (long format), selecting only the specified columns
df_long <- reshape2::melt(df, measure.vars = c("Phy_ME", "RF_ME", "NNM2_ME", "Phy_MAE", "RF_MAE", "NNM2_MAE"))
# Rename columns for easier plotting
df_long$Model <- factor(ifelse(grepl("Phy", df_long$variable), "Physical",
                               ifelse(grepl("RF", df_long$variable), "Random Factor", "DNN")),
                        levels = c("Physical", "Random Factor", "DNN"))
df_long$Metric <- ifelse(grepl("ME", df_long$variable), "ME", "MAE")

# Base plot for both ME and MAE with no legend and no X-axis title, jitter with transparency
p_me <- ggplot(df_long[df_long$Metric=="ME",], aes(x = Model, y = value, color = Model)) +
  geom_boxplot(outlier.shape = NA) +   # Add boxplot
  geom_jitter(width = 0.2, alpha = 0.3) +  # Add jitter with transparency (alpha)
  scale_color_manual(values = c("Physical" = "orange", "Random Factor" = "#006400", "DNN" = "red")) +  # Set custom colors
  geom_hline(yintercept = 0, linetype = "dashed") +  # Add horizontal line at 0
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA),  # Add box around each plot
        axis.text = element_text(size = 12),  # Adjust axis text size
        axis.ticks = element_line(size = 1),  # Add thicker axis ticks
        panel.grid.major = element_blank(),   # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines,
        legend.position = "none",  # Remove legends from individual plots
        axis.title.x = element_blank(),  # Remove X-axis title
        plot.title = element_blank()) + # Remove panel titles
  labs(y = "ME (℃)") + 
  coord_cartesian(ylim = c(min(df_long$value[df_long$Metric == "ME"]), 
                           max(df_long$value[df_long$Metric == "ME"])))

p_mae <- ggplot(df_long[df_long$Metric=="MAE",], aes(x = Model, y = value, color = Model)) +
  geom_boxplot(outlier.shape = NA) +   # Add boxplot
  geom_jitter(width = 0.2, alpha = 0.3) +  # Add jitter with transparency (alpha)
  scale_color_manual(values = c("Physical" = "orange", "Random Factor" = "#006400", "DNN" = "red")) +  # Set custom colors
  geom_hline(yintercept = 0, linetype = "dashed") +  # Add horizontal line at 0
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA),  # Add box around each plot
        axis.text = element_text(size = 12),  # Adjust axis text size
        axis.ticks = element_line(size = 1),  # Add thicker axis ticks
        panel.grid.major = element_blank(),   # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines,
        legend.position = "none",  # Remove legends from individual plots
        axis.title.x = element_blank(),  # Remove X-axis title
        plot.title = element_blank()) + # Remove panel titles
  labs(y = "MAE (℃)") + 
  coord_cartesian(ylim = c(min(df_long$value[df_long$Metric == "MAE"]), 
                           max(df_long$value[df_long$Metric == "MAE"])))

# Now add a shared legend below using ggarrange
legend <- get_legend(ggplot(df_long, aes(x = Model, y = value, color = Model)) +
                       geom_boxplot() +
                       scale_color_manual(values = c("Physical" = "orange", "Random Factor" = "green", "DNN" = "red")) +
                       theme(legend.position = "bottom"))

# Arrange the plots with the shared legend and no panel titles or x-axis title
jpeg(filename = "Figure S4.jpg", width=3000, height = 1400, res=300)
ggarrange(p_me, p_mae, ncol = 2, common.legend = TRUE, legend = "bottom",labels = c("(a)", "(b)"),
          label.y = 0.95 ,label.x = 0.08, # Add panel labels
          align = "hv")   +                    # Place legend at the bottom
  theme(plot.margin = margin(0.1,0.5,1,0.1, "cm"))
dev.off()

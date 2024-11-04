library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)

data = fread("/home/ofir/Downloads/SamplePixels_all.csv", select = c("beforeML", "afterML", "afterML+mPE"))


# Transform data from wide to long format using melt from data.table
long_data <- melt(data, id.vars = NULL, variable.name = "Stage", value.name = "Error")

# Bin the data by 0.25 intervals
long_data[, Bin := floor(Error / 0.25) * 0.25]
long_data[, AbsBin := floor(abs(Error) / 0.25) * 0.25]


# Summarize data: calculate counts and totals for percentages
summary_data <- long_data[, .(Count = .N), by = .(Stage, Bin)]
summary_data[, Total := sum(Count), by = Stage]
summary_data[, Percentage := (Count / Total) * 100]
summary_data[, CumulativePercentage := cumsum(Percentage), by = Stage]

# Summarize data for absolute errors: calculate counts and totals for percentages
summary_data_abs <- long_data[, .(Count = .N), by = .(Stage, AbsBin)]
summary_data_abs = summary_data_abs[order(Stage, AbsBin)]
summary_data_abs[, Total := sum(Count), by = Stage]
summary_data_abs[, Percentage := (Count / Total) * 100]
summary_data_abs[, CumulativePercentage := cumsum(Percentage), by = Stage]

a_panel = ggplot() +
  geom_area(data = subset(summary_data, Stage == "beforeML"),
            aes(x = Bin, y = Percentage, fill = Stage), alpha = 0.3) +
  geom_line(data = subset(summary_data, Stage != "beforeML"),
            aes(x = Bin, y = Percentage, color = Stage, linetype = Stage), size=1) +  # Include linetype in aes
  scale_color_manual(name="Stage", values = c("beforeML" = NA, "afterML" = "black", "afterML+mPE" = "black")) +  # Custom colors for line plots
  scale_linetype_manual(name="Stage", values = c("beforeML" = NA, "afterML" = "solid", "afterML+mPE" = "dotted")) +  # Correctly map linetypes
  scale_fill_manual(name="Stage", values = c("beforeML" = "grey30", "afterML" = NA, "afterML+mPE" = NA)) +  # Custom color for area plot
  scale_x_continuous(breaks = seq(-18, 18, by = 3)) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), legend.text = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position = "right",  # Adjust legend position to show it
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add box around the plot
        axis.line = element_line(color = "black"),  # Add axis lines
        axis.ticks = element_line(color = "black"),  # Add axis ticks
        axis.text.x = element_text(color = "black"),  # Ensure x-axis text is displayed
        axis.text.y = element_text(color = "black")   # Ensure y-axis text is displayed
  ) +
  labs(title = "",
       x = "Error (°C)",
       y = "Percentage (%)",
       color = "Stage", linetype = "Stage")  # Ensure legend labels match aes mappings


b_panel = ggplot() +
  geom_area(data = subset(summary_data_abs, Stage == "beforeML"),
            aes(x = AbsBin, y = CumulativePercentage, fill = Stage), alpha = 0.3) +
  geom_line(data = subset(summary_data_abs, Stage != "beforeML"),
            aes(x = AbsBin, y = CumulativePercentage, color = Stage, linetype = Stage), size=1) +  # Include linetype in aes
  scale_fill_manual(name="", values = c("beforeML" = "grey30", "afterML" = NA, "afterML+mPE" = NA), labels = c("Physical model", "", "")) +  # Custom color for area plot
  scale_color_manual(name="Stage", values = c("beforeML" = NA, "afterML" = "black", "afterML+mPE" = "black")) +  # Custom colors for line plots
  scale_linetype_manual(name="Stage", values = c("beforeML" = NA, "afterML" = "solid", "afterML+mPE" = "dotted")) +  # Correctly map linetypes
  scale_x_continuous(breaks = seq(-18, 18, by = 2))  +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), legend.text = element_text(size = 14), axis.title = element_text(size = 14),
        legend.position = "right",  # Adjust legend position to show it
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add box around the plot
        axis.line = element_line(color = "black"),  # Add axis lines
        axis.ticks = element_line(color = "black"),  # Add axis ticks
        axis.text.x = element_text(color = "black"),  # Ensure x-axis text is displayed
        axis.text.y = element_text(color = "black"),   # Ensure y-axis text is displayed
  ) +
  labs(title = "",
       x = "Absoute Error (°C)",
       y = "Percentage (%)",
       color = "Stage", linetype = "Stage")  # Ensure legend labels match aes mappings

library(ggpubr)

jpeg(filename = "Figure 4.jpg", width=3000, height = 1400, res=300)
ggarrange(a_panel, b_panel, ncol = 2, legend = FALSE, labels = c("(a)", "(b)"),
          label.y = 0.9 ,label.x = 0.12, # Add panel labels
          align = "hv")   +                    # Place legend at the bottom
  theme(plot.margin = margin(0.1,0.3,1,0.1, "cm"))
dev.off()

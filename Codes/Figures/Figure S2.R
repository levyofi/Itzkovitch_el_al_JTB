library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(ggpubr)

# "Zeelim_18.9.19_1300", "Zeelim_30.1.20_0810"
"/home/ofir/Dropbox/eclipse workspace/lab/Alon/cropped_maps//Zeelim_23.9.19_0950" 
#[14] "/home/ofir/Dropbox/eclipse workspace/lab/Alon/cropped_maps//Zeelim_23.9.19_1100" 
# [18] "/home/ofir/Dropbox/eclipse workspace/lab/Alon/cropped_maps//Zeelim_23.9.19_1700" 
# [19] "/home/ofir/Dropbox/eclipse workspace/lab/Alon/cropped_maps//Zeelim_23.9.19_1730" 
# [34] "/home/ofir/Dropbox/eclipse workspace/lab/Alon/cropped_maps//Zeelim_31.05.21_0611"
# [35] "/home/ofir/Dropbox/eclipse workspace/lab/Alon/cropped_maps//Zeelim_31.05.21_0707"
# [36] "/home/ofir/Dropbox/eclipse workspace/lab/Alon/cropped_maps//Zeelim_31.05.21_0819"
# [37] "/home/ofir/Dropbox/eclipse workspace/lab/Alon/cropped_maps//Zeelim_31.05.21_0924"
# [38] "/home/ofir/Dropbox/eclipse workspace/lab/Alon/cropped_maps//Zeelim_31.05.21_1033"
# [39] "/home/ofir/Dropbox/eclipse workspace/lab/Alon/cropped_maps//Zeelim_31.05.21_1140"

# Define the train maps
train_maps <- c("Zeelim_29.5.19_0830", "Zeelim_29.5.19_1650", "Zeelim_29.5.19_1730", 
                "Zeelim_30.5.19_0600", "Zeelim_30.5.19_0630", "Zeelim_18.9.19_0900", 
                "Zeelim_18.9.19_1200", "Zeelim_18.9.19_1400", 
                "Zeelim_18.9.19_1500", "Zeelim_18.9.19_1720", "Zeelim_7.11.19_1030", 
                "Zeelim_7.11.19_1100", "Zeelim_7.11.19_1310", "Zeelim_7.11.19_1550", 
                "Zeelim_7.11.19_1640", "Zeelim_30.1.20_0920", 
                "Zeelim_30.1.20_0950", "Zeelim_30.1.20_1050", "Zeelim_30.1.20_1200", 
                "Zeelim_30.1.20_1300", "Zeelim_30.1.20_1350", "Zeelim_30.1.20_1449", 
                "Zeelim_30.1.20_1523")
test_maps <- c("Zeelim_12.04.21_1118", "Zeelim_23.9.19_0610", "Zeelim_23.9.19_0700", "Zeelim_23.9.19_0800", "Zeelim_23.9.19_1410", 
"Zeelim_23.9.19_1510", "Zeelim_23.9.19_1610", "Zeelim_31.05.21_1516", "Zeelim_31.05.21_1712",
"Zeelim_31.05.21_1805")

# Parameter name mapping for more descriptive labels
param_labels <- list(
  TG = "Ground Temperature (°K)",
  TAIR = "Air Temperature (°K)",
  Albedo = "Albedo (%)",
  Cloud_cover = "Cloud Cover (%)",
  Wind = "Wind Speed (m/s)",
  P = "Pressure (hPa)",
  QAIR = "Specific Humidity (kg/kg)"
)

# Initialize the data frame
dct <- data.frame(DataSet = character(), Value = numeric(), Type = character(), stringsAsFactors = FALSE)

# Process the CSV files
files <- list.files("/home/ofir/Dropbox/eclipse workspace/lab/Alon/met_data/", pattern = "Zeelim.*\\_new.csv", full.names = TRUE)

for (file in files) {
  name <- strsplit(basename(file), "_input")[[1]][1]
  if (name %in% c(train_maps, test_maps)){
    ds <- ifelse(name %in% train_maps, "Train", "Test")
    mdf <- read_csv(file)
    
    # Add the parameters to the dataframe
    for (param in c('TG', 'TAIR', 'Albedo', 'Cloud_cover', 'Wind', 'P', "QAIR")) {
      p1 <- if (param != 'P') mdf[[param]][1] else mdf[[param]][1] / 100000
      dct <- rbind(dct, data.frame(DataSet = ds, Value = p1, Type = param, stringsAsFactors = FALSE))
    }
  }
}

# Plot Figure S2 with informative x-axis titles
plot_list <- list()
for (type in unique(dct$Type)) {
  p <- ggplot(dct %>% filter(Type == type), aes(x = Value)) +
    geom_histogram(aes(y = after_stat(100 * count / sum(count))), bins = 20, fill = "blue", alpha = 0.7) +
    facet_wrap(~Type, scales = "free") +
    xlab(param_labels[[type]]) +         # Use descriptive x-axis label
    ylab("Percentage (%)") +             # Set y-axis label to "Percentage (%)"
    theme_minimal() +
    theme(strip.text = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),  # Add box around each plot
          axis.text = element_text(size = 12),  # Adjust axis text size
          axis.ticks = element_line(size = 1),  # Add thicker axis ticks
          panel.grid.major = element_blank(),   # Remove major gridlines
          panel.grid.minor = element_blank())  # Remove the facet strip text (top axis text)
  
  plot_list[[type]] <- p
}

# Arrange the plots in grid using ggarrange
final_plot = ggarrange(plotlist = plot_list, ncol = 3, nrow = ceiling(length(plot_list) / 3))

# Display the final plot
jpeg(filename = "Figure S2.jpg", width=2700, height = 2500, res=300)
print(final_plot)
dev.off()

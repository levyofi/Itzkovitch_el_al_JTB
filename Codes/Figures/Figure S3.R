library(raster)
library(dplyr)
library(tidyr)

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

# Initialize a list to store the data
dct <- list(TGI = numeric(), Height = numeric(), Shade = numeric(), RealSolar = numeric(), Skyview = numeric(), Set = character())

# Directory where the TIFF files are located (adjust the path accordingly)
flight_files <- list.files('/../Example data/Input data', full.names = TRUE)

#keep only folders included in this paper
base_names <- unlist(strsplit(basename(flight_files), "//"))
flight_files = flight_files[base_names %in% c(train_maps, test_maps) ]

# Loop over each flight
for (flight in flight_files) {
  for (i in 1:5) {
    # Generate random indices
    rand_indices <- sample(1:(1024^2), 100)  # Random 100 indices
    
    # Read the data from each corresponding TIFF file using the raster package
    tgi_raster <- raster(paste0(flight, '/TGI_', i, '.tif'))
    height_raster <- raster(paste0(flight, '/height_', i, '.tif'))
    shade_raster <- raster(paste0(flight, '/shade_', i, '.tif'))
    real_solar_raster <- raster(paste0(flight, '/real_solar_', i, '.tif'))
    skyview_raster <- raster(paste0(flight, '/skyview_', i, '.tif'))
    
    # Extract the random samples
    tgi_data <- as.vector(tgi_raster)[rand_indices]
    height_data <- as.vector(height_raster)[rand_indices]
    shade_data <- as.vector(shade_raster)[rand_indices]
    real_solar_data <- as.vector(real_solar_raster)[rand_indices]
    skyview_data <- as.vector(skyview_raster)[rand_indices]
    
    # Append the random samples to the list
    dct$TGI <- c(dct$TGI, tgi_data)
    dct$Height <- c(dct$Height, height_data)
    dct$Shade <- c(dct$Shade, shade_data)
    dct$RealSolar <- c(dct$RealSolar, real_solar_data)
    dct$Skyview <- c(dct$Skyview, skyview_data)
    
    # Assign the dataset label (TrainSet or TestSet)
    name <- basename(flight)
    ds <- if (name %in% train_maps) "TrainSet" else "TestSet"
    dct$Set <- c(dct$Set, rep(ds, 100))
  }
}

# Convert the list to a data frame
df <- as.data.frame(dct)

# Filter for Height < 10 and melt the data for plotting
mdf <- df %>%
  filter(Height < 10) %>%
  pivot_longer(cols = c("TGI", "Shade", "RealSolar", "Skyview", "Height"), names_to = "Type", values_to = "Value")

# mdf now contains the melted data, ready for further analysis or visualization.
library(ggplot2)
library(ggpubr)

# Define more informative labels for each variable, including units
param_labels <- list(
  TGI = "TGI",
  Height = "Height (m)",
  Shade = "Shade (yes/no)",
  RealSolar = "Real Solar Radiation (W/m^2)",
  Skyview = "Skyview Factor (%)"
)

# Calculate percentages for Shade within each group (TrainSet and TestSet)
shade_mdf <- mdf %>%
  filter(Type == "Shade") %>%
  mutate(Value = ifelse(Value > 0, "Yes", "No")) %>%
  group_by(Set, Value) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = 100 * Count / sum(Count))

# Create another data frame for the continuous variables
numeric_mdf <- mdf %>%
  filter(Type != "Shade") %>%
  mutate(Value = as.numeric(Value))

# Initialize a list to store the plots
plot_list <- list()

# Plot for continuous variables using geom_histogram
for (type in unique(numeric_mdf$Type)) {
  p <- ggplot(numeric_mdf %>% filter(Type == type), aes(x = Value, fill = Set)) +
    geom_histogram(aes(y = after_stat(100 * count / sum(count))), bins = 20, alpha = 0.7) +  # Adjust alpha for more transparency
    xlab(param_labels[[type]]) +         # Use descriptive x-axis label
    ylab("Percentage (%)") +             # Set y-axis label to "Percentage (%)"
    theme_minimal() +
    theme(strip.text = element_blank(),  # Remove facet titles
          panel.border = element_rect(colour = "black", fill = NA),  # Add border around each plot
          axis.text = element_text(size = 12),  # Adjust axis text size
          axis.ticks = element_line(size = 1),  # Thicker axis ticks
          panel.grid.major = element_blank(),   # Remove major gridlines
          panel.grid.minor = element_blank()) + # Remove minor gridlines
    scale_fill_manual(values = c("TrainSet" = "blue", "TestSet" = "red"),  # Color for Train and Test sets
                      labels = c("TestSet" = "Test", "TrainSet" = "Train"))  # Custom labels for the legend
  
  plot_list[[type]] <- p
}

# Plot for Shade using geom_bar with pre-calculated percentages
p_shade <- ggplot(shade_mdf, aes(x = Value, y = Percentage, fill = Set)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  xlab(param_labels[["Shade"]]) +         # Use descriptive x-axis label
  ylab("Percentage (%)") +                # Set y-axis label to "Percentage (%)"
  theme_minimal() +
  theme(strip.text = element_blank(),  # Remove facet titles
        panel.border = element_rect(colour = "black", fill = NA),  # Add border around each plot
        axis.text = element_text(size = 12),  # Adjust axis text size
        axis.ticks = element_line(size = 1),  # Thicker axis ticks
        panel.grid.major = element_blank(),   # Remove major gridlines
        panel.grid.minor = element_blank()) + # Remove minor gridlines
  scale_fill_manual(values = c("TrainSet" = "blue", "TestSet" = "red"),  # Color for Train and Test sets
                    labels = c("TestSet" = "Test", "TrainSet" = "Train"))  # Custom labels for the legend
# Add the Shade plot to the plot list
plot_list[["Shade"]] <- p_shade

# Arrange the plots in a grid using ggarrange
final_plot <- ggarrange(plotlist = plot_list, ncol = 3, nrow = ceiling(length(plot_list) / 3), common.legend = TRUE, legend = "top")

# Save the final plot as a JPEG image
jpeg(filename = "Figure_S3.jpg", width = 2700, height = 1800, res = 300)
print(final_plot)
dev.off()

library(shiny)
library(shinythemes)
library(leaflet)
library(bslib)
library(sf)
library(dplyr)
library(raster)
library(terra)
library(caret)
library(ranger)
library(randomForest)
library(viridis)
library(RColorBrewer)
library(ggplot2)
library(plotly)
library(patchwork)

# Read the datasets
MMTemp_19 <- raster("resampled_MMTemp_19.tif")
Rainfall_2020 <- raster("Rainfall_2020.tif")
RHU_2020 <- raster("resampled_RHU_2020.tif")

Cholera_Histo_Cases <- st_read("Cholera_Histo_Cases.shp")

DumpingSites_buf500m <- sf::st_read("DumpingSites_buf500m.shp")

LULC_2019 <- sf::st_read("LULC_2019.shp")

Nairobi_Pop2019 <- sf::read_sf("Nairobi_Pop2019.shp")
Nairobi_Pop2019 <- st_transform(Nairobi_Pop2019, crs = 4326)

Rivers_Bufzone_100m <- sf::read_sf("Rivers-Bufzone_100m.shp")
Rivers_Bufzone_100m <- st_transform(Rivers_Bufzone_100m, crs = 4326)

################################################################################
#population
bins <- c(185777, 189536, 197489, 206564, 210423, 268276, 308854, 434208, 780656, 1615290)
pal <- colorBin("Oranges", domain = Nairobi_Pop2019$Pop_2019, bins = bins)


### Mapping LULC
# Your data
filtered_LULC_2019 <- LULC_2019

# Define a function to assign colors to land use types
assign_color <- function(landuse) {
  switch(landuse,
         "industrial" = "#E41A1C",
         "institutional" = "#377EB8",
         "residential" = "#4DAF4A",
         "commercial" = "#984EA3",
         "transportation" = "#FF7F00",
         "res_slum" = "#FFFF33",
         "#808080" # Default color for missing or unknown land use types
  )
}

# Add a new column with colors for each land use type
filtered_LULC_2019$color <- sapply(filtered_LULC_2019$LANDUSE, assign_color)


############
#Mapping Temperature
pal0 = colorNumeric(c("#7f007f", "#0000ff",  "#007fff", "#00ffff", "#00bf00", "#7fdf00",
                      "#ffff00", "#ff7f00", "#ff3f00", "#ff0000", "#bf0000"), values(MMTemp_19),  na.color = "transparent")


#mapping Rainfall
# Mapping Rainfall
pal1 <- colorNumeric(viridis(40), values(Rainfall_2020), na.color = "transparent")


####Mapping Humidity
pal2 <- colorNumeric(c("RdYlBu"), values(RHU_2020 ),
                     na.color = "transparent")


####Mapping Distance to rivers in leaflet
bins1 <- c(100, 200, 300, 400)
pal3 <- colorBin("YlGnBu", domain = Rivers_Bufzone_100m$distance, bins = bins1)


###Mapping Distance from dumping sites
bins2 <- c(500, 1000, 1500, 2000)
pal4 <- colorBin("YlOrRd", domain = DumpingSites_buf500m$distance, bins = bins2)

IIR <- read.csv("Data2.csv")

# Convert CFR values to numeric
IIR$CFR_2009 <- as.numeric(sub("%", "", IIR$CFR_2009)) / 100
IIR$CFR_2015 <- as.numeric(sub("%", "", IIR$CFR_2015)) / 100
IIR$CFR_2019 <- as.numeric(sub("%", "", IIR$CFR_2019)) / 100

# Create a new data frame for plotting
to_plot <- data.frame(
  x = rep(IIR$subcounty, 3),
  year = rep(c("2009", "2015", "2019"), each = nrow(IIR)),
  cases = c(IIR$Number_of_Cases_2009, IIR$Number_of_Cases_2015, IIR$Number_of_Cases_2019),
  incidence_rate = c(IIR$Incidence_Rate_2009, IIR$Incidence_Rate_2015, IIR$Incidence_Rate_2019),
  cfr = c(IIR$CFR_2009, IIR$CFR_2015, IIR$CFR_2019),
  dths = c(IIR$Deaths_2009, IIR$Deaths_2015, IIR$Deaths_2019)
)

##########################################################################
###Prediction
# Define the target CRS (UTM Zone 37 South)
target_crs <- "+proj=utm +zone=37 +south +datum=WGS84 +units=m +no_defs"

# Load your raster datasets
# Load the Cholera_KDE raster
Cholera_KDE1 <- raster("Cholera_kde_m_arc.tif")
# Get the minimum and maximum values of the raster
min_val <- minValue(Cholera_KDE1)
max_val <- maxValue(Cholera_KDE1)
# Rescale the raster values to the range 0 to 200
Cholera_KDE <- (Cholera_KDE1 - min_val) / (max_val - min_val) * 200

lulc_rasterr <- raster("resampled_lulc.tif")
pop_rasterr <- raster("resampled_population.tif")
Temp_19 <- raster("resampled_temperature.tif")
Rainfall_20 <- raster("resampled_rainfall.tif")
RHU_20 <- raster("resampled_relative_humidity.tif")

# Combine predictors into a raster stack
rr <- stack(Cholera_KDE, pop_rasterr, lulc_rasterr, Temp_19, Rainfall_20, RHU_20)

# Read cholera case locations
points <- st_read("Cholera_Histo_Cases.shp")
# Reproject to the target CRS
points <- st_transform(points, crs = target_crs)

# Extract values from rasters at cholera case locations
cholera_data <- data.frame(
  Cholera_KDE = extract(Cholera_KDE, points),
  LULC = extract(lulc_rasterr, points),
  Population = extract(pop_rasterr, points),
  MeanTemp_19 = extract(Temp_19, points),
  Rainfall = extract(Rainfall_20, points),
  RHU = extract(RHU_20, points)
)

# Rename the layer to match the variable name in the random forest model
names(rr)[names(rr) == "Cholera_kde_m_arc"] <- "Cholera_KDE"
names(rr)[names(rr) == "resampled_lulc"] <- "LULC"
names(rr)[names(rr) == "resampled_population"] <- "Population"
names(rr)[names(rr) == "resampled_temperature"] <- "MeanTemp_19"
names(rr)[names(rr) == "resampled_rainfall"] <- "Rainfall"
names(rr)[names(rr) == "resampled_relative_humidity"] <- "RHU"

# Prepare response variable (cholera density)
cholera_density <- extract(Cholera_KDE, points)

# Combine predictors and response into a data frame
cholera_df <- cbind(cholera_data, Cholera_Density = cholera_density)

# Split data into training and validation sets
set.seed(123)
train_indices <- sample(1:nrow(cholera_df), 0.7 * nrow(cholera_df))
train_data <- cholera_df[train_indices, ]
validation_data <- cholera_df[-train_indices, ]

# Train Random Forest model (Regression)
rf_model <- randomForest(
  Cholera_Density ~ .,
  data = train_data,
  ntree = 500,  # Adjust as needed
  mtry = 3,     # Number of variables to consider at each split
  importance = TRUE  # Calculate variable importance
)

# Make predictions on validation set
validation_preds <- predict(rf_model, newdata = validation_data)

# Evaluate model performance
# Mean Absolute Error (MAE)
mae <- mean(abs(validation_preds - validation_data$Cholera_Density))

# Adjusted R-squared
rsquared <- cor(validation_preds, validation_data$Cholera_Density)^2
n <- nrow(validation_data)
p <- ncol(validation_data) - 1
adjusted_rsquared <- 1 - (1 - rsquared) * ((n - 1) / (n - p - 1))

# Feature importance
var_importance <- importance(rf_model)

# Save the random forest model
#saveRDS(rf_model, file = "rf_model_regression.rds")

# Spatial prediction
# Predict cholera density across the entire study area
#Included inthe server part
#predicted_cholera <- predict(rr, rf_model)

# Define a custom color palette
custom_palette <- c("#7f007f", "#0000ff", "#007fff", "#00ffff", "#00bf00", "#7fdf00",
                    "#ffff00", "#ff7f00", "#ff3f00", "#ff0000", "#bf0000")

# Convert colors to RGB
custom_palette_rgb <- col2rgb(custom_palette)

# Convert RGB to HEX
custom_palette_hex <- rgb(custom_palette_rgb[1, ], custom_palette_rgb[2, ], custom_palette_rgb[3, ], maxColorValue = 255)

#######################################################################
# Define UI 
ui <- shinyUI(fluidPage(
  theme = bs_theme(bg = "#F0F7FF", 
                   fg = "#2A3F54",
                   primary = "#557",
                   font_scale = 0.7, bootswatch = "flatly"),
  
  titlePanel("Nairobi Cholera Hotspots Prediction Model"),
  
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("environmental_variables",
                         label = "Environmental Variables:",
                         choices = c("Temperature",
                                     "Rainfall",
                                     "Relative Humidity",
                                     "Distance from River",
                                     "Distance from Dumping Site"),
                         selected = "Temperature",  # Set "Temperature" as the default selected option
                         inline = TRUE),
      
      checkboxGroupInput("socio_ecological_variables",
                         label = "Socio-Ecological Variables:",
                         choices = c("Population",
                                     "Land Use Land Cover"),
                         inline = TRUE),
      
      checkboxGroupInput("historical_hotspot_variable",
                         label = "Historical Cholera Hotspots:",
                         choices = c("Historical Cholera Hotspots"),
                         inline = TRUE),
      
      # Styling the action button
      tags$div(
        actionButton("predict_hotspots_button", 
                     icon = icon("search"),    # Add a search icon
                     "Predict Cholera Hotspot",
                     style = "color: #fff; background-color: #ff7f0e; border-color: #ff7f0e; padding: 10px; border-radius: 5px;")
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Maps",
                 leafletOutput("map")
        ),
        tabPanel("Incidence_Rate_Plots",
                 plotOutput("incidence_rate_graph")
        ),
        tabPanel("CFR_Plots",
                 plotOutput("case_fatality_rate_graph"))
        
      )
    )
  )
))

# Define server logic
server <- function(input, output) {
  
  observeEvent(input$historical_hotspot_variable, {
    if ("Historical Cholera Hotspots" %in% input$historical_hotspot_variable) {
      output$map <- renderLeaflet({
        leaflet() %>%
          addTiles() %>%
          addMarkers(data = Cholera_Histo_Cases,
                     popup = ~glue::glue("<b>Local name:</b> {Local_name}"))
        
      })
    }
  })
  
  
  observeEvent(input$environmental_variables, {
    if ("Temperature" %in% input$environmental_variables) {
      output$map <- renderLeaflet({
        leaflet() %>%
          addTiles() %>%
          addRasterImage(MMTemp_19, colors = pal0, opacity = 0.8) %>%
          addLegend("bottomright" ,pal = pal0, values = values(MMTemp_19),
                    opacity = 0.5,
                    title = "Temperature (°C)") 
      })
    }
    if ("Rainfall" %in% input$environmental_variables) {
      output$map <- renderLeaflet({
        leaflet() %>%
          addTiles() %>%
          addRasterImage(Rainfall_2020, colors = pal1, opacity = 0.8) %>%
          addLegend("bottomright", pal = pal1, values = values(Rainfall_2020),
                    opacity = 0.8,
                    title = "Rainfall (mm)")
      })
    }
    if ("Relative Humidity" %in% input$environmental_variables) {
      output$map <- renderLeaflet({
        leaflet() %>%
          addTiles() %>%
          addRasterImage(RHU_2020, colors = pal2, opacity = 0.8) %>%
          addLegend("bottomright" ,pal = pal2, values = values(RHU_2020),
                    opacity = 0.8,
                    title = "Relative Humidity")
      })
    }
    if ("Distance from River" %in% input$environmental_variables) {
      output$map <- renderLeaflet({
        leaflet() %>%
          addTiles() %>%
          addPolygons(data = Rivers_Bufzone_100m, 
                      stroke = TRUE,
                      color = "rgba(0, 0, 0, 0.5)", 
                      opacity = 1,
                      weight = 1,
                      fillOpacity = 0.7, smoothFactor = 0.9,
                      fillColor = ~pal3(distance),
                      popup = ~glue::glue("<b>distance :</b> {distance}")) %>%  
          addLegend("bottomright", pal = pal3, values = Rivers_Bufzone_100m$distance, 
                    opacity = 0.9,
                    title = "Distance From River (m)")
      })
    }
    if ("Distance from Dumping Site" %in% input$environmental_variables) {
      output$map <- renderLeaflet({
        leaflet() %>%
          addTiles() %>%
          addPolygons(data = DumpingSites_buf500m, 
                      stroke = TRUE,
                      color = "rgba(0, 0, 0, 0.5)", 
                      opacity = 1,
                      weight = 1,
                      fillOpacity = 0.7, smoothFactor = 0.9,
                      fillColor = ~pal4(distance),
                      popup = ~glue::glue("<b>Distance :</b> {distance}")) %>%  
          addLegend("bottomright", pal = pal4, values = DumpingSites_buf500m$distance, 
                    opacity = 0.9,
                    title = "Distance From Dumping Site (m)")
      })
    }
  })
  
  observeEvent(input$socio_ecological_variables, {
    if ("Population" %in% input$socio_ecological_variables) {
      output$map <- renderLeaflet({
        
        leaflet() %>%
          addTiles() %>%
          addPolygons(data = Nairobi_Pop2019, 
                      stroke = TRUE,
                      color = "rgba(0, 0, 0, 0.5)",
                      opacity = 1,
                      weight = 1,
                      fillOpacity = 0.7, smoothFactor = 0.5,
                      fillColor = ~pal(Pop_2019),
                      popup = ~glue::glue("<b>Subcounty:</b> {subcounty}<br><b>Population:</b> {Pop_2019}")) %>%
          addLegend("bottomright", pal = pal, values = Nairobi_Pop2019$Pop_2019, 
                    title = "Population")
        
        
      })
    }
    if ("Land Use Land Cover" %in% input$socio_ecological_variables) {
      output$map <- renderLeaflet({
        leaflet() %>%
          addTiles() %>%
          addPolygons(data = filtered_LULC_2019, 
                      stroke = TRUE,
                      color = "rgba(0, 0, 0, 0.5)",  # Adjust the transparency of the boundary
                      opacity = 1,
                      weight = 1,
                      fillOpacity = 0.8,
                      fillColor = ~color,
                      popup = ~glue::glue("<b>Land Use:</b> {LANDUSE}")) %>%
          addLegend("bottomright", 
                    colors = unique(filtered_LULC_2019$color),
                    labels = unique(filtered_LULC_2019$LANDUSE),
                    title = "Land Use Types")
      })
    }
  })
  
  # Plot the bar chart with separate legends for bar plots and line graph
  output$incidence_rate_graph <- renderPlot({    
    ggplot(to_plot, aes(x = x, y = cases, fill = year)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_line(aes(y = incidence_rate * 10, color = year, group = year), size = 1) +
      geom_point(aes(y = incidence_rate * 10, color = year), size = 3) +
      labs(title = "Number of Cases and Incidence Rate per Year by Subcounty",
           x = "Subcounty", y = "Number of cases") +
      scale_fill_manual(values = c("2009" = "lightblue", "2015" = "lightgray", "2019" = "salmon"),
                        name = "Number of cases/Year") +   # Legend for bar plots
      scale_color_manual(values = c("2009" = "lightblue", "2015" = "lightgrey", "2019" = "salmon"),
                         name = "Incidence Rate (IR)") +  # Legend for line graph
      theme_minimal() +
      scale_y_continuous(sec.axis = sec_axis(~./10, name = "Incidence Rate")) +
      guides(color = guide_legend(title = "Incidence Rate (IR)"),           # Legend for line graph
             fill = guide_legend(title = "Number of cases/Year"))  # Legend for bar plots
  })
  
  # Plot the Case Fatality Rate (CFR) graph
  output$case_fatality_rate_graph <- renderPlot({
    # Plot the CFR and number of deaths on the same graph
    ggplot(IIR, aes(x = subcounty)) +
      geom_bar(aes(y = Deaths_2009, fill = "2009"), stat = "identity", position = "dodge") +
      geom_bar(aes(y = Deaths_2015, fill = "2015"), stat = "identity", position = "dodge") +
      geom_bar(aes(y = Deaths_2019, fill = "2019"), stat = "identity", position = "dodge") +
      geom_line(aes(y = CFR_2009 * 100, color = "2009"), size = 1, group = 1) +
      geom_line(aes(y = CFR_2015 * 100, color = "2015"), size = 1, group = 1) +
      geom_line(aes(y = CFR_2019 * 100, color = "2019"), size = 1, group = 1) +
      geom_text(aes(y = Deaths_2009, label = Deaths_2009), vjust = -0.5, hjust = 0.5, size = 3) + # Add labels for Deaths_2009
      geom_text(aes(y = Deaths_2015, label = Deaths_2015), vjust = -0.5, hjust = 0.5, size = 3) + # Add labels for Deaths_2015
      geom_text(aes(y = Deaths_2019, label = Deaths_2019), vjust = -0.5, hjust = 0.5, size = 3) + # Add labels for Deaths_2019
      geom_text(aes(y = CFR_2009 * 100, label = sprintf("%.2f", CFR_2009 * 100)), vjust = 1.5, hjust = 0.5, size = 3) + # Add labels for CFR_2009
      geom_text(aes(y = CFR_2015 * 100, label = sprintf("%.2f", CFR_2015 * 100)), vjust = 1.5, hjust = 0.5, size = 3) + # Add labels for CFR_2015
      geom_text(aes(y = CFR_2019 * 100, label = sprintf("%.2f", CFR_2019 * 100)), vjust = 1.5, hjust = 0.5, size = 3) + # Add labels for CFR_2019
      labs(title = "Case Fatality Rate (CFR) per Year by Subcounty",
           x = "Subcounty", y = "Number of Deaths / CFR (%)") +
      scale_fill_manual(name = "Number of Deaths/Year", values = c("2009" = "blue", "2015" = "red", "2019" = "green")) +
      scale_color_manual(name = "CFR/Year(%)", values = c("2009" = "blue", "2015" = "red", "2019" = "green")) +
      theme_minimal()
  })
  # Add an observer for the prediction button click event
  observeEvent(input$predict_hotspots_button, {
    # Perform prediction using the trained random forest model
    predicted_cholera <- predict(rr, rf_model)  
    
    # Use the custom palette in your plot
    pal_custom <- colorNumeric(custom_palette_hex, values(predicted_cholera), na.color = "transparent")
    
    # Render the leaflet map with the predicted cholera density
    output$map <- renderLeaflet({
      leaflet() %>%
        addTiles() %>%
        addRasterImage(predicted_cholera, colors = pal_custom, opacity = 0.7) %>%
        addLegend("bottomright" ,pal = pal_custom, values = values(predicted_cholera),
                  opacity = 0.7,
                  title = "Predicted Cholera Hotspot")
    })
  })
}

# Run application
shinyApp(ui = ui, server = server)


##### Create species occurrence data ####
# Import Amazonian shapefile
vectorShape <- readOGR("Maps/wgs1984.shp")


# Convert region shapefile into a raster (make grid of 1x1 cells)
r <- raster(nrow = 180, ncol = 360)
regionRaster <- rasterize(vectorShape, r)

# Extract coordinates from the centroid of each cell
regionCoords <- coordinates(regionRaster[[1]]) # Get all cells in the globe
regionCoords <- regionCoords[!is.na(regionRaster[[1]]@data@values), ] # Filter only cells in region

### Extract occurrence data for each cell
# Folder with the shapefiles of all passeriformes birds
shapeFolder <- "Maps/Passeriformes/"

# List bird shapefiles
Vectorfiles <- list.files(shapeFolder, pattern = "*.shp", recursive = TRUE, full.names = TRUE)
shapenames <- gsub("(.*?)/(.*?)_(.*?).shp","\\2", Vectorfiles)

makePAM <- function(Vectorfiles){
  # Set counter to zero
  n <- 0
  # Crete empty matrix
  community <- NULL
  for(i in Vectorfiles){
    vectorShape <- readOGR(paste(i, sep = ""))  
    tmpRaster <- rasterize(vectorShape, regionRaster,mask=TRUE)  
    values <- extract(tmpRaster, regionCoords)
    community <- cbind(community,(!is.na(values))*1)
    n <- n + 1
    print(n)
    }
  colnames(community) <- shapenames
}

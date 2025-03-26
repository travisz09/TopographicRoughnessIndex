# Topographic Roughness Index
# Travis Zalesky
# 3/26/25

 # Calculate the Topographic Roughness Index of a digital elevation model

# Import libraries
# Requires terra package. If not already installed run install.packages('terra') to install
library(tera)

# TRI function
calcTRI <- function(r) {
  start <- Sys.time()
  cat(paste0('TRI starting at ', format(start, "%X"), '...\n'))
  cat('\tExtracting Coordinate Reference System...\n')  # Progress Statement
  # Get crs of input rast
  crs <- crs(r)
  cat('\tExtracting Focal Values...\n')  # Progress Statement
  # Get matrix of focal values
  # Neighborhood = 3x3, Center = 5
  focal <- focalValues(r)
  # Get center values as 1-col matrix
  center <- focal[,5]

  cat('\tCalculating Difference Squared...\n')  # Progress Statement
  # Calculate difference squared for all cells in neighborhood (relative to center)
  # Returns a nested list
  sqrDiff <- lapply(seq_len(ncol(focal)), function(i) {
    (focal[,i] - center)^2
  })
  # Convert list back into matrix
  sqrDiff <- matrix(unlist(sqrDiff), ncol = ncol(focal), byrow = F)

  cat('\tCalculating Topographic Roughness Index...\n')  # Progress Statement
  # Reduce matrix to sum of rows
  tri <- rowSums(sqrDiff, na.rm = T)^(1/2)  # sqrt row sums
  cat('\tFinishing Up...\n')  # Progress Statement
  # Convert tri list to spatRast matching input rast
  m <- matrix(tri, ncol = ncol(r), byrow = T)
  tri <- rast(m, crs = crs)

  # Use input raster to preserve raster structure
  # Swap tri values into r
  values(r) <- values(tri)
  
  end <- Sys.time()
  duration <- end - start
  cat(paste0('TRI ended at ', format(end, "%X"), '. Duration, ', duration, ' seconds'))
  
  return(r)
}

# Some testing data
testRast <- rast(
  matrix(c(0,0,0,0,1, 0,0,0,0,0, 0,0,0,1,0, 1,1,5,0,1, 1,0,0,1,1), nrow = 5), 
  crs = 'EPSG:32612'
) 
plot(testRast)

# Test run
calcTRI(testRast)
plot(calcTRI(testRast))

# Very large DEM files may fill up available memory and result in a failure. If this happens then use batch processing to compute TRI in smaller chunks.

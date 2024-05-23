# Load libraries 
require(pacman)
pacman::p_load(tidyverse, sf, fs, rgbif, dplyr, terra, raster, dismo)

# Prepare data input

# Bioclimatic variable calculation
prec <- list.files(path = "raster/prec/", pattern='.tif', 
                       all.files=TRUE, full.names=TRUE)
prec <- stack(prec)

tmax <- list.files(path = "raster/tmax/", pattern='.tif', 
                   all.files=TRUE, full.names=TRUE)
tmax <- stack(tmax)

tmin <- list.files(path = "raster/tmin/", pattern='.tif', 
                   all.files=TRUE, full.names=TRUE)
tmin <- stack(tmin)

bio <- biovars(prec, tmin, tmax)
stackSave(bio, "raster/bioclim")

# Load environmental layers 

elev <- raster(paste0("raster/elevation.tif"))
popden <- raster(paste0("raster/popden.tif"))
urbacc <- raster(paste0("raster/urbacc.tif"))
LC <- raster(paste0("raster/landc.tif"))

### Resampled to 2.5m

LC <- resample(LC, elev, method='ngb')

# Crop rasters

NEA <- read_sf("shp/Sudamérica.shp")

bio <- terra::crop(bio, NEA)
elev <- terra::crop(elev, NEA)
popden <- terra::crop(popden, NEA)
urbacc <- terra::crop(urbacc, NEA)
LC <- terra::crop(LC, NEA)

# Mask

bio <- terra::mask(bio, NEA)
elev <- terra::mask(elev, NEA)
popden <- terra::mask(popden, NEA)
urbacc <- terra::mask(urbacc, NEA)
LC <- terra::mask(LC, NEA)

# Separate Bioclimatic variables

bio1 <- bio[[1]]
bio2 <- bio[[2]]
bio3 <- bio[[3]]
bio4 <- bio[[4]]
bio5 <- bio[[5]]
bio6 <- bio[[6]]
bio7 <- bio[[7]]
bio8 <- bio[[8]]
bio9 <- bio[[9]]
bio10 <- bio[[10]]
bio11 <- bio[[11]]
bio12 <- bio[[12]]
bio13 <- bio[[13]]
bio14 <- bio[[14]]
bio15 <- bio[[15]]
bio16 <- bio[[16]]
bio17 <- bio[[17]]
bio18 <- bio[[18]]
bio19 <- bio[[19]]

# Write rasters

terra::writeRaster(x = bio1, filename = 'raster/bio1.tif')
terra::writeRaster(x = bio2, filename = 'raster/bio2.tif')
terra::writeRaster(x = bio3, filename = 'raster/bio3.tif')
terra::writeRaster(x = bio4, filename = 'raster/io4.tif')
terra::writeRaster(x = bio5, filename = 'raster/bio5.tif')
terra::writeRaster(x = bio6, filename = 'raster/bio6.tif')
terra::writeRaster(x = bio7, filename = 'raster/bio7.tif')
terra::writeRaster(x = bio8, filename = 'raster/bio8.tif')
terra::writeRaster(x = bio9, filename = 'raster/bio9.tif')
terra::writeRaster(x = bio10, filename = 'raster/bio10.tif')
terra::writeRaster(x = bio11, filename = 'raster/bio11.tif')
terra::writeRaster(x = bio12, filename = 'raster/bio12.tif')
terra::writeRaster(x = bio13, filename = 'raster/bio13.tif')
terra::writeRaster(x = bio14, filename = 'raster/bio14.tif')
terra::writeRaster(x = bio15, filename = 'raster/bio15.tif')
terra::writeRaster(x = bio16, filename = 'raster/bio16.tif')
terra::writeRaster(x = bio17, filename = 'raster/bio17.tif')
terra::writeRaster(x = bio18, filename = 'raster/bio18.tif')
terra::writeRaster(x = bio19, filename = 'raster/bio19.tif')
terra::writeRaster(x = elev, filename = 'raster/elevation.tif')
terra::writeRaster(x = LC, filename = 'raster/land_cover.tif')
terra::writeRaster(x = popden, filename = 'raster/population_density.tif')
terra::writeRaster(x = urbacc, filename = 'raster/urban_accessibility.tif')

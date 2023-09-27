# Load libraries 
require(pacman)
pacman::p_load(geodata, terra, glue, rgeos, gtools, ggspatial, sf, 
               colourpicker, tidyverse, fs, rnaturalearthdata, 
               rnaturalearth, glue, outliers, stringr, raster)

### Load environmental data 

bio1 <- raster(paste0("raster/bio1.tif"))
bio2 <- raster(paste0("raster/bio2.tif"))
bio3 <- raster(paste0("raster/bio3.tif"))
bio4 <- raster(paste0("raster/bio4.tif"))
bio5 <- raster(paste0("raster/bio5.tif"))
bio6 <- raster(paste0("raster/bio6.tif"))
bio7 <- raster(paste0("raster/bio7.tif"))
bio8 <- raster(paste0("raster/bio8.tif"))
bio9 <- raster(paste0("raster/bio9.tif"))
bio10 <- raster(paste0("raster/bio10.tif"))
bio11 <- raster(paste0("raster/bio11.tif"))
bio12 <- raster(paste0("raster/bio12.tif"))
bio13 <- raster(paste0("raster/bio13.tif"))
bio14 <- raster(paste0("raster/bio14.tif"))
bio15 <- raster(paste0("raster/bio15.tif"))
bio16 <- raster(paste0("raster/bio16.tif"))
bio17 <- raster(paste0("raster/bio17.tif"))
bio18 <- raster(paste0("raster/bio18.tif"))
bio19 <- raster(paste0("raster/bio19.tif"))
elev <- raster(paste0("raster/elevation.tif"))
PopDen <- raster(paste0("raster/populationdensity.tif"))
UrbAcc <- raster(paste0("raster/urbanaccessibility.tif"))
LC <- raster(paste0("raster/landcover.tif"))

NEA <- read_sf("shp/Sudamérica.shp")

# stacking the variables to process them at one go 
var <- raster::stack(bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10,
                     bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19,
                     elev, LC, PopDen, UrbAcc)

### Load occurrence data 

occ_albo_raw <- read.csv("data/albopictus_occ.csv") #10536 obs.
occ_albo_raw <- select(occ_albo_raw, -scientific_name, -layer, -path)

occ_aeg_raw <- read.csv("data/aegypti_occ.csv") #64110 obs.
occ_aeg_raw <- select(occ_aeg_raw, -scientific_name, -layer, -path)


### Clean occurrence data

# remove erroneous coordinates, where either the latitude or longitude is missing

occ_albo_clean <- subset(occ_albo_raw,(!is.na(latitude))&(!is.na(longitude)))
cat(nrow(occ_albo_raw)-nrow(occ_albo_clean), "records are removed") #0 records are removed

occ_aeg_clean <- subset(occ_aeg_raw,(!is.na(latitude))&(!is.na(longitude)))
cat(nrow(occ_aeg_raw)-nrow(occ_aeg_clean), "records are removed") #0 records are removed

# remove duplicated data based on latitude and longitude

dups_albo <- duplicated(occ_albo_clean[c("latitude","longitude")])
occ_albo <- occ_albo_clean[!dups_albo,]
cat(nrow(occ_albo_clean)-nrow(occ_albo), "records are removed") #7031 records are removed


dups_aeg <- duplicated(occ_aeg_clean[c("latitude","longitude")])
occ_aeg <- occ_aeg_clean[!dups_aeg,]
cat(nrow(occ_aeg_clean)-nrow(occ_aeg), "records are removed") #47830 records are removed

# only one occurrence point per pixel

coordinates(occ_albo) <- ~ longitude + latitude
cells1 <- cellFromXY(var[[1]],occ_albo)
dups1 <- duplicated(cells1)
occ_final1 <- occ_albo[!dups1,]
cat(nrow(occ_albo)-nrow(occ_final1), "records are removed") #69 records are removed

coordinates(occ_aeg) <- ~ longitude + latitude
cells2 <- cellFromXY(var[[1]],occ_aeg)
dups2 <- duplicated(cells2)
occ_final2 <- occ_aeg[!dups2,]
cat(nrow(occ_aeg)-nrow(occ_final2), "records are removed") #10844 records are removed

# Write the result as table

write.csv(occ_final1, 'data/albopictus_final.csv', row.names = FALSE) 
write.csv(occ_final2, 'data/aegypti_final.csv', row.names = FALSE) 

albo <- read.csv("data/albopictus_final.csv")
aeg <- read.csv("data/aegypti_final.csv")

# Extract the values for the presences 

vles1 <- terra::extract(var, albo[,1:2])
vles1 <- as_tibble(cbind(albo, vles1))

vles2 <- terra::extract(var, aeg[,1:2])
vles2 <- as_tibble(cbind(aeg, vles2))

# Add environmental values

pnts1 <- cbind(vles1[,1:2], raster::extract(var, vles1[,1:2]))
names(pnts1) <- c("x", "y", "bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19","elev","LC", "popden","urbacc")
pnts1 <- as_tibble(pnts1)

write.csv(pnts1, 'data/albopictus_clean.csv', row.names = FALSE)


pnts2 <- cbind(vles2[,1:2], raster::extract(var, vles2[,1:2]))
names(pnts2) <- c("x", "y", "bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19","elev","LC", "popden","urbacc")
pnts2 <- as_tibble(pnts2)

write.csv(pnts2, 'data/aegypti_clean.csv', row.names = FALSE)

#albopictus_clean: 3411 obs.
#aegypti_clean: 5051 obs.

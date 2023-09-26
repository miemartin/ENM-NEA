library(raster)
library(rstudioapi)
library(dismo)

# set working directory to file location

setwd("C:/Users/Giga/Desktop/aedes potential maps/South America/current climate/")


# list all binary rasters created in Threshold

aeg1 <- raster(paste0("binary maps/current_all_aeg_p10.tif"))
aeg2 <- raster(paste0("binary maps/current_bio_aeg_p10.tif"))
albo1 <- raster(paste0("binary maps/current_all_albo_p10.tif"))
albo2 <- raster(paste0("binary maps/current_bio_albo_p10.tif"))

# create empty matrix to be filled with suitable area values 

suitability <- matrix(data = NA, nrow = 4, ncol = 3)

# create object that lists the frequency of 0 and 1 pixels 

f <- freq(aeg1)

# calculate percentage of area suitable by dividing number of 1 pixels by the number of 0 + 1 pixels

area1 <- round( ( f[2,2] / ( f[1,2] + f[2,2] ) ) * 100, 2)

# write results and label to the matrix 

suitability[1,1] <- paste0("Ae. aegypti")
suitability[1,2] <- paste0("all")
suitability[1,3] <- area1

f <- freq(aeg2)
area2 <- round( ( f[2,2] / ( f[1,2] + f[2,2] ) ) * 100, 2)
suitability[2,1] <- paste0("Ae. aegypti")
suitability[2,2] <- paste0("bioc")
suitability[2,3] <- area2

f <- freq(albo1)
area3 <- round( ( f[2,2] / ( f[1,2] + f[2,2] ) ) * 100, 2)
suitability[3,1] <- paste0("Ae. albopictus")
suitability[3,2] <- paste0("all")
suitability[3,3] <- area3

f <- freq(albo2)
area4 <- round( ( f[2,2] / ( f[1,2] + f[2,2] ) ) * 100, 2)
suitability[4,1] <- paste0("Ae. albopictus")
suitability[4,2] <- paste0("bioc")
suitability[4,3] <- area4


# format table and write .csv to the Binary Maps folder 

colnames(suitability) <- c("Specie","Variables","SuitableArea")
write.csv(suitability, file = "SuitableArea.csv")


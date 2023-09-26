library(spocc)
library(spThin)
library(dismo)
library(rgeos)
library(raster)
library(rgdal)
library(ENMeval)
library(wallace)
library(sf)
library(stats)

### Selected best model: fc.LQ_rm.2

### Partition occurrence data

groups_Aa2 <- part_partitionOccs(
  occs = occs_Aa ,
  bg =  bgSample_Aa, 
  method = "rand",
  kfolds = 25) 

# Running modelling

model_Aa_best <- model_maxent(
  occs = occs_Aa,
  bg = bgEnvsVals_Aa,
  user.grp = groups_Aa2, 
  bgMsk = bgMask_Aa,
  rms = c(2, 3), 
  rmsStep =  1.5,
  fcs = c("LQ"),
  clampSel = TRUE,
  algMaxent = "maxent.jar",
  catEnvs = "landcover",
  parallel = FALSE)

evalTbl_b <- model_Aa_best@results
write.csv(evalTbl_b, "best_model_aegypti_all.csv", row.names=FALSE)
evalMods <- model_Aa_best@models
evalPred <- model_Aa_best@predictions
evalVar <- model_Aa_best@variable.importance
write.csv(evalVar, "var.imp_aegypti_all.csv", row.names=FALSE)

# view response curves for environmental variables

response(evalMods[["fc.LQ_rm.2"]])


### Transfer model
#Transfering the model to a new area

NEA <- read_sf("shp/NEA2.shp")

# Generate a transfer of the model to the desired area
xfer_area_Aa <- xfer_area(
  evalOut = model_Aa_best,
  curModel = "fc.LQ_rm.2",
  envs = envs_Aa , 
  outputType = "cloglog",
  alg = "maxent.jar",
  clamp = TRUE,
  xfExt = NEA) 

# store the cropped transfer variables

xferExt_Aa <- xfer_area_Aa$xferExt

plot(xfer_area_Aa$xferArea)
writeRaster(xfer_area_Aa$xferArea, "current_all_aeg.tif")


### Binary map

#Generate a map of the MaxEnt generated model with a “p10” threshold rule.
#10 percentile training presence: 0.4684694

# Select current model and obtain raster prediction

m_Aa <- model_Aa_best@models[["fc.LQ_rm.2"]] 
predSel_Aa <- dismo::predict(m_Aa, bgMask_Aa,
                             args = c(paste0("outputformat=", "cloglog"), 
                                      paste0("doclamp=", tolower(as.character(TRUE)))), 
                             na.rm = TRUE)

# determine the threshold based on the current prediction

occPredVals_Aa <- raster::extract(predSel_Aa, occs_xy_Aa)

# define probability of quantile based on selected threshold

thresProb_Aa <- switch("p10", "mtp" = 0, "p10" = 0.1, "qtp" = 0)

# define threshold value

thres_Aa <- stats::quantile(occPredVals_Aa, probs = thresProb_Aa, na.rm = TRUE)

# applied selected threshold

predSel_Aa <- predSel_Aa > thres_Aa


### Transfer model

# Transfering the model to a new area using a “p10” threshold

# add threshold if specified 

xfer_area_Aa <- xfer_area_Aa$xferArea > thres_Aa

# plot and save the transfer raster

xferExt_Aa <- xfer_area_Aa

# plot and save the transfer raster

plot(xfer_area_Aa)
writeRaster(xfer_area_Aa, "current_all_aeg_p10.tif")

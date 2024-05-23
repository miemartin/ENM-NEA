# Load libraries
require(pacman)
pacman::p_load(spocc, spThin, dismo, raster, rgeos, rgdal, ENMeval, wallace, 
               sf, stats, SDMtune)


### Analysis for Aedes aegypti (Aa)

# Load occurrence data

occs_path <- "data"
occs_path <- file.path(occs_path, "wallace_aeg.csv")

# get a list of species occurrence data

userOccs_Aa <- occs_userOccs(
  txtPath = occs_path, 
  txtName = "wallace_aeg.csv", 
  txtSep = ";", 
  txtDec = ".")
occs_Aa <- userOccs_Aa$Aedes_aegypti$cleaned


# Load environmental data

dir_envs_Aa <- "raster"
envs_path <- file.path(dir_envs_Aa, c("bio2.tif", "bio4.tif", "bio8.tif", "bio9.tif", "bio15.tif", "bio17.tif", "bio18.tif", "bio19.tif", "elevation.tif", "land_cover.tif", "population_density.tif", "urban_accessibility.tif"))

# Create environmental object 

envs_Aa <- envs_userEnvs(
  rasPath = envs_path,
  rasName = c("bio2.tif", "bio4.tif", "bio8.tif", "bio9.tif", "bio15.tif", "bio17.tif", "bio18.tif", "bio19.tif", "elevation.tif", "land_cover.tif", "population_density.tif", "urban_accessibility.tif"),
  doBrick = FALSE)

occs_xy_Aa <- occs_Aa[c("longitude", "latitude")]
occs_vals_Aa <- as.data.frame(raster::extract(envs_Aa, occs_xy_Aa, cellnumbers = TRUE))

# remove duplicated same cell values

occs_Aa <- occs_Aa[!duplicated(occs_vals_Aa[, 1]), ]
occs_vals_Aa <- occs_vals_Aa[!duplicated(occs_vals_Aa[, 1]), -1]

# remove occurrence records with NA environmental values

occs_Aa <- occs_Aa[!(rowSums(is.na(occs_vals_Aa)) >= 1), ]

# also remove variable value rows with NA environmental values

occs_vals_Aa <- na.omit(occs_vals_Aa)

# add columns for env variable values for each occurrence record

occs_Aa <- cbind(occs_Aa, occs_vals_Aa)


# Process environmental data
# Sampling of 10000 background points and corresponding environmental data using a “point buffers” method with a 2 degree buffer.

# Generate background extent 

bgExt_Aa <- penvs_bgExtent(
  occs = occs_Aa,
  bgSel = "point buffers",
  bgBuf = 2)

plot(bgExt_Aa)

# Mask environmental data to provided extent

bgMask_Aa <- penvs_bgMask(
  occs = occs_Aa,
  envs = envs_Aa,
  bgExt = bgExt_Aa)

plot(bgMask_Aa)

# Sample background points from the provided area

bgSample_Aa <- penvs_bgSample(
  occs = occs_Aa,
  bgMask =  bgMask_Aa,
  bgPtsNum = 10000)

# Extract values of environmental layers for each background point

bgEnvsVals_Aa <- as.data.frame(raster::extract(bgMask_Aa,  bgSample_Aa))

# Add extracted values to background points table

bgEnvsVals_Aa <- cbind(scientific_name = paste0("bg_", "Aedes aegypti"), bgSample_Aa,
                       occID = NA, year = NA, institution_code = NA, country = NA,
                       state_province = NA, locality = NA, elevation = NA,
                       record_type = NA, bgEnvsVals_Aa)


# Partition occurrence data
# Partition occurrences and background points for model training and validation using random k-fold, a non-spatial partition method.

groups_Aa <- part_partitionOccs(
  occs = occs_Aa ,
  bg =  bgSample_Aa, 
  method = "rand",
  kfolds = 2) 


# Build and Evaluate Niche Model
# Generating a species distribution model using the maxent.jar algorithm as implemented in ENMeval V2.0 (with clamping = TRUE). 
# For tuning using L, LQ, LQH, LQHP feature classes.

# Run maxent model for the selected species

# Regularization multipliers in the 0.1-1 range increasing by 0.1. 

model_Aa <- model_maxent(
  occs = occs_Aa,
  bg = bgEnvsVals_Aa,
  user.grp = groups_Aa, 
  bgMsk = bgMask_Aa,
  rms = c(0.1, 1), 
  rmsStep =  0.1,
  fcs = c('L', 'LQ', "H", "LQH", "LQHP"),
  clampSel = TRUE,
  algMaxent = "maxent.jar",
  catEnvs = "land_cover",
  parallel = FALSE)

evalTbl <- model_Aa@results #50 models
write.csv(evalTbl, "candidate models/candidate_models_aegypti_all1.csv", row.names=FALSE)

# Regularization multipliers in the 2-6 range increasing by 1. 

model_Aa2 <- model_maxent(
  occs = occs_Aa,
  bg = bgEnvsVals_Aa,
  user.grp = groups_Aa, 
  bgMsk = bgMask_Aa,
  rms = c(2, 6), 
  rmsStep =  1,
  fcs = c('L', 'LQ', "H", "LQH", "LQHP"),
  clampSel = TRUE,
  algMaxent = "maxent.jar",
  catEnvs = "landcover",
  parallel = FALSE)

evalTbl2 <- model_Aa2@results #25 models
write.csv(evalTbl2, "candidate models/candidate_models_aegypti_all2.csv", row.names=FALSE)

# Regularization multipliers 8 and 10

model_Aa3 <- model_maxent(
  occs = occs_Aa,
  bg = bgEnvsVals_Aa,
  user.grp = groups_Aa, 
  bgMsk = bgMask_Aa,
  rms = c(8, 10), 
  rmsStep =  2,
  fcs = c('L', 'LQ', "H", "LQH", "LQHP"),
  clampSel = TRUE,
  algMaxent = "maxent.jar",
  catEnvs = "landcover",
  parallel = FALSE)

evalTbl3 <- model_Aa3@results #10 models
write.csv(evalTbl3, "candidate models/candidate_models_aegypti_all3.csv", row.names=FALSE)


# The same code was used for the model that includes only the bioclimatic variables and also for the Ae. albopictus models. 

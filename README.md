# ENM-NEA
This repository contains the scripts that we have developed to develop Ecological Niche Models and potential geographical distribution of Aedes aegypti and Aedes albopictus in Northeastern Argentina.

The order in which the scripts should be used is as follows:

1a. Environmental Data: This script imports environmental data, assigns projections, and clips variables to the study area.

1b. occurrence data: This script imports the occurrence data and does all the processing for it. It eliminates duplicate data, assigns only one data per pixel and finally creates a .csv file where the values of the environmental variables are obtained for each location.

2. VIF analysis: this is the script where the selection of environmental variables is carried out based on the VIF values and Pearson correlation.

3a. Modeling all var Ae. aegypti: this is the script to perform the calibration of candidate models from the combination of different feature classes and regularization multipliers.

5a. Best model Ae. aegypti all fc.LQ_rm.2: this is the script that runs the best selected model. Applies 25-fold cross validation. The model metrics, the average predictive map, and the response curves of the environmental variables used are obtained. Then the model is transferred to the new study area (in this case, Northeast Argentina) and finally the binary map is obtained by using the “p10” threshold rule.

7. SuitableArea: this is the script where the area (in percentages) that is suitable for the species is calculated from the binary maps obtained previously.
# Code Author - Suraj Reddy Rodda
# This code is generated as part of the study 
# LiDAR-based reference aboveground biomass maps for tropical forests of South Asia and Central Africa
# Submitted to Nature Scientific Data

# Install the BIOMASS package for biomass analysis
#install.packages("BIOMASS") - Install if this package is not available

# Load the required packages
require(BIOMASS)
require(dplyr)  # Required for data manipulation

# Part-1 (Plot-AGB Computation) ----

# Since the plot data in the current study is not shared, 
# we use a sample dataset made available as part of the "BIOMASS" package 
# for computing plot-level AGB (Aboveground Biomass) and its uncertainty.

# Load the field_inv_dataset dataset and display its structure
data("NouraguesHD")
str(NouraguesHD)

# Assign the NouraguesHD dataset to the field_inv_dataset variable
field_inv_dataset <- NouraguesHD

# Correct the taxonomy of the species in the field_inv_dataset dataset
Taxo <- correctTaxo(genus = field_inv_dataset$genus, species = field_inv_dataset$species)
# field_inv_dataset$genusCorr <- Taxo$genusCorrected
# field_inv_dataset$speciesCorr <- Taxo$speciesCorrected

field_inv_dataset$genusCorr <- field_inv_dataset$genus
field_inv_dataset$speciesCorr <- field_inv_dataset$species


# Retrieve taxonomic information (family and order) based on the corrected genus
APG <- getTaxonomy(field_inv_dataset$genusCorr, findOrder = TRUE)
field_inv_dataset$familyAPG <- APG$family
field_inv_dataset$orderAPG <- APG$order

# Retrieve wood density information based on the corrected genus, species, and plot information
dataWD <- getWoodDensity(genus = field_inv_dataset$genusCorr, species = field_inv_dataset$speciesCorr, stand = field_inv_dataset$plotId)

# Compute AGB using the modelHD function with the provided dataset
result <- modelHD(
  D = field_inv_dataset$D,
  H = field_inv_dataset$H,
  useWeight = TRUE
)

# Print the result of the modelHD computation
print(result) #select the HD model with lowest Bias (in this case "log2" model)

# Compute the HDmodel using the modelHD function with the provided dataset
HDmodel <- modelHD(
  D = field_inv_dataset$D,
  H = field_inv_dataset$H,
  method = "log2",
  useWeight = TRUE
)

# Retrieve H values using the retrieveH function based on the HDmodel
dataHlocal <- retrieveH(
  D = field_inv_dataset$D,
  model = HDmodel
)

# Assign wood density values and standard deviation to the field_inv_dataset
field_inv_dataset$WD <- dataWD$meanWD
field_inv_dataset$sdWD <- dataWD$sdWD

# Since few trees in the dataset have directly measured height on ground, we assume 
# error in height to be 0.5m for these trees and use the modeled height for 
# trees with missing height measurements. 

# Prepare the mixed H values and RSE values in the field_inv_dataset
field_inv_dataset$Hmix <- field_inv_dataset$H
field_inv_dataset$RSEmix <- 0.5
filt <- is.na(field_inv_dataset$Hmix)
field_inv_dataset$Hmix[filt] <- dataHlocal$H[filt]
field_inv_dataset$RSEmix[filt] <- HDmodel$RSE

# Perform AGB monte carlo simulation using the AGBmonteCarlo function with the field_inv_dataset
resultMC <- AGBmonteCarlo(
  D = field_inv_dataset$D,
  WD = field_inv_dataset$WD,
  errWD = field_inv_dataset$sdWD,
  H = field_inv_dataset$Hmix,
  errH = field_inv_dataset$RSEmix,
  Dpropag = "chave2004"
)

# Create a data frame to store the simulated AGB values and assign plot IDs
AGB_simulated <- data.frame(resultMC$AGB_simu)
AGB_simulated$PlotID <- field_inv_dataset$plotId

# Group the AGB_simulated data frame by PlotID and compute the sum of AGB values
AGB1000_plotlevel <- AGB_simulated %>%
  group_by(PlotID) %>%
  summarise_all(sum)

# Rename the column names of AGB1000_plotlevel to AGB1, AGB2, ..., AGB1000
# each column indicating the total AGB per one MonteCarlo simulation
names(AGB1000_plotlevel)[2:1001] <- paste0("AGB",1:1000)


# Part-2 (LiDAR-Metrics Extraction) ----

# CHM and Plot-shapefiles are not shared in this study
# The code is provided to extract LiDAR metrics using the shapefile of the plots.

chm <- raster("path_to_the_chmfile") #in GeoTiff Format
plots100m <- shapefile("path_to_the_plots_shapefile") #can be 40m/100m plots

#merge plot-shapefile and AGB simulations
plots100m <- merge(plots100m,AGB1000_plotlevel,by="PlotID")

#Function to extract CHM metrics for the study.
StructureEstimates <- function(plots100m_shp) {
  
  for (i in 1:length(plots100m_shp$PlotID)) {
    print(i)
    tmp <- plots100m_shp[i,]
    
    # CHM Led Structure -
    tmp_vals <- unlist(extract(chm,tmp))
    tmp_vals <- tmp_vals[!is.na(tmp_vals)]
    
    plots100m_shp$Count[i] <- length(tmp_vals)
    
    plots100m_shp$RH40[i] <- quantile(tmp_vals,0.4)
    plots100m_shp$RH50[i] <- quantile(tmp_vals,0.5)
    plots100m_shp$RH60[i] <- quantile(tmp_vals,0.6)
    plots100m_shp$RH70[i] <- quantile(tmp_vals,0.7)
    plots100m_shp$RH80[i] <- quantile(tmp_vals,0.8)
    plots100m_shp$RH90[i] <- quantile(tmp_vals,0.9)
    plots100m_shp$RH98[i] <- quantile(tmp_vals,0.98)
    
    plots100m_shp$meanH[i] <- mean(tmp_vals)
    plots100m_shp$sdH[i] <- sd(tmp_vals) # rugosity.
    plots100m_shp$CV[i] <- cv(tmp_vals)
    
    plots100m_shp$CCF2[i] <- sum(tmp_vals>2)/length(tmp_vals)
    plots100m_shp$CCF5[i] <- sum(tmp_vals>5)/length(tmp_vals)
    plots100m_shp$CCF10[i] <- sum(tmp_vals>10)/length(tmp_vals)
    
    plots100m_shp$QMCH_chm[i] <- QMCH(tmp_vals)
    
    tryCatch({
      chm_clip <- crop(chm,tmp)
      plots100m_shp$rumple[i] <- rumple_index(chm_clip)
      
    }, error=function(e){})
  }
  return(plots100m_shp)
}

#Extract the Structure Metrics using the custom function
plots100m <- StructureEstimates(plots100m)

# Part-3 (LiDAR-AGB Modelling) ----
tmpdata <- plots100m@data
AGBcolumns <- paste0("AGB",1:1000)
tmpdata$AGB <- mean(tmpdata[,AGBcolumns])

tmpdata <- tmpdata[,c("AGB","meanH")]
tmpdata <- log(tmpdata)

mod_log <-  lm(AGB~meanH,data = tmpdata) 
summary(mod_log)
rse <- summary(mod_log)$sigma
rse
predAGB <- predict(mod_log)
predAGB_unloggued = exp(predAGB) * exp(rse^2/2) # correction of back-transformed predictions

y <- exp(tmpdata$AGB) # reference field-derived AGB in natural unit
resid = y - predAGB_unloggued 
RSE=sqrt(sum(resid^2)/(length(resid)-2))
R2 = cor(predAGB, tmpdata$AGB)^2 # model R2, in log unit
RSE
R2
RMSE <- sqrt( mean(resid^2))
RMSE
round(mod_log$coefficients,3) #modelled coefficients
RMSE/mean(y) # Relative RMSE

# Part-4 (Error Propagation and Output Map Generation) ----
source("metropolis_algo.R") #sourced from the original paper - circulated in the current respository
#see Appendix S1 of Réjou-Méchain et al.18 for codes and details
# Réjou-Méchain, M., Tanguy, A., Piponiot, C., Chave, J. & Hérault, B. 
#biomass: an r package for estimating above-ground biomass and 
#its uncertainty in tropical forests. Methods Ecol. Evol. 8, 1163-1167 (2017).

tmpdata <- plots100m@data
pred=c()
res = 100 #for 100m maps using 100m plot model
Rast_Hm = aggregate(chm, fact = res, fun = nofun, expand = F)
crs(Rast_Hm) <- crs(chm)

newdataMap=as.data.frame(Rast_Hm)
names(newdataMap)=c("Hm")
RSEsim_out  <- NA

for (i in 1:1000){
  print(paste0(i,"/",1000))
  # Build the model with simulated AGB: Error field
  AGBcolid <- paste0("AGB",i)
  datasim=data.frame(AGB=tmpdata[,AGBcolid],Hm=tmpdata$meanH)
  modsim=lm(log(AGB)~log(Hm),data=datasim)

  param4_site <- param4(datasim$AGB,datasim$Hm) 
  #generating the coefficients of LiDAR-AGB model using Markov chain Monte Carlo algorithm
  
  sel_id <- sample(1:nrow(param4_site),1) #only one sample is selected.
  Ea <- param4_site[sel_id, "intercept"] #selecting the respective intercept of the sample
  Eb <- param4_site[sel_id, "logagbt"] #selecting the respective slope of the sample
  RSE <- param4_site[sel_id, "sd"] #selecting the respective RSE of the sample
  
  RSEsim <- rnorm(sd = RSE, n = 1) #generating a random error within the RSE
  
  # Predict values for maps by propagating a random error from the LiDAR-AGB model.
  predsim <- (log(newdataMap$Hm) * Eb) + Ea
  predsim <- predsim + RSEsim
  predsim <- exp(predsim)
  
  RSEsim_out[i] <- RSEsim
  pred=cbind(pred,predsim)
}

AGBmean=(apply(pred,1,mean,na.rm=T)) #The mean AGB of 1000 simulated AGB maps
AGBsd=(apply(pred,1,sd,na.rm=T)) # Standard deviation of 1000 simulated AGB Maps
E3_Total=mean((100*AGBsd/AGBmean),na.rm=T)

Rast_Hm$AGBpred <- AGBmean
Rast_Hm$AGBpredsd <- AGBsd
crs(Rast_Hm) <- crs(chm)

AGB_outMap <- Rast_Hm$AGBpred
AGB_outMap$AGBpredsd <- Rast_Hm$AGBpredsd

writeRaster(AGB_outMap,"path_to_output_file",overwrite=T)


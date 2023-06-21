# LIDAR_AGB
The provided script is a part of the study titled "LiDAR-based reference aboveground biomass maps for tropical forests of South Asia and Central Africa," which has been submitted to Nature Scientific Data.

The script focuses on processing tree-level field inventory data to estimate plot-level aboveground biomass (AGB) along with uncertainties using the Monte Carlo Error simulation model. This is achieved by utilizing the BIOMASS R package.

After obtaining the plot-level AGB estimates, the script enables the linking of these estimates to LiDAR Canopy Height Model (CHM) metrics at the desired resolution. Although the study primarily used 100m and 40m resolutions, the code is adaptable to other resolutions as well.

The subsequent step involves modelling the plot-level AGB using a chosen LiDAR CHM metric, specifically the "meanH" metric in this case.

To generate AGB maps, the script utilizes Monte Carlo simulation to propagate the plot-level errors and the LiDAR-AGB model error. The script has been set up to perform 1000 simulations for this purpose.

Finally, the computed 1000 AGB maps can be merged to obtain site-level meanAGB and sdAGB (standard deviation of AGB) values. This consolidation allows for a comprehensive understanding of the AGB distribution at the site level.
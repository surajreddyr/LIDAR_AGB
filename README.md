# LIDAR_AGB
The provided script is a part of the study titled "LiDAR-based reference aboveground biomass maps for tropical forests of South Asia and Central Africa," which has been submitted to Nature Scientific Data.

The script focuses on processing tree-level field inventory data to estimate plot-level aboveground biomass (AGB) along with uncertainties using the Markov chain Monte Carlo Error simulation model. This is achieved by utilizing the BIOMASS R package.

After obtaining the plot-level AGB estimates, the script enables the linking of these estimates to LiDAR Canopy Height Model (CHM) metrics at the desired resolution. Although the study primarily used 100m and 40m resolutions, the code is adaptable to other resolutions as well.

The subsequent step involves modelling the plot-level AGB using a chosen LiDAR CHM metric, specifically the "meanH" metric in this case.

To generate AGB maps, the script mimics the propagation algorithm of Biomass package to propagate the plot-level errors and the LiDAR-AGB model errors to generate LiDAR-AGB maps. The script has been set up to perform 1000 simulations for this purpose.

Finally, the computed 1000 AGB maps can be merged to obtain site-level meanAGB and sdAGB (standard deviation of AGB) values. This consolidation allows for a comprehensive understanding of the AGB distribution at the site level.

The files "McRoberts2022_implementation_04Dec23.R" and "Saarela2020_3pHMB_3pHHY_MC-based_valid_inputs.R" are the implementation of respective papers!

To be updated more after the acceptance of the research paper.

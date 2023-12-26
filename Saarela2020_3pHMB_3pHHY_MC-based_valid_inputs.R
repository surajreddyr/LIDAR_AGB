rm(list=ls())
setwd("K:/AMAP_2022/LIDAR_AGB/Paper_2022/Hierarchical_error_model")

fulldata <- readRDS("final_data_100m_saarela_04Dec23.rds")
fulldata$DBH <- sqrt((fulldata$BA/fulldata$Density)/pi)

siteids <- unique(fulldata$siteID)

for (i in 1:length(siteids)) {
  tmpdata <- fulldata[fulldata$siteID==siteids[i],]
  # X -> DBH | Y -> AGB | Z -> meanH
  mux = mean(tmpdata$DBH*100)
  sdx = sd(tmpdata$DBH*100)
  
  muy = mean(tmpdata$AGB)
  sdy = sd(tmpdata$AGB)
  
  muz = mean(tmpdata$meanH)
  sdz = sd(tmpdata$meanH)
  
  rxy = cor(tmpdata$DBH,tmpdata$AGB)
  rxz = cor(tmpdata$DBH,tmpdata$meanH)
  ryz = cor(tmpdata$AGB,tmpdata$meanH)
  
  n <- length(tmpdata$PlotID)
  
  tmp_out <- c(siteids[i],n,mux,sdx,muy,sdy,muz,sdz,rxy,rxz,ryz)
  
  
  if (i==1){
    params_out <- tmp_out
  } else {
    params_out <- rbind(params_out,tmp_out)
  }
  
  print(siteids[i])
  print(sum(tmpdata$Density))
}

params_out <- data.frame(params_out)
names(params_out) <- c("SiteID","N","Mu_X","SD_X","Mu_Y","SD_Y","Mu_Z","SD_Z",
                       "R_XY","R_XZ","R_YZ")

write.csv(params_out,"tbl_out_100m_plots_saarela.csv",row.names=F)

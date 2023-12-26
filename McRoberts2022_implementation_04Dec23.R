# ALS - AGB Relationships Script
rm(list=ls())
setwd("K://AMAP_2022//LIDAR_AGB/Paper_2022/Hierarchical_error_model/McRoberts2022/")
library(raster)
library(svMisc)

# Input Data ----
res = 40
tmpdata <- readRDS(paste0("../final_data_",res,"m_saarela_04Dec23.rds")) #field samples - AGB and meanH for the site
tmpdata <- tmpdata[tmpdata$siteID=="ACK",]
#tmpdata <- tmpdata[tmpdata$AGB<1100,]

chm <- raster("CHMs_all/Achanakmar_CHM_Apr2021_1m_extend.tif")

nofun <- function(x,...){if(all(is.na(x))) NA else mean(x,na.rm=TRUE)}
Rast_Hm = aggregate(chm, fact = res, fun = nofun, expand = F)

names(Rast_Hm) <- "CHMval"
Rast_Hm=as.data.frame(Rast_Hm)

# power model function
backlogmod=function(y,x){
  m=lm(log(y)~log(x))
  #print(summary(m))
  rse=summary(m)$sigma
  coeffa=exp(m$coefficients[1]+0.5*rse^2) #applying baskerville correction
  coeffb=m$coefficients[2]
  yi.hat <- coeffa*x^coeffb
  residback=y- yi.hat #equation 7; epsilon.hat; yi.hat is predicted values of y from x using the power model.
  rseback=sqrt(sum(residback^2)/(length(residback)-2))
  AICback = extractAIC(m)+sum(2*log(y))
  resu=list(coeffa=coeffa,coeffb=coeffb,residback=residback,rseback=rseback,AIC=AICback,yi.hat=yi.hat)
}

# power model fit.
modlogAGB=backlogmod(tmpdata$AGB,tmpdata$meanH)
beta0 <- modlogAGB$coeffa #Equation 6a
beta1 <- modlogAGB$coeffb #Equation 6a
episilon.hat <- modlogAGB$residback #equation 7 (epsilon hat)

mu.mean.b <- mean(modlogAGB$yi.hat) #actual mean of the predicted AGB values.

#print(mu.mean.b) 
plot(tmpdata$meanH,tmpdata$AGB)

## Residual Variability (Section 3.2) ----
#Step-1: For sample units i = 1,2,.. n, sort the pairs (xi, ei.hat) from least to greatest with respect to xi
data1 <- data.frame(xi = tmpdata$meanH,ei.hat = episilon.hat,yi = tmpdata$AGB,yi.hat=modlogAGB$yi.hat)
data1 <- data1[order(data1$xi),]
data1$ID <- 1:nrow(data1)

#Step-2: Aggregate the sorted pairs into groups, g, of size ng >= 10, depending on the overall sample size
ng = 10 #initially for testing
data1$groups <- as.integer(data1$ID/3) + 1
#data1$groups[1] <- 1

#step-3 : For each group, g, calculate the mean within-group value of the indepdent variable xg.bar and within-group variance sigma.g.sq
library(dplyr)
data1_summary <- data1 %>% group_by(groups) %>%
  summarise(xg.bar = mean(xi),
            sigma.g.sq = var(ei.hat))

#step-4: fir a model h, often a power or exponential model, sigma.g.sq = h(xg.bar; alpha) + delta
moddata1 <- backlogmod(data1_summary$sigma.g.sq,data1_summary$xg.bar)

#step-5 : Estimate the residual variance for any particular sample or population unit as sigma.i.sq = h(xi,alpha.hat) 
alpha0 <- moddata1$coeffa
alpha1 <- moddata1$coeffb
sigma.i.sq <- alpha0 * ((data1$xi)^alpha1) 


# Sampling Variability : Bootstrap Variance Estimator (Section 3.3.1) ----

#Step-1: For bootstrap replication,b, select a random sample with replacement from the orignal sample and the same size as the orginal sample
b_sample <- sample(1:nrow(tmpdata),nrow(tmpdata),replace=T)
b_sample <- tmpdata[b_sample,]

#Step-2: Estimate or optimize the parameters associated with the prediction technique using the resample
mod_sample <- backlogmod(b_sample$AGB,b_sample$meanH)

#Step-3: Using the predicted values from step-2, calculate model-based prediction of the population mean
mu.mean.b <- mean(mod_sample$yi.hat)
refpred <- mu.mean.b
#Step-4: Repeat steps 1-3 nboot times.

# Fixing of nboot? through iteration based on stopping criteria.
#nboot = 5000

N=as.data.frame(Rast_Hm) #N is the total population of the CHM map
names(N)=c("N")
N <- N[!is.na(N)]
N <- N[N>0]

mu.mean.b <- NA
var.pop.iterations <- NA

anyflag = 9999
nboot_i <- 1
while (anyflag !=0) {
  
  b_sample <- sample_n(tmpdata,nrow(tmpdata),replace=T)
  mod_sample <- backlogmod(b_sample$AGB,b_sample$meanH)
  yi.b.hat.temp <- mod_sample$coeffa * (N ^ mod_sample$coeffb)
  mu.mean.b[nboot_i] <- mean(yi.b.hat.temp)
  mu.mean.b[mu.mean.b=="Inf"] <- refpred
  var.pop.iterations[nboot_i] <- var(mu.mean.b)
  
  if (nboot_i > 1000) { #stopping criterion - for which I posed the question.
    #var_final_flag <- 0.01*var(mu.mean.b) #1% change
    tail.50per <- tail(var.pop.iterations[-1],nboot_i/2)
    tail.50per <- tail.50per[-length(tail.50per)]
    tail.deviation <- abs(tail.50per-var(mu.mean.b))
    tail.deviation_per <- 100*(tail.deviation/var(mu.mean.b)) #variability in percent
    
    anyflag <- sum(tail.deviation_per > 2)
    
  }
  
  nboot_i <- nboot_i+1
}

plot(var.pop.iterations,ylim=c(0,1500),type='l')
print(nboot_i)
nboot <- nboot_i

mu.boot = mean(mu.mean.b)
var.pop.sam = var(mu.mean.b)
print(sqrt(var.pop.sam)) #check - 1

# Equivalent Bootstrap variance estimator (Section 3.3.2) ---- 
N=as.data.frame(Rast_Hm) #N is the total population of the CHM map
names(N)=c("N")

Nsamp <- N
Nsamp <- Nsamp[!is.na(Nsamp)]
Nsamp <- Nsamp[Nsamp>0]

#nboot = 5
yi.b.hat <- matrix(data=NA,nrow=length(Nsamp),ncol=nboot) #predictions from N
for (i in 1:nboot){
  b_sample <- sample(1:nrow(tmpdata),nrow(tmpdata),replace=T) #step1 - select a random sample with replacement
  b_sample <- tmpdata[b_sample,] #step1 - select a random sample with replacement
  mod_sample <- backlogmod(b_sample$AGB,b_sample$meanH) #step2 - estimate model parameters
  yi.b.hat.temp <- mod_sample$coeffa * (Nsamp ^ mod_sample$coeffb)
  yi.b.hat[,i] <- (yi.b.hat.temp)
}

cov.matrix <- cov(t(yi.b.hat)) # it is similar as the below code.

# cov.matrix <- matrix(data=NA,nrow=length(Nsamp),ncol=length(Nsamp)) #predictions from N
# 
# # Generate population unit covariance matrices - slowest step?
# for (i in 1:length(Nsamp)){
#   print(i)
#   for (j in 1:length(Nsamp)) {
#     tmp1 <- NA
#     
#     #covariance of yi.b.hat.temp,yj.b.hat.temp
#     cov.matrix[i,j] <- cov(yi.b.hat[i,],yi.b.hat[j,])
#     
#     # covariance original code from paper - direct function from R is used. 
#     # for (b in 1:nboot){
#     #   yi.b.hat.temp <- yi.b.hat[i,b]
#     #   yj.b.hat.temp <- yi.b.hat[j,b]
#     #   yi.hat.bar <- mean(yi.b.hat[i,])
#     #   yj.hat.bar <- mean(yi.b.hat[j,])
#     #   tmp1[b] <- (yi.b.hat.temp-yi.hat.bar)*(yj.b.hat.temp-yj.hat.bar)
#     # }
#     # cov.matrix[i,j] <- sum(tmp1)/(nboot-1)
#   }
# }

var.unit.sam <- sum(cov.matrix)/(length(Nsamp)*length(Nsamp))
print(sqrt(var.unit.sam)) #check - 2

# Binned Equivalent Bootstrap variance estimator (Section 3.3.3) ----
N=as.data.frame(Rast_Hm) #N is the total population of the CHM map
names(N)=c("N")

#Nsamp <- sample(t(N),size = 1000)
Nsamp <- N
Nsamp <- Nsamp[!is.na(Nsamp)]
Nsamp <- Nsamp[Nsamp>0]

binID <- as.integer(Nsamp/3)+1 # generated ~5m height bins - so that we have 10x10 cov matrix

binID_count <- aggregate(binID, list(num=binID), length)
binID_count

binID[binID>9] <- 9
#binID[binID<3] <- 3

#nboot = 1000
yi.b.hat <- matrix(data=NA,nrow=length(Nsamp),ncol=nboot) #predictions from N
for (i in 1:nboot){
  b_sample <- sample(1:nrow(tmpdata),nrow(tmpdata),replace=T) #step1 - select a random sample with replacement
  b_sample <- tmpdata[b_sample,] #step1 - select a random sample with replacement
  mod_sample <- backlogmod(b_sample$AGB,b_sample$meanH) #step2 - estimate model parameters
  yi.b.hat.temp <- mod_sample$coeffa * (Nsamp ^ mod_sample$coeffb)
  yi.b.hat[,i] <- (yi.b.hat.temp)
}

unique_bins <- sort(unique(binID))
no_bins <- length(unique_bins)
Cpq <- matrix(data=NA,nrow=no_bins,ncol=no_bins) #predictions from N
Np <- matrix(data=NA,nrow=no_bins,ncol=no_bins) #predictions from N
Nq <- matrix(data=NA,nrow=no_bins,ncol=no_bins) #predictions from N

# Generate population unit covariance matrices
for (p in 1:no_bins){
  progress(p,no_bins)
  for (q in 1:no_bins) {
    sample.p <- Nsamp[binID==unique_bins[p]]
    sample.q <- Nsamp[binID==unique_bins[q]]
    yi.b.hat.binp <- yi.b.hat[binID==unique_bins[p],]
    yj.b.hat.binq <- yi.b.hat[binID==unique_bins[q],]
    
    if (class(yi.b.hat.binp)=="numeric") {
      yi.b.hat.binp <- as.matrix(t(yi.b.hat.binp))
    }
    if (class(yj.b.hat.binq)=="numeric") {
      yj.b.hat.binq <- as.matrix(t(yj.b.hat.binq))
    }
    
    Np.tmp <- nrow(yi.b.hat.binp)
    Nq.tmp <- nrow(yj.b.hat.binq)
    
    cov.matrix <- cov(t(yi.b.hat.binp),t(yj.b.hat.binq))
    
    # cov.matrix <- matrix(data=NA,nrow=Np.tmp,ncol=Nq.tmp) #predictions from N
    # for (i in 1:Np.tmp){
    #   for (j in 1:Nq.tmp){
    #     cov.matrix[i,j] <- cov(yi.b.hat.binp[i,],yj.b.hat.binq[j,])
    #     # 
    #     # for (b in 1:nboot){
    #     #   yi.b.hat.temp <- yi.b.hat.binp[i,b]
    #     #   yj.b.hat.temp <- yj.b.hat.binq[j,b]
    #     #   yi.hat.bar <- mean(yi.b.hat.binp[i,])
    #     #   yj.hat.bar <- mean(yj.b.hat.binq[j,])
    #     #   tmp1[b] <- (yi.b.hat.temp-yi.hat.bar)*(yj.b.hat.temp-yj.hat.bar)
    #     # }
    #     # cov.matrix[i,j] <- sum(tmp1)/(nboot-1)
    #   }
    # }
    
    Cpq[p,q] <- sum(cov.matrix)
    Np[p,q] <- Np.tmp
    Nq[p,q] <- Nq.tmp
    
    rm(cov.matrix)
    gc()
  }
}

Cpq.final <- Cpq/(Np*Nq)
var.bin.sam <- sum(Cpq.final)/(no_bins^2)
print(sqrt(var.bin.sam))

# Verification ----
# For the population : Var.pop.sam = var.unit.sam = var.bin.sam [Section 3.3.3]
print(sqrt(var.pop.sam))
print(sqrt(var.unit.sam))
print(sqrt(var.bin.sam))

# I observe there are slight differences on the above estimates! Is this okay?
print(sqrt(var.bin.sam)/sqrt(var.unit.sam)) # this should be close!!!
print(sqrt(var.bin.sam)/sqrt(var.pop.sam))

#write.csv(Cpq.final,"2_cov_matrix_ACK_40m.csv",row.names=F)

#******************************************************************
#######**MULTI-season OCCUPANCY/ DYNAMIC ANALYSIS########*******###
#******************************************************************

#####**Packages**####
library(camtrapR)
library(lubridate)
library(unmarked)
library(AICcmodavg)
library(stats)  


#####Cam-op-MAtrix####
##**Creation of camera operation matrix**##
SA_TE <- read.csv("SA_TE1.csv")
summary(SA_TE)
class(SA_TE$Setup_date)
#Changing format of columns to date#
SA_TE$Setup_date <- as.Date(SA_TE$Setup_date)
SA_TE$Retrieval_date <- as.Date(SA_TE$Retrieval_date)
SA_TE$Problem1_from <- as.Date(SA_TE$Problem1_from)
SA_TE$Problem1_to <- as.Date(SA_TE$Problem1_to)

##**yearly and site covariates**####
# Create yearlySitecovs
years1 <- matrix(unique(as.character(SA_TE$Session)),
                 ncol = length(unique(as.character(SA_TE$Session))),
                 nrow = length(unique(SA_TE$Station)),
                 byrow = TRUE )

colnames(years1)=c('year1', 'year2','year3')

## create Site covs for multiseason
site_covariates <- read.csv("Sitecov_St.csv", header = TRUE, row.names="Station")

#####DetectionHist-MAtrix####
###**Panthera onca***###
DetHist_Jag <- read.csv("Panthera_10d.csv",header = TRUE, row.names="X")
effort <- read.csv("Panthera_10dEf.csv",header = TRUE, row.names="X")

DetHistJagx10 <- list(detection_history = DetHist_Jag, effort = effort)


#####Create unmarkedMult frame####
UMFJag <- unmarked::unmarkedMultFrame(y = DetHistJagx10$detection_history,
                                      siteCovs = site_covariates[,1:14],
                                      yearlySiteCovs = list(year = years1, 
                                                            NDVI = site_covariates[,c(3:5)],
                                                            NBR = site_covariates[,c(6:8)],
                                                            AB = site_covariates[,c(9:11)],
                                                            BP = site_covariates[,c(12:14)]),
                                      obsCovs = list(effort = DetHistJagx10$effort), 
                                      numPrimary = 3)


####RUNNING MODELS - Jaguars####
##NuLL Model
JaNull <- unmarked::colext(psiformula = ~1,         
                           gammaformula = ~1,          
                           epsilonformula = ~1,      
                           pformula = ~1,         
                           data = UMFJag ,
                           method="BFGS")

##Running detection#
##**DETECTION (p)**
Jap1 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~Trail ,         
                         data = UMFJag ,
                         method="BFGS")
Jap2 <- unmarked::colext(psiformula = ~ 1,            
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~Trail + effort ,         
                         data = UMFJag ,
                         method="BFGS")
Jap3 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~Trail + year ,
                         data = UMFJag ,
                         method="BFGS")

Jap4 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~effort ,
                         data = UMFJag ,
                         method="BFGS")

Jap5 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~year ,
                         data = UMFJag ,
                         method="BFGS")

Jap6 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~effort + year ,
                         data = UMFJag ,
                         method="BFGS")

Jap7 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~Trail + effort + year ,
                         data = UMFJag ,
                         method="BFGS")

modlistJa_p <-list(Jap.null=JaNull ,Jap.Trail=Jap1, Jap.Trail.Effort=Jap2,Jap.Trail.year=Jap3,
                   Jap.Effort=Jap4,Jap.year=Jap5,Jap.Effort.year=Jap6,Jap.Trail.Effort.year=Jap7)
aictab(modlistJa_p)

#**OCCUPANCY (psi)**#
Japsi1 <- unmarked::colext(psiformula = ~Dwater,           
                           gammaformula = ~ 1,          
                           epsilonformula = ~ 1,      
                           pformula = ~effort ,         
                           data = UMFJag ,
                           method="BFGS")

Japsi2 <- unmarked::colext(psiformula = ~NDVI1,            
                           gammaformula = ~ 1,          
                           epsilonformula = ~ 1,      
                           pformula = ~effort ,         
                           data = UMFJag ,
                           method="BFGS")

Japsi3 <- unmarked::colext(psiformula = ~Dwater+NDVI1,            
                           gammaformula = ~ 1,          
                           epsilonformula = ~ 1,      
                           pformula = ~effort ,         
                           data = UMFJag ,
                           method="BFGS")

modlistJa_psi <-list(JaNull= JaNull,Japsi1.p.Eff=Jap4,Japsi.Dwater.p.Eff=Japsi1,
                     Japsi.NDVI.p.Eff=Japsi2,
                     Japsi.Dwater.NDVI.p.Eff=Japsi3)
aictab(modlistJa_psi)

##**Now Gamma (g)**##
Jag1 <- unmarked::colext(psiformula = ~1,            
                         gammaformula = ~AB,            
                         epsilonformula = ~ 1,       
                         pformula = ~effort ,         
                         data = UMFJag ,
                         method="BFGS")


ModlistJa_g1 <- list(JaNull= JaNull,Japsi.1.p.E.g1=Jap4,Japsi1.p.E.g.AB=Jag1)
aictab(ModlistJa_g1)

#**Now Epsilon (e)**#

Jae1 <- unmarked::colext(psiformula = ~1,            
                         gammaformula = ~1,             
                         epsilonformula = ~AB,  
                         pformula = ~effort ,         
                         data = UMFJag ,
                         method="BFGS")


ModlistJa_e1 <- list(JaNull= JaNull,Japsi1.p.E.g1.e1=Jap4,
                     Japsi1.p.E.g1.e.AB=Jae1)
aictab(ModlistJa_e1)

#GOF
mb.bootJ <- AICcmodavg::mb.gof.test(Jae1, nsim = 1000)
print(mb.bootJ, digit.vals = 4, digits.chisq = 4)


# Predicted occupancy in each year
mJ <- nonparboot(Jae1,  B = 1000)
predicted_occupancyJ <- data.frame(year = c(1:3),
                                   smoothed_occ = smoothed(Jae1)[2,], 
                                   SE = mJ@smoothed.mean.bsse[2,])
predicted_occupancyJ$year<-as.character(predicted_occupancyJ$year)

# Calculate the critical value for a desired confidence level
confidence_level <- 0.95
critical_value <- qnorm((1 + confidence_level) / 2)  # Assuming normal distribution and 95% confidence level
# Calculate the confidence interval for each row
predicted_occupancyJ$lcl <- predicted_occupancyJ$smoothed_occ - (critical_value * predicted_occupancyJ$SE)
predicted_occupancyJ$ucl <- predicted_occupancyJ$smoothed_occ + (critical_value * predicted_occupancyJ$SE)

predicted_occupancyJ

#*#*#*#**#*#*#*#*#*#*#**#*#*#*#*#*#*#*#**#*#*#*#*#*#*#**#*#*#
###**Puma concolor***###
DetHist_Pum <- read.csv("Puma_15d.csv",header = TRUE, row.names="X")
effortP <- read.csv("Puma_15dEf.csv",header = TRUE, row.names="X")

DetHistPumx15 <- list(detection_history = DetHist_Pum, effort = effortP)

UMFPu <- unmarked::unmarkedMultFrame(y = DetHistPumx15$detection_history,
                                     siteCovs = site_covariates[,1:14],
                                     yearlySiteCovs = list(year = years1, 
                                                           NDVI = site_covariates[,c(3:5)],
                                                           NBR = site_covariates[,c(6:8)],
                                                           AB = site_covariates[,c(9:11)],
                                                           BP = site_covariates[,c(12:14)]),
                                     obsCovs = list(effort = DetHistPumx15$effort), 
                                     numPrimary = 3)

##NuLL Model
PuNull <- unmarked::colext(psiformula = ~1,           
                           gammaformula = ~1,          
                           epsilonformula = ~1,      
                           pformula = ~1,         
                           data = UMFPu ,
                           method="BFGS")
##DETECTION (p)
Pup1 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~Trail ,         
                         data = UMFPu ,
                         method="BFGS")
Pup2 <- unmarked::colext(psiformula = ~ 1,            
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~Trail + effort ,         
                         data = UMFPu ,
                         method="BFGS")
Pup3 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~Trail + year ,
                         data = UMFPu ,
                         method="BFGS")

Pup4 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~effort ,
                         data = UMFPu ,
                         method="BFGS")

Pup5 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~year ,
                         data = UMFPu ,
                         method="BFGS")

Pup6 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~effort + year ,
                         data = UMFPu ,
                         method="BFGS")

Pup7 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~Trail + effort + year ,
                         data = UMFPu ,
                         method="BFGS")

modlistPu_p <-list(Pup.null=PuNull ,Pup.Trail=Pup1, Pup.Trail.Effort=Pup2,Pup.Trail.year=Pup3,
                   Pup.Effort=Pup4,Pup.year=Pup5,Pup.Effort.year=Pup6,Pup.Trail.Effort.year=Pup7)
aictab(modlistPu_p)

##GOF test to choose best covariate
mb.bootp1 <- AICcmodavg::mb.gof.test(Pup4, nsim = 1000)
print(mb.bootp1, digit.vals = 4, digits.chisq = 4)
mb.bootp2 <- AICcmodavg::mb.gof.test(Pup2, nsim = 1000)
print(mb.bootp2, digit.vals = 4, digits.chisq = 4)

#**OCCUPANCY (psi)**#
Pupsi1 <- unmarked::colext(psiformula = ~Dwater,          
                           gammaformula = ~ 1,          
                           epsilonformula = ~ 1,      
                           pformula = ~effort+Trail ,         
                           data = UMFPu ,
                           method="BFGS")

Pupsi2 <- unmarked::colext(psiformula = ~NDVI1,          
                           gammaformula = ~ 1,          
                           epsilonformula = ~ 1,      
                           pformula = ~effort+Trail ,         
                           data = UMFPu ,
                           method="BFGS")

Pupsi3 <- unmarked::colext(psiformula = ~Dwater+NDVI1,            
                           gammaformula = ~ 1,          
                           epsilonformula = ~ 1,      
                           pformula = ~effort+Trail ,         
                           data = UMFPu ,
                           method="BFGS")

modlistPu_psi1 <-list(PuNull= PuNull,Pupsi1.p.E.T=Pup2,Pupsi.Dwater.p.E.T=Pupsi1,
                      Pupsi.NDVI.p.E.T=Pupsi2,
                      Pupsi.Dwater.NDVI.p.E.T=Pupsi3)

aictab(modlistPu_psi1)

#**Gamma (g)**#
Pug1 <- unmarked::colext(psiformula = ~NDVI1,            
                         gammaformula = ~AB,           
                         epsilonformula = ~ 1,       
                         pformula = ~effort+Trail ,         
                         data = UMFPu ,
                         method="BFGS")

ModlistPu_g1 <- list(PuNull= PuNull, Pupsi.NDVI.p.E.T.g1=Pupsi2,
                     Pupsi.NDVI.p.E.T.g.AB=Pug1)
aictab(ModlistPu_g1)

#**Epsilon (e)**#
Pue1 <- unmarked::colext(psiformula = ~NDVI1,            
                         gammaformula = ~1,             
                         epsilonformula = ~AB,       
                         pformula = ~effort+Trail ,         
                         data = UMFPu ,
                         method="BFGS")
ModlistPu_e1 <- list(PuNull= PuNull, Pupsi.NDVI.p.E.T.g1.e1=Pupsi2,
                     Pupsi.NDVI.p.E.T.g1.e.AB=Pue1)
aictab(ModlistPu_e1)

#GOF
mb.bootP <- AICcmodavg::mb.gof.test(Pupsi2, nsim = 1000)
print(mb.bootP, digit.vals = 4, digits.chisq = 4)


# Predicted occupancy in each year
mP <- nonparboot(Pupsi2,  B = 1000)
predicted_occupancyP <- data.frame(year = c(1:3),
                                   smoothed_occ = smoothed(Pupsi2)[2,], 
                                   SE = mP@smoothed.mean.bsse[2,])
predicted_occupancyP$year<-as.character(predicted_occupancyP$year)

# Calculate Confidence levels
predicted_occupancyP$lcl <- predicted_occupancyP$smoothed_occ - (critical_value * predicted_occupancyP$SE)
predicted_occupancyP$ucl <- predicted_occupancyP$smoothed_occ + (critical_value * predicted_occupancyP$SE)

predicted_occupancyP

#*#*#*#**#*#*#*#*#*#*#**#*#*#*#*#*#*#*#**#*#*#*#*#*#*#**#*#*#
###**Leopardus pardalis***###
DetHist_Oc <- read.csv("Leopardus_21d.csv",header = TRUE, row.names="X")
effortO <- read.csv("Leopardus_21dEf.csv",header = TRUE, row.names="X")

DetHistOcex21 <- list(detection_history = DetHist_Oc, effort = effortO)

UMFOce <- unmarked::unmarkedMultFrame(y = DetHistOcex21$detection_history,
                                      siteCovs = site_covariates[,1:14],
                                      yearlySiteCovs = list(year = years1, 
                                                            NDVI = site_covariates[,c(3:5)],
                                                            NBR = site_covariates[,c(6:8)],
                                                            AB = site_covariates[,c(9:11)],
                                                            BP = site_covariates[,c(12:14)]),
                                      obsCovs = list(effort = DetHistOcex21$effort), 
                                      numPrimary = 3)

##NuLL Model
OcNull <- unmarked::colext(psiformula = ~1,           
                           gammaformula = ~1,          
                           epsilonformula = ~1,      
                           pformula = ~1,         
                           data = UMFOce ,
                           method="BFGS")
##DETECTION (p)
Ocp1 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~Trail ,         
                         data = UMFOce ,
                         method="BFGS")
Ocp2 <- unmarked::colext(psiformula = ~ 1,        
                         gammaformula = ~ 1,        
                         epsilonformula = ~ 1,       
                         pformula = ~Trail + effort ,         
                         data = UMFOce ,
                         method="BFGS")
Ocp3 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~Trail + year ,
                         data = UMFOce ,
                         method="BFGS")

Ocp4 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~effort ,
                         data = UMFOce ,
                         method="BFGS")

Ocp5 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~year ,
                         data = UMFOce ,
                         method="BFGS")

Ocp6 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~effort + year ,
                         data = UMFOce ,
                         method="BFGS")

Ocp7 <- unmarked::colext(psiformula = ~ 1,           
                         gammaformula = ~ 1,          
                         epsilonformula = ~ 1,        
                         pformula = ~Trail + effort + year ,
                         data = UMFOce ,
                         method="BFGS")

modlistOc_p <-list(Ocp.null=OcNull ,Ocpsi1.p.T.g1.e1=Ocp1, Ocp.Trail.Effort=Ocp2,Ocp.Trail.year=Ocp3,
                   Ocp.E=Ocp4,Ocp.year=Ocp5,Ocp.E.y=Ocp6,Ocp.Trail.Effort.year=Ocp7)
aictab(modlistOc_p)

#**OCCUPANCY (psi)**#
Ocpsi1 <- unmarked::colext(psiformula = ~Dwater,           
                           gammaformula = ~ 1,          
                           epsilonformula = ~ 1,      
                           pformula = ~effort+year ,         
                           data = UMFOce ,
                           method="BFGS")

Ocpsi2 <- unmarked::colext(psiformula = ~NDVI1,          
                           gammaformula = ~ 1,          
                           epsilonformula = ~ 1,      
                           pformula = ~effort+year ,         
                           data = UMFOce ,
                           method="BFGS")

Ocpsi3 <- unmarked::colext(psiformula = ~Dwater+NDVI1,            
                           gammaformula = ~ 1,          
                           epsilonformula = ~ 1,      
                           pformula = ~effort+year ,         
                           data = UMFOce ,
                           method="BFGS")
modlistOc_psi <-list(OcNull= OcNull,Ocpsi1.p.E.y=Ocp6,Ocpsi.Dwater.p.E.y=Ocpsi1,
                     Ocpsi.NDVI.p.E.y=Ocpsi2,
                     Ocpsi.Dwater.NDVI.p.E.y=Ocpsi3)

aictab(modlistOc_psi)

#**Gamma (g)**#
Ocg1 <- unmarked::colext(psiformula = ~1,            
                         gammaformula = ~AB,            
                         epsilonformula = ~ 1,       
                         pformula = ~effort+year ,         
                         data = UMFOce ,
                         method="BFGS")
ModlistOc_g <- list(OcNull= OcNull, Ocpsi1.P.E.y.g1=Ocp6,Ocpsi1.P.E.y.g.AB=Ocg1)
aictab(ModlistOc_g)

##Top model did not converge so we chose the next best one

#**Epsilon (e)**#
Oce1 <- unmarked::colext(psiformula = ~1,            
                         gammaformula = ~1,             
                         epsilonformula = ~AB,     
                         pformula = ~effort+year ,         
                         data = UMFOce ,
                         method="BFGS")
ModlistOc_e <- list(OcNull= OcNull,Ocpsi1.p.E.y=Ocp6,
                    Ocpsi1.p.E.y.gAB.e.AB=Oce1)
aictab(ModlistOc_e)

#GOF
mb.bootO <- AICcmodavg::mb.gof.test(Ocp6, nsim = 1000)
print(mb.bootO, digit.vals = 4, digits.chisq = 4)


# Predicted occupancy in each year
mO <- nonparboot(Ocp6,  B = 1000)
predicted_occupancyO <- data.frame(year = c(1:3),
                                   smoothed_occ = smoothed(Ocp6)[2,], 
                                   SE = mO@smoothed.mean.bsse[2,])
predicted_occupancyO$year<-as.character(predicted_occupancyO$year)

# Calculate Confidence levels
predicted_occupancyO$lcl <- predicted_occupancyO$smoothed_occ - (critical_value * predicted_occupancyO$SE)
predicted_occupancyO$ucl <- predicted_occupancyO$smoothed_occ + (critical_value * predicted_occupancyO$SE)

predicted_occupancyO

#*#*#*#**#*#*#*#*#*#*#**#*#*#*#*#*#*#*#**#*#*#*#*#*#*#**#*#*#
###**Tapirus terrestris***###
DetHist_Tap <- read.csv("Tapirus_14d.csv",header = TRUE, row.names="X")
effortT <- read.csv("Tapirus_14dEf.csv",header = TRUE, row.names="X")

DetHistTapx14 <- list(detection_history = DetHist_Tap, effort = effortT)

UMFTap <- unmarked::unmarkedMultFrame(y = DetHistTapx14$detection_history,
                                      siteCovs = site_covariates[,1:14],
                                      yearlySiteCovs = list(year = years1, 
                                                            NDVI = site_covariates[,c(3:5)],
                                                            NBR = site_covariates[,c(6:8)],
                                                            AB = site_covariates[,c(9:11)],
                                                            BP = site_covariates[,c(12:14)]),
                                      obsCovs = list(effort = DetHistTapx14$effort), 
                                      numPrimary = 3)

##NuLL Model
TapNull <- unmarked::colext(psiformula = ~1,           
                            gammaformula = ~1,          
                            epsilonformula = ~1,      
                            pformula = ~1,         
                            data = UMFTap ,
                            method="BFGS")

##DETECTION (p)
Tapp1 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~Trail ,         
                          data = UMFTap ,
                          method="BFGS")
Tapp2 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~Trail + effort ,         
                          data = UMFTap ,
                          method="BFGS")
Tapp3 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~Trail + year ,
                          data = UMFTap ,
                          method="BFGS")

Tapp4 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~effort ,
                          data = UMFTap ,
                          method="BFGS")

Tapp5 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~year ,
                          data = UMFTap ,
                          method="BFGS")

Tapp6 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~effort + year ,
                          data = UMFTap ,
                          method="BFGS")

Tapp7 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~Trail + effort + year ,
                          data = UMFTap ,
                          method="BFGS")

modlistTap_p <-list(Tapp.null=TapNull ,Tapp.Trail=Tapp1, Tapp.Trail.Effort=Tapp2,Tapp.Trail.year=Tapp3,
                    Tapp.Effort=Tapp4,Tapp.year=Tapp5,Tapp.Effort.year=Tapp6,Tapp.Trail.Effort.year=Tapp7)

aictab(modlistTap_p)

#**OCCUPANCY (psi)**#
Tappsi1 <- unmarked::colext(psiformula = ~Dwater,          
                            gammaformula = ~ 1,          
                            epsilonformula = ~ 1,      
                            pformula = ~Trail + effort + year ,         
                            data = UMFTap ,
                            method="BFGS")

Tappsi2 <- unmarked::colext(psiformula = ~NDVI1,          
                            gammaformula = ~ 1,          
                            epsilonformula = ~ 1,      
                            pformula = ~Trail + effort + year ,         
                            data = UMFTap ,
                            method="BFGS")
Tappsi3 <- unmarked::colext(psiformula = ~Dwater+NDVI1,            
                            gammaformula = ~ 1,          
                            epsilonformula = ~ 1,      
                            pformula = ~Trail + effort + year ,         
                            data = UMFTap ,
                            method="BFGS")
modlistTap_psi <-list(TapNull= TapNull,Tap.psi1.p.T.E.y=Tapp7,Tap.psi.Dw.p.T.E.y=Tappsi1,
                      Tap.psi.NDVI.p.T.E.y=Tappsi2,
                      Tap.psi.Dw.NDVI.p.T.E.y=Tappsi3)

aictab(modlistTap_psi)

#**Gamma (g)**#
Tapg1 <- unmarked::colext(psiformula = ~NDVI1,            
                          gammaformula = ~AB,            
                          epsilonformula = ~ 1,       
                          pformula = ~Trail + effort + year ,         
                          data = UMFTap ,
                          method="BFGS")
ModlistTap_g <- list(TapNull= TapNull,Tappsi.NDVI.p.E.T.Y.g1=Tappsi2,Tapsi.NDVI.p.E.T.Y.g.AB=Tapg1)
aictab(ModlistTap_g)

#**Epsilon (e)**#
Tape1 <- unmarked::colext(psiformula = ~NDVI1,            
                          gammaformula = ~1,             
                          epsilonformula = ~AB,       
                          pformula = ~Trail + effort + year ,         
                          data = UMFTap ,
                          method="BFGS")
ModlistTap_e1 <- list(TapNull= TapNull,Tappsi.NDVI.p.T.E.Y.g1.e1=Tappsi2,
                      Tappsi.NDVI.p.T.E.Y.g1.e.AB=Tape1)
aictab(ModlistTap_e1)

#GOF
mb.bootT <- AICcmodavg::mb.gof.test(Tappsi2, nsim = 1000)
print(mb.bootT, digit.vals = 4, digits.chisq = 4)


# Predicted occupancy in each year
mT <- nonparboot(Tappsi2,  B = 1000)
predicted_occupancyT <- data.frame(year = c(1:3),
                                   smoothed_occ = smoothed(Tappsi2)[2,], 
                                   SE = mT@smoothed.mean.bsse[2,])
predicted_occupancyT$year<-as.character(predicted_occupancyT$year)

# Calculate Confidence levels
predicted_occupancyT$lcl <- predicted_occupancyT$smoothed_occ - (critical_value * predicted_occupancyT$SE)
predicted_occupancyT$ucl <- predicted_occupancyT$smoothed_occ + (critical_value * predicted_occupancyT$SE)

predicted_occupancyT



#*#*#*#**#*#*#*#*#*#*#**#*#*#*#*#*#*#*#**#*#*#*#*#*#*#**#*#*#
###**Mazama americana***###
DetHist_Maz <- read.csv("Mazama_16d.csv",header = TRUE, row.names="X")
effortM <- read.csv("Mazama_16dEf.csv",header = TRUE, row.names="X")

DetHistMazx16 <- list(detection_history = DetHist_Maz, effort = effortM)

UMFMaz <- unmarked::unmarkedMultFrame(y = DetHistMazx16$detection_history,
                                      siteCovs = site_covariates[,1:14],
                                      yearlySiteCovs = list(year = years1, 
                                                            NDVI = site_covariates[,c(3:5)],
                                                            NBR = site_covariates[,c(6:8)],
                                                            AB = site_covariates[,c(9:11)],
                                                            BP = site_covariates[,c(12:14)]),
                                      obsCovs = list(effort = DetHistMazx16$effort), 
                                      numPrimary = 3)

##NuLL Model
MazNull <- unmarked::colext(psiformula = ~1,            
                            gammaformula = ~1,          
                            epsilonformula = ~1,      
                            pformula = ~1,         
                            data = UMFMaz ,
                            method="BFGS")
##Now DETECTION (p)
Mazp1 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~Trail ,         
                          data = UMFMaz ,
                          method="BFGS")
Mazp2 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~Trail + effort ,         
                          data = UMFMaz ,
                          method="BFGS")
Mazp3 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~Trail + year ,
                          data = UMFMaz ,
                          method="BFGS")

Mazp4 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~effort ,
                          data = UMFMaz ,
                          method="BFGS")

Mazp5 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~year ,
                          data = UMFMaz ,
                          method="BFGS")

Mazp6 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~effort + year ,
                          data = UMFMaz ,
                          method="BFGS")

Mazp7 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~Trail + effort + year ,
                          data = UMFMaz ,
                          method="BFGS")

modlistMaz_p <-list(Mazp.null=MazNull ,Mazp.Trail=Mazp1, Mazp.Trail.Effort=Mazp2,Mazp.Trail.year=Mazp3,
                    Mazp.Effort=Mazp4,Mazp.year=Mazp5,Mazp.Effort.year=Mazp6,Mazp.Trail.Effort.year=Mazp7)

aictab(modlistMaz_p)

#**OCCUPANCY (psi)**#
Mazpsi1 <- unmarked::colext(psiformula = ~Dwater,         
                            gammaformula = ~ 1,          
                            epsilonformula = ~ 1,      
                            pformula = ~effort+year ,         
                            data = UMFMaz ,
                            method="BFGS")

Mazpsi2 <- unmarked::colext(psiformula = ~NDVI1,      
                            gammaformula = ~ 1,          
                            epsilonformula = ~ 1,      
                            pformula = ~effort+year ,         
                            data = UMFMaz ,
                            method="BFGS")

Mazpsi3 <- unmarked::colext(psiformula = ~Dwater+NDVI1,            
                            gammaformula = ~ 1,          
                            epsilonformula = ~ 1,      
                            pformula = ~effort+year ,         
                            data = UMFMaz ,
                            method="BFGS")

modlistMaz_psi <-list(MazNull= MazNull,Mazpsi1.p.Eff.y=Mazp6,Mazpsi.Dwater.p.Eff.y=Mazpsi1,
                      Mazpsi.NDVI.p.Eff.y=Mazpsi2,
                      Mazpsi.Dwater.NDVI.p.Eff.y=Mazpsi3)

aictab(modlistMaz_psi)

#**Gamma (g)**#
Mazg1 <- unmarked::colext(psiformula = ~1,            
                          gammaformula = ~AB,            
                          epsilonformula = ~ 1,       
                          pformula = ~effort+year ,         
                          data = UMFMaz ,
                          method="BFGS")


ModlistMaz_g <- list(MazNull= MazNull,Mazpsi.1.p.E.y.g1=Mazp6,Mazpsi1.p.E.y.g.AB=Mazg1)
aictab(ModlistMaz_g)

#**Epsilon (e)**#
Maze1 <- unmarked::colext(psiformula = ~1,            
                          gammaformula = ~1,             
                          epsilonformula = ~AB,    
                          pformula = ~effort+year ,         
                          data = UMFMaz ,
                          method="BFGS")

ModlistMaz_e <- list(MazNull= MazNull,Mazpsi1.p.E.y.g1.e1=Mazp6,Mazpsi1.p.E.y.g1.e.AB=Maze1)
aictab(ModlistMaz_e)

#GOF
mb.bootM <- AICcmodavg::mb.gof.test(Mazp6, nsim = 1000)
print(mb.bootM, digit.vals = 4, digits.chisq = 4)


# Predicted occupancy in each year
mM <- nonparboot(Mazp6,  B = 1000)
predicted_occupancyM <- data.frame(year = c(1:3),
                                   smoothed_occ = smoothed(Mazp6)[2,], 
                                   SE = mM@smoothed.mean.bsse[2,])
predicted_occupancyM$year<-as.character(predicted_occupancyM$year)

# Calculate Confidence levels
predicted_occupancyM$lcl <- predicted_occupancyM$smoothed_occ - (critical_value * predicted_occupancyM$SE)
predicted_occupancyM$ucl <- predicted_occupancyM$smoothed_occ + (critical_value * predicted_occupancyM$SE)

predicted_occupancyM


#*#*#*#**#*#*#*#*#*#*#**#*#*#*#*#*#*#*#**#*#*#*#*#*#*#**#*#*#
###**Priodontes maximus***###
DetHist_Pri <- read.csv("Priodontes_25d.csv",header = TRUE, row.names="X")
effortPr <- read.csv("Priodontes_25dEf.csv",header = TRUE, row.names="X")

DetHistPrix25<- list(detection_history = DetHist_Pri, effort = effortPr)

UMFPri <- unmarked::unmarkedMultFrame(y = DetHistPrix25$detection_history,
                                      siteCovs = site_covariates[,1:14],
                                      yearlySiteCovs = list(year = years1, 
                                                            NDVI = site_covariates[,c(3:5)],
                                                            NBR = site_covariates[,c(6:8)],
                                                            AB = site_covariates[,c(9:11)],
                                                            BP = site_covariates[,c(12:14)]),
                                      obsCovs = list(effort = DetHistPrix25$effort), 
                                      numPrimary = 3)

##NuLL Model
PriNull <- unmarked::colext(psiformula = ~1,           
                            gammaformula = ~1,          
                            epsilonformula = ~1,      
                            pformula = ~1,         
                            data = UMFPri ,
                            method="BFGS")
##**DETECTION (p)*#
Prip1 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,         
                          epsilonformula = ~ 1,        
                          pformula = ~Trail ,        
                          data = UMFPri ,
                          method="BFGS")
Prip2 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,         
                          epsilonformula = ~ 1,        
                          pformula = ~Trail + effort ,         
                          data = UMFPri ,
                          method="BFGS")
Prip3 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~Trail + year ,
                          data = UMFPri ,
                          method="BFGS")

Prip4 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~effort ,
                          data = UMFPri ,
                          method="BFGS")

Prip5 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~year ,
                          data = UMFPri ,
                          method="BFGS")

Prip6 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~effort + year ,
                          data = UMFPri ,
                          method="BFGS")

Prip7 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~Trail + effort + year ,
                          data = UMFPri ,
                          method="BFGS")

modlistPri_p <-list(Prip.null=PriNull ,Prip.Trail=Prip1, Prip.Trail.Effort=Prip2,Prip.Trail.year=Prip3,
                    Prip.Effort=Prip4,Prip.year=Prip5,Prip.Effort.year=Prip6,Prip.Trail.Effort.year=Prip7)

aictab(modlistPri_p)

#**OCCUPANCY (psi)**#
Pripsi1 <- unmarked::colext(psiformula = ~Dwater,          
                            gammaformula = ~ 1,          
                            epsilonformula = ~ 1,      
                            pformula = ~effort ,         
                            data = UMFPri ,
                            method="BFGS")

Pripsi2 <- unmarked::colext(psiformula = ~NDVI1,        
                            gammaformula = ~ 1,          
                            epsilonformula = ~ 1,      
                            pformula = ~effort ,         
                            data = UMFPri ,
                            method="BFGS")

Pripsi3 <- unmarked::colext(psiformula = ~Dwater+NDVI1,            
                            gammaformula = ~ 1,          
                            epsilonformula = ~ 1,      
                            pformula = ~effort ,         
                            data = UMFPri ,
                            method="BFGS")


modlistPri_psi <-list(PriNull= PriNull,Prip.Effort=Prip4,Pripsi.Dwater.pE=Pripsi1,
                      Pripsi.NDVI.pE=Pripsi2,
                      Pripsi.Dwater.NDVI.pE=Pripsi3)

aictab(modlistPri_psi)

#**Gamma (g)**#
Pri1 <- unmarked::colext(psiformula = ~Dwater,            
                         gammaformula = ~AB,           
                         epsilonformula = ~ 1,       
                         pformula = ~effort ,         
                         data = UMFPri ,
                         method="BFGS")

ModlistPri_g <- list(PriNull= PriNull,Pripsi.Dw.pE.g1=Pripsi1,Dapsi.Dw.pE.g.AB=Pri1)
aictab(ModlistPri_g)

#**Epsilon (e)**#
Prie1 <- unmarked::colext(psiformula = ~Dwater,            
                          gammaformula = ~1,             
                          epsilonformula = ~AB,     
                          pformula = ~effort ,         
                          data = UMFPri ,
                          method="BFGS")

ModlistPri_e1 <- list(PriNull= PriNull,Pripsi.Dw.pE.g1.e1=Pripsi1,Pripsi.Dw.pE.g1.e.AB=Prie1)
aictab(ModlistPri_e1)

###GOF
mb.bootPr <- AICcmodavg::mb.gof.test(Pripsi1, nsim = 1000)
print(mb.bootPr, digit.vals = 4, digits.chisq = 4)

# Predicted occupancy in each year
mPr <- nonparboot(Pripsi1,  B = 1000)
predicted_occupancyPr <- data.frame(year = c(1:3),
                                    smoothed_occ = smoothed(Pripsi1)[2,], 
                                    SE = mPr@smoothed.mean.bsse[2,])
predicted_occupancyPr$year<-as.character(predicted_occupancyPr$year)

# Calculate Confidence levels
predicted_occupancyPr$lcl <- predicted_occupancyPr$smoothed_occ - (critical_value * predicted_occupancyPr$SE)
predicted_occupancyPr$ucl <- predicted_occupancyPr$smoothed_occ + (critical_value * predicted_occupancyPr$SE)

predicted_occupancyPr



#*#*#*#**#*#*#*#*#*#*#**#*#*#*#*#*#*#*#**#*#*#*#*#*#*#**#*#*#
###**Dicotyles tajacu***###
DetHist_Dic <- read.csv("Dicotyles_16d.csv",header = TRUE, row.names="X")
effortDi <- read.csv("Dicotyles_16dEf.csv",header = TRUE, row.names="X")

DetHistDicx16<- list(detection_history = DetHist_Dic, effort = effortDi)

UMFDic <- unmarked::unmarkedMultFrame(y = DetHistDicx16$detection_history,
                                      siteCovs = site_covariates[,1:14],
                                      yearlySiteCovs = list(year = years1, 
                                                            NDVI = site_covariates[,c(3:5)],
                                                            NBR = site_covariates[,c(6:8)],
                                                            AB = site_covariates[,c(9:11)],
                                                            BP = site_covariates[,c(12:14)]),
                                      obsCovs = list(effort = DetHistDicx16$effort), 
                                      numPrimary = 3)

##NuLL Model
DicNull <- unmarked::colext(psiformula = ~1,           
                            gammaformula = ~1,          
                            epsilonformula = ~1,      
                            pformula = ~1,         
                            data = UMFDic,
                            method="BFGS")
##Now DETECTION (p)
Dicp1 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~Trail ,       
                          data = UMFDic,
                          method="BFGS")
Dicp2 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~Trail + effort ,        
                          data = UMFDic,
                          method="BFGS")
Dicp3 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~Trail + year ,
                          data = UMFDic,
                          method="BFGS")

Dicp4 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~effort ,
                          data = UMFDic,
                          method="BFGS")

Dicp5 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~year ,
                          data = UMFDic,
                          method="BFGS")

Dicp6 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~effort + year ,
                          data = UMFDic,
                          method="BFGS")

Dicp7 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~Trail + effort + year ,
                          data = UMFDic,
                          method="BFGS")

modlistDic_p <-list(Dicp.null=DicNull ,Dicp.Trail=Dicp1, Dicp.Trail.Effort=Dicp2,Dicp.Trail.year=Dicp3,
                    Dicp.Effort=Dicp4,Dicp.year=Dicp5,Dicp.Effort.year=Dicp6,Dicp.Trail.Effort.year=Dicp7)
aictab(modlistDic_p)

#**OCCUPANCY (psi)**#
Dicpsi1 <- unmarked::colext(psiformula = ~Dwater,            
                            gammaformula = ~ 1,          
                            epsilonformula = ~ 1,      
                            pformula = ~Trail + effort + year ,         
                            data = UMFDic ,
                            method="BFGS")

Dicpsi2 <- unmarked::colext(psiformula = ~NDVI1,            
                            gammaformula = ~ 1,          
                            epsilonformula = ~ 1,      
                            pformula = ~Trail + effort + year ,         
                            data = UMFDic ,
                            method="BFGS")

Dicpsi3 <- unmarked::colext(psiformula = ~Dwater+NDVI1,            
                            gammaformula = ~ 1,          
                            epsilonformula = ~ 1,      
                            pformula = ~Trail + effort + year ,         
                            data = UMFDic ,
                            method="BFGS")

modlistDic_psi <-list(DicNull= DicNull,Dicpsi1.p.T.E.y=Dicp2,Dicpsi.Dwater.p.T.E.y=Dicpsi1,
                      Dicpsi.NDVI.p.T.E.y=Dicpsi2,
                      Dicpsi.Dwater.NDVI.p.T.E.y=Dicpsi3)

aictab(modlistDic_psi)

#**Gamma (g)**#
Dicg1 <- unmarked::colext(psiformula = ~Dwater,            
                          gammaformula = ~AB,           
                          epsilonformula = ~ 1,       
                          pformula = ~Trail + effort ,         
                          data = UMFDic ,
                          method="BFGS")

ModlistDic_g <- list(Dic.Null= DicNull,Dic.psi.Dw.p.T.E.y.g1=Dicpsi1,Dic.psi.Dw.p.T.E.y.g.AB=Dicg1)
aictab(ModlistDic_g)

#**Now Epsilon (e)**#
Dice1 <- unmarked::colext(psiformula = ~Dwater,            
                          gammaformula = ~1,             
                          epsilonformula = ~AB,     
                          pformula = ~Trail + effort ,         
                          data = UMFDic ,
                          method="BFGS")

ModlistDic_e <- list(Dic.Null= DicNull,Dic.psi.Dw.p.T.E.y.g1=Dicpsi1,Dic.psi.Dw.p.T.E.y.g1.e.AB=Dice1)
aictab(ModlistDic_e)

###GOF
mb.bootDi <- AICcmodavg::mb.gof.test(Dicpsi1, nsim = 1000)
print(mb.bootDi, digit.vals = 4, digits.chisq = 4)

# Predicted occupancy in each year
mDi <- nonparboot(Dicpsi1,  B = 1000)
predicted_occupancyDi <- data.frame(year = c(1:3),
                                    smoothed_occ = smoothed(Dicpsi1)[2,], 
                                    SE = mDi@smoothed.mean.bsse[2,])
predicted_occupancyDi$year<-as.character(predicted_occupancyDi$year)

# Calculate Confidence levels
predicted_occupancyDi$lcl <- predicted_occupancyDi$smoothed_occ - (critical_value * predicted_occupancyDi$SE)
predicted_occupancyDi$ucl <- predicted_occupancyDi$smoothed_occ + (critical_value * predicted_occupancyDi$SE)

predicted_occupancyDi


#*#*#*#**#*#*#*#*#*#*#**#*#*#*#*#*#*#*#**#*#*#*#*#*#*#**#*#*#
###**Daziprocta azarae***###
DetHist_Daz <- read.csv("Dasyprocta_25d.csv",header = TRUE, row.names="X")
effortDa <- read.csv("Dasyprocta_25dEf.csv",header = TRUE, row.names="X")

DetHistDazx25<- list(detection_history = DetHist_Daz, effort = effortDa)

UMFDaz <- unmarked::unmarkedMultFrame(y = DetHistDazx25$detection_history,
                                      siteCovs = site_covariates[,1:14],
                                      yearlySiteCovs = list(year = years1, 
                                                            NDVI = site_covariates[,c(3:5)],
                                                            NBR = site_covariates[,c(6:8)],
                                                            AB = site_covariates[,c(9:11)],
                                                            BP = site_covariates[,c(12:14)]),
                                      obsCovs = list(effort = DetHistDazx25$effort), 
                                      numPrimary = 3)

##NuLL Model
DazNull <- unmarked::colext(psiformula = ~1,          
                            gammaformula = ~1,          
                            epsilonformula = ~1,      
                            pformula = ~1,         
                            data = UMFDaz ,
                            method="BFGS")
##**DETECTION (p)**
Dazp1 <- unmarked::colext(psiformula = ~ 1,          
                          gammaformula = ~ 1,         
                          epsilonformula = ~ 1,       
                          pformula = ~Trail ,       
                          data = UMFDaz ,
                          method="BFGS")
Dazp2 <- unmarked::colext(psiformula = ~ 1,        
                          gammaformula = ~ 1,        
                          epsilonformula = ~ 1,      
                          pformula = ~Trail + effort ,        
                          data = UMFDaz ,
                          method="BFGS")
Dazp3 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~Trail + year ,
                          data = UMFDaz ,
                          method="BFGS")

Dazp4 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~effort ,
                          data = UMFDaz ,
                          method="BFGS")

Dazp5 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~year ,
                          data = UMFDaz ,
                          method="BFGS")

Dazp6 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~effort + year ,
                          data = UMFDaz ,
                          method="BFGS")

Dazp7 <- unmarked::colext(psiformula = ~ 1,           
                          gammaformula = ~ 1,          
                          epsilonformula = ~ 1,        
                          pformula = ~Trail + effort + year ,
                          data = UMFDaz ,
                          method="BFGS")

modlistDaz_p <-list(Dazp.null=DazNull,Dazp.Trail=Dazp1, Dazp.Trail.Effort=Dazp2,Dazp.Trail.year=Dazp3,
                    Dazp.Effort=Dazp4,Dazp.year=Dazp5,Dazp.Effort.year=Dazp6,Dazp.Trail.Effort.year=Dazp7)

aictab(modlistDaz_p)

#**OCCUPANCY (psi)**#
Dazpsi1 <- unmarked::colext(psiformula = ~Dwater,           
                            gammaformula = ~ 1,          
                            epsilonformula = ~ 1,      
                            pformula = ~effort+year ,         
                            data = UMFDaz ,
                            method="BFGS")

Dazpsi2 <- unmarked::colext(psiformula = ~NDVI1,           
                            gammaformula = ~ 1,          
                            epsilonformula = ~ 1,      
                            pformula = ~effort+year ,         
                            data = UMFDaz ,
                            method="BFGS")

Dazpsi3 <- unmarked::colext(psiformula = ~Dwater+NDVI1,            
                            gammaformula = ~ 1,          
                            epsilonformula = ~ 1,      
                            pformula = ~effort+year ,         
                            data = UMFDaz ,
                            method="BFGS")

modlistDaz_psi <-list(DazNull= DazNull,Dazpsi1.p.E.y=Dazp6,Dazpsi.Dwater.p.E.y=Dazpsi1,
                      Dazpsi.NDVI.p.E.y=Dazpsi2,
                      Dazpsi.Dwater.NDVI.p.E.y=Dazpsi3)
aictab(modlistDaz_psi)

#**Gamma (g)**#
Daz1 <- unmarked::colext(psiformula = ~1,            
                         gammaformula = ~AB,            
                         epsilonformula = ~ 1,       
                         pformula = ~effort+year ,         
                         data = UMFDaz ,
                         method="BFGS")

ModlistDaz_g <- list(DazNull= DazNull,Dazpsi1.p.E.y=Dazp6,Dapsi1.p.E.y.g.AB=Daz1)
aictab(ModlistDaz_g)

#**Epsilon (e)**#
Daze1 <- unmarked::colext(psiformula = ~1,            
                          gammaformula = ~1,             
                          epsilonformula = ~AB,    
                          pformula = ~effort ,         
                          data = UMFDaz ,
                          method="BFGS")

ModlistDaz_e <- list(DazNull= DazNull,Dazpsi1.p.E.y=Dazp6,Dazpsi1.p.E.y.g1.e.AB=Daze1)
aictab(ModlistDaz_e)

###GOF
mb.bootDa <- AICcmodavg::mb.gof.test(Dazp6, nsim = 1000)
print(mb.bootDa, digit.vals = 4, digits.chisq = 4)

# Predicted occupancy in each year
mDa <- nonparboot(Dazp6,  B = 1000)
predicted_occupancyDa <- data.frame(year = c(1:3),
                                    smoothed_occ = smoothed(Dazp6)[2,], 
                                    SE = mDa@smoothed.mean.bsse[2,])
predicted_occupancyDa$year<-as.character(predicted_occupancyDi$year)

# Calculate Confidence levels
predicted_occupancyDa$lcl <- predicted_occupancyDa$smoothed_occ - (critical_value * predicted_occupancyDa$SE)
predicted_occupancyDa$ucl <- predicted_occupancyDa$smoothed_occ + (critical_value * predicted_occupancyDa$SE)

predicted_occupancyDa
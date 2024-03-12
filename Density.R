#******************************************************************
#######**MULTI-season DENSITY ANALYSIS########*******###
#**************Jaguars & Ocelots************************************

#####**Packages**####
library(secr)

###Getting data-Jaguars####
Sess1 <- read.capthist("CaptHist_jags1.csv", "Effort_s1.csv", covnames = c("Sex", "Site", "Prey"), detector = "count", fmt = "trapID")                                              
Sess2 <- read.capthist("CaptHist_jags2.csv", "Effort_s2.csv", covnames = c("Sex", "Site", "Prey"), detector = "count", fmt = "trapID")                                              
Sess3 <- read.capthist("CaptHist_jags3.csv", "Effort_s3.csv", covnames = c("Sex", "Site", "Prey"), detector = "count", fmt = "trapID")

mask_1 <- read.mask("Mask1s.csv", spacing = 500, header = TRUE, columns = c("NDVI", "dWater", "NBR"))
mask_2 <- read.mask("Mask2s.csv", spacing = 500, header = TRUE, columns = c("NDVI", "dWater", "NBR"))
mask_3 <- read.mask("Mask3s.csv", spacing = 500, header = TRUE, columns = c("NDVI", "dWater", "NBR"))

mask <- list('Sess1' = mask_1, 'Sess2' = mask_2, 'Sess3' = mask_3)
lapply(mask, function(mask_i) {covariates(mask_i) <- scale(covariates(mask_i))})

All_Jag <- MS.capthist(Sess1, Sess2, Sess3)
summary(All_Jag, terse = TRUE)

####running models-Jaguars####
#1° step
fits3a <- secr.fit(All_Jag,model=list(D~1,g0~Site,sigma~session),mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead', CL=TRUE)
fits3b <- secr.fit(All_Jag,model=list(D~1,g0~b,sigma~session),mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead', CL=TRUE)
fits3c <- secr.fit(All_Jag,model=list(D~1,g0~1,sigma~session),mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead', CL=TRUE)
fits3d <- secr.fit(All_Jag,model=list(D~1,g0~Site+b,sigma~session),mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead', CL=TRUE)

AICsJag_ms1=AIC(fits3a,fits3b,fits3c,fits3d)
AICsJag_ms1

##2°step
fits3b1 <- secr.fit(All_Jag,model=list(D~NDVI,g0~1,sigma~session),mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead')
fits3c1 <- secr.fit(All_Jag,model=list(D~dWater,g0~1,sigma~session),mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead')
fits3d1 <- secr.fit(All_Jag,model=list(D~NDVI+dWater,g0~1,sigma~session),mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead')
fits3a1 <- secr.fit(All_Jag,model=list(D~NBR,g0~1,sigma~session),mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead')
fits3e1 <- secr.fit(All_Jag,model=list(D~session,g0~1,sigma~session),mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead')
fits3f1 <- secr.fit(All_Jag,model=list(D~NBR+dWater,g0~1,sigma~session),mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead')
fits3g1 <- secr.fit(All_Jag,model=list(D~NBR+dWater+NDVI,g0~1,sigma~session),mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead')

AicsJagMS=AIC(fits3e1,fits3a1,fits3b1,fits3c1,fits3d1,fits3f1)
AicsJagMS

#test against null
fitsNULL2 <- secr.fit(All_Jag,model=list(D~1,g0~1,sigma~session),mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead')
AicsJag3=AIC(fits3b1,fits3c1,fitsNULL2)
AicsJag3

#getting values
density <- DENSITY(fits3c1, method = "average")
derived(fits3c1)
predict(fits3c1)
beta <- coef(fits3c1)[1:length(coef(fits3c1))]

###Getting data-Ocelots####
Sess1b=read.capthist("CaptHist_ocs1.csv","Effort_s1.csv",covnames =c("Sex","Site"),detector="count",fmt="trapID")                                              
Sess2b=read.capthist("CaptHist_ocs2.csv", "Effort_s2.csv",covnames =c("Sex","Site"),detector="count",fmt="trapID")                                              
Sess3b=read.capthist("CaptHist_ocs3.csv", "Effort_s3.csv",covnames =c("Sex","Site"),detector="count",fmt="trapID")

mask_1<-read.mask("Mask1s.csv", spacing = 500, header=TRUE,columns=c("NDVI","dWater","NBR"))
mask_2<-read.mask("Mask2s.csv", spacing = 500, header=TRUE,columns=c("NDVI","dWater","NBR"))
mask_3<-read.mask("Mask3s.csv", spacing = 500, header=TRUE,columns=c("NDVI","dWater","NBR"))
mask_Oc=list('Sess1b' =mask_1,'Sess2b' =mask_2,'Sess3b' =mask_3)

lapply(mask_Oc, function(mask_i) {covariates(mask_i) <- scale(covariates(mask_i))})

All_Oc=MS.capthist(Sess1b,Sess2b,Sess3b)
summary(All_Oc, terse = TRUE)

####running models - Ocelots####
fits1a <- secr.fit(All_Oc,model=list(D~1,g0~Site,sigma~session),mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead', CL=TRUE)
fits1b <- secr.fit(All_Oc,model=list(D~1,g0~b,sigma~session),mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead', CL=TRUE)
fits1c <- secr.fit(All_Oc,model=list(D~1,g0~1,sigma~session),mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead', CL=TRUE)
fits1d <- secr.fit(All_Oc,model=list(D~1,g0~Site+b,sigma~session),mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead', CL=TRUE)

AICsOc_ms1=AIC(fits1a,fits1b,fits1c,fits1d)
AICsOc_ms1

fits2b1 <- secr.fit(All_Oc,model=list(D~NDVI,g0~h2,sigma~session),hcov="Site",mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead')
fits2c1 <- secr.fit(All_Oc,model=list(D~dWater,g0~h2,sigma~session),hcov="Site",mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead')
fits2d1 <- secr.fit(All_Oc,model=list(D~NDVI+dWater,g0~h2,sigma~session),hcov="Site",mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead')
fits2a1 <- secr.fit(All_Oc,model=list(D~NBR,g0~h2,sigma~session),hcov="Site",mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead')
fits2e1 <- secr.fit(All_Oc,model=list(D~session,g0~h2,sigma~session),hcov="Site",mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead')
fits2f1 <- secr.fit(All_Oc,model=list(D~NBR+dWater,g0~h2,sigma~session),hcov="Site",mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead')

AicsOcMS=AIC(fits2a1,fits2b1,fits2c1,fits2d1,fits2e1,fits2f1)
AicsOcMS

fitsNULL2 <- secr.fit(All_Oc,model=list(D~1,g0~h2,sigma~session),hcov="Site",mask=mask, trace = FALSE,verify=FALSE,method = 'Nelder-Mead')
AicsOc3=AIC(fits2b1,fitsNULL2)
AicsOc3 

coef(fits2b1)
derived(fits2b1)
predict(fits2b1)
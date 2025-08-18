
####Code for landscape compositions withing buffers of 300 m radius 

library(spOccupancy)
library(datasets)
library(base)
library(MCMCvis)


# 1. Defining matrices Preparing the data -----------------------------------------------------------

y<-as.matrix(occurrance_300c)
y<-data.matrix(occurrance_300c)
str(y)
occ.covs<- as.matrix(cov_occ_300c)
#occ.covs<- as.data.frame(cov_occ_300)
str(occ.covs)
precip<- as.matrix(precip_300_c)
tempmax<- as.matrix(temp_300_b)
tempmin<- as.matrix(MinTemp)
tempmean<- as.matrix(MeanTemp)

det.covs<-list(precip = precip, 
               tempmax = tempmax,
               tempmin= tempmin,
               tempmean=tempmean)
str(det.covs)

coords<- as.matrix(cord_UTM_300)
str(coords)

data.list <- list(y = y, 
                  occ.covs = occ.covs, 
                  det.covs = det.covs, 
                  coords = coords)

plot(data.list$coords, pch=19)

str(data.list)



# 2. Models -----------------------------------------------------------

####Detection

##Linear equation: logit(pij​)=α0​+α1​⋅tempmean​+α2​⋅precip
​
btbw.occ.formula1.300<- ~ 1
btbw.det.formula <- ~ scale(tempmean)+scale(precip)
out1.300 <- spPGOcc(occ.formula = btbw.occ.formula1.300,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 4000, 
                    batch.length = 25,
                    NNGP = TRUE, 
                    cov.model = "exponential", 
                    n.neighbors = 15, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out1.300)

####Single variables

###Successional forest

####linear equation: logit(ψi​)=β0​+β1​⋅Successional_forest​+β2​⋅Successional_forest^2


btbw.occ.formula1.300<- ~ I(scale(`Sucecional Forest`)^2)
btbw.det.formula <- ~ scale(tempmax)
out1.300 <- spPGOcc(occ.formula = btbw.occ.formula1.300,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 4000, 
                    batch.length = 25,
                    NNGP = TRUE, 
                    cov.model = "exponential", 
                    n.neighbors = 15, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out1.300)
ppc.out1.300 <- ppcOcc(out1.300, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out1.300)
waicOcc(out1.300)
names(out1.300)
str(out1.300$beta.samples)
MCMCplot(out1.300$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "successional forest"),ci=c(50, 95))

##No spatial, Successional forest 

btbw.occ.formula1b <- ~ scale(`Sucecional Forest`)
btbw.det.formula <- ~ scale(tempmax)
n.samples <- 5000
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72))
inits.list <- list(alpha = 0, beta = 0,
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
n.report <- 1000
out1b <- PGOcc(occ.formula = btbw.occ.formula1b, 
               det.formula = btbw.det.formula, 
               data = data.list, 
               inits = inits.list,
               n.samples = n.samples,
               priors = prior.list,
               n.omp.threads = 1,
               verbose = TRUE,
               n.report = n.report, 
               n.burn = 1000, 
               n.thin = 1, 
               n.chains = 40)
summary(out1b)
ppc.out1b <- ppcOcc(out1b, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out1b)
waicOcc(out1b)

##### Native forest
####linear equation: logit(ψi​)=β0​+β1​⋅Native_forest+β2​⋅Native_forest^2

btbw.occ.formula2.300<- ~ I(scale(`Native forest`)^2)
btbw.det.formula <- ~ scale(tempmax)
out2.300 <- spPGOcc(occ.formula = btbw.occ.formula2.300,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 4000, 
                    batch.length = 25,
                    NNGP = TRUE, 
                    cov.model = "exponential", 
                    n.neighbors = 15, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out2.300)
ppc.out2.300 <- ppcOcc(out2.300, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out2.300)
waicOcc(out2.300)
names(out2.300)
str(out2.300$beta.samples)
MCMCplot(out2.300$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "Native forest"),ci=c(50, 95))

##No spatial, Native forest 

btbw.occ.formula2b <- ~ scale(`Native forest`)
btbw.det.formula <- ~ scale(tempmax)
n.samples <- 5000

prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72))

inits.list <- list(alpha = 0, beta = 0,
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
n.report <- 1000

out2b <- PGOcc(occ.formula = btbw.occ.formula2b, 
               det.formula = btbw.det.formula, 
               data = data.list, 
               inits = inits.list,
               n.samples = n.samples,
               priors = prior.list,
               n.omp.threads = 1,
               verbose = TRUE,
               n.report = n.report, 
               n.burn = 1000, 
               n.thin = 1, 
               n.chains = 40)

summary(out2b)
ppc.out2b <- ppcOcc(out2b, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out2b)
waicOcc(out2b)

### Grassland
####linear equation: logit(ψi​)=β0​+β1​⋅Grassland+β2​⋅Grassland^2

btbw.occ.formula3.300<- ~ I(scale(Grassland)^2)
btbw.det.formula <- ~ scale(tempmax)
out3.300 <- spPGOcc(occ.formula = btbw.occ.formula3.300,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 4000, 
                    batch.length = 25,
                    NNGP = TRUE, 
                    cov.model = "exponential", 
                    n.neighbors = 15, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out3.300)
ppc.out3.300 <- ppcOcc(out3.300, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out3.300)
waicOcc(out3.300)
names(out3.300)
str(out3.300$beta.samples)
MCMCplot(out3.300$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "Grassland"),ci=c(50, 95))

##No spatial, Grassland 

btbw.occ.formula3b <- ~ scale(Grassland)
btbw.det.formula <- ~ scale(tempmax)
n.samples <- 5000

prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72))

inits.list <- list(alpha = 0, beta = 0,
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
n.report <- 1000

out3b <- PGOcc(occ.formula = btbw.occ.formula3b, 
               det.formula = btbw.det.formula, 
               data = data.list, 
               inits = inits.list,
               n.samples = n.samples,
               priors = prior.list,
               n.omp.threads = 1,
               verbose = TRUE,
               n.report = n.report, 
               n.burn = 1000, 
               n.thin = 1, 
               n.chains = 40)

summary(out3b)
ppc.out3b <- ppcOcc(out3b, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out3b)
waicOcc(out3b)

### Cropland
####linear equation: logit(ψi​)=β0​+β1​⋅Cropland+β2​⋅Cropland^2

btbw.occ.formula4.300<- ~ scale(Cropland)
btbw.det.formula <- ~ scale(tempmax)
out4.300 <- spPGOcc(occ.formula = btbw.occ.formula4.300,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 4000, 
                    batch.length = 25,
                    NNGP = TRUE, 
                    cov.model = "exponential", 
                    n.neighbors = 15, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out4.300)
ppc.out4.300 <- ppcOcc(out4.300, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out4.300)
waicOcc(out4.300)
names(out4.300)
str(out4.300$beta.samples)
MCMCplot(out4.300$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "Cropland"),ci=c(50, 95))

### Elevation
####linear equation: logit(ψi​)=β0​+β1​⋅Elevation+β2​⋅Elevation^2

btbw.occ.formula5.300<- ~ I(scale(Elevation)^2)
btbw.det.formula <- ~ scale(tempmax)
out5.300 <- spPGOcc(occ.formula = btbw.occ.formula5.300,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 4000, 
                    batch.length = 25,
                    NNGP = TRUE, 
                    cov.model = "exponential", 
                    n.neighbors = 15, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out5.300)
ppc.out5.300 <- ppcOcc(out5.300, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out5.300)
waicOcc(out5.300)
names(out5.300)
str(out5.300$beta.samples)
MCMCplot(out5.300$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "Elevation"),ci=c(50, 90))

##No spatial, Cropland 

btbw.occ.formula5b <- ~ I(scale(Cropland)^2)
btbw.det.formula <- ~ scale(tempmax)
n.samples <- 5000

prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72))

inits.list <- list(alpha = 0, beta = 0,
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
n.report <- 1000

out5b <- PGOcc(occ.formula = btbw.occ.formula5b, 
               det.formula = btbw.det.formula, 
               data = data.list, 
               inits = inits.list,
               n.samples = n.samples,
               priors = prior.list,
               n.omp.threads = 1,
               verbose = TRUE,
               n.report = n.report, 
               n.burn = 1000, 
               n.thin = 1, 
               n.chains = 40)

summary(out5b)
ppc.out5b <- ppcOcc(out5b, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out5b)
waicOcc(out5b)


###
####linear equation: logit(ψi​)=β0​+β1​⋅Stem height+β2​⋅Stem height^2

btbw.occ.formula6.300<- ~ I(scale(`Stem height`)^2)
btbw.det.formula <- ~ scale(tempmax)
out6.300 <- spPGOcc(occ.formula = btbw.occ.formula6.300,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 4000, 
                    batch.length = 25,
                    NNGP = TRUE, 
                    cov.model = "exponential", 
                    n.neighbors = 15, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out6.300)
ppc.out6.300 <- ppcOcc(out6.300, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out6.300)
waicOcc(out6.300)
names(out6.300)
str(out6.300$beta.samples)
MCMCplot(out6.300$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "Stem height"),ci=c(50, 89))

###
####linear equation: logit(ψi​)=β0​+β1​⋅DOM score+β2​⋅DOM score^2
btbw.occ.formula7.300<- ~ scale(`DOM score`)
btbw.det.formula <- ~ scale(tempmax)
out7.300 <- spPGOcc(occ.formula = btbw.occ.formula7.300,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 4000, 
                    batch.length = 25,
                    NNGP = TRUE, 
                    cov.model = "exponential", 
                    n.neighbors = 15, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out7.300)
ppc.out7.300 <- ppcOcc(out7.300, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out7.300)
waicOcc(out7.300)
names(out7.300)
str(out7.300$beta.samples)
MCMCplot(out7.300$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "DOM score"),ci=c(50, 89))

###
####linear equation: logit(ψi​)=β0​+β1​⋅infrutescence+β2​⋅infrutescence^2
btbw.occ.formula8.300<- ~ I(scale(infrutescence)^2)
btbw.det.formula <- ~ scale(tempmax)
out8.300 <- spPGOcc(occ.formula = btbw.occ.formula8.300,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 4000, 
                    batch.length = 25,
                    NNGP = TRUE, 
                    cov.model = "exponential", 
                    n.neighbors = 15, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out8.300)
ppc.out8.300 <- ppcOcc(out8.300, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out8.300)
waicOcc(out8.300)
names(out8.300)
str(out8.300$beta.samples)
MCMCplot(out8.300$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "infrutescence"),ci=c(50, 90))

###
####linear equation: logit(ψi​)=β0​+β1​⋅connected_trees+β2​⋅connected_trees^2
btbw.occ.formula9.300<- ~ scale(connected_trees)
btbw.det.formula <- ~ scale(precip)
out9.300 <- spPGOcc(occ.formula = btbw.occ.formula9.300,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 4000, 
                    batch.length = 25,
                    NNGP = TRUE, 
                    cov.model = "exponential", 
                    n.neighbors = 15, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out9.300)
ppc.out9.300 <- ppcOcc(out9.300, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out9.300)
waicOcc(out9.300)
names(out9.300)
str(out9.300$beta.samples)
MCMCplot(out9.300$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "Connected Trees"),ci=c(50, 90))

###DOM_score

btbw.occ.formula10.300<- ~ I(scale(DOM_score)^2)
btbw.det.formula <- ~ scale(tempmax)
out10.300 <- spPGOcc(occ.formula = btbw.occ.formula10.300,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 4000, 
                    batch.length = 25,
                    NNGP = TRUE, 
                    cov.model = "exponential", 
                    n.neighbors = 15, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out10.300)
ppc.out10.300 <- ppcOcc(out10.300, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out10.300)
waicOcc(out10.300)
names(out10.300)
str(out10.300$beta.samples)
MCMCplot(out10.300$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "DOM score"),ci=c(50, 90))

waicOcc(out1.300)
waicOcc(out2.300)
waicOcc(out3.300)
waicOcc(out4.300)
waicOcc(out5.300)
waicOcc(out6.300)
waicOcc(out7.300)
waicOcc(out8.300)
waicOcc(out9.300)
waicOcc(out4a)


######Additive models

###Linear equation: logit(ψi​)=β0​+β1​⋅Sucecional Forest^2​+β2​⋅Grassland^2​+β3​⋅infrutescence

btbw.occ.formula10.300 <- ~ I(scale(`Sucecional Forest`)^2) + I(scale(Grassland))^2 + scale(infrutescence)
btbw.det.formula <- ~ scale(precip)
out10.300  <- spPGOcc(occ.formula = btbw.occ.formula10.300 ,
                      det.formula = btbw.det.formula,
                      data = data.list, 
                      n.batch = 400, 
                      batch.length = 25,
                      NNGP = TRUE, 
                      n.neighbors = 15, 
                      n.burn = 2000, 
                      n.thin = 4, 
                      n.chains = 3, 
                      verbose = FALSE, 
                      k.fold = 2)

summary(out10.300)
ppc.out10.300 <- ppcOcc(out10.300, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out10.300)
waicOcc(out10.300)
names(out10.300)
str(out10.300$beta.samples)
MCMCplot(out10.300$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), ci=c(50, 95))
MCMCplot(out10.300$alpha.samples, ref_ovl=TRUE, ci=c(50, 95))

####### Models with only land Use covers

###Linear equation: logit(ψi​)=β0​+β1​⋅Sucecional Forest^2​+β2​⋅Grassland^2​+β3​⋅Cropland^2​+β4​⋅Artificial^2​
btbw.occ.formula11.300 <- ~ I(scale(`Sucecional Forest`)^2) + I(scale(Grassland))^2 + I(scale(Cropland)^2) + I(scale(Artificial)^2) 
btbw.det.formula <- ~ scale(precip)
out11.300  <- spPGOcc(occ.formula = btbw.occ.formula11.300 ,
                      det.formula = btbw.det.formula,
                      data = data.list, 
                      n.batch = 400, 
                      batch.length = 25,
                      NNGP = TRUE, 
                      n.neighbors = 15, 
                      n.burn = 2000, 
                      n.thin = 4, 
                      n.chains = 3, 
                      verbose = FALSE, 
                      k.fold = 2)

summary(out11.300)
ppc.out11.300 <- ppcOcc(out11.300, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out11.300)
waicOcc(out11.300)
names(out11.300)
str(out11.300$beta.samples)
MCMCplot(out11.300$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), ci=c(50, 95))
MCMCplot(out11.300$alpha.samples, ref_ovl=TRUE, ci=c(50, 95))



####
####Linear Equation: logit(ψi​)=β0​+β1​⋅Sucecional Forest^2​+β2​⋅Grassland^2​+β3​⋅Cropland^2​+β4​⋅Artificial^2​+β5​⋅infrutescence​
btbw.occ.formula12.300 <- ~ I(scale(`Sucecional Forest`)^2) + I(scale(Grassland))^2 + I(scale(Cropland)^2) + I(scale(Artificial)^2) + scale(infrutescence)
btbw.det.formula <- ~ scale(precip)
out12.300  <- spPGOcc(occ.formula = btbw.occ.formula12.300 ,
                      det.formula = btbw.det.formula,
                      data = data.list, 
                      n.batch = 400, 
                      batch.length = 25,
                      NNGP = TRUE, 
                      n.neighbors = 15, 
                      n.burn = 2000, 
                      n.thin = 4, 
                      n.chains = 3, 
                      verbose = FALSE, 
                      k.fold = 2)

summary(out12.300)
ppc.out12.300 <- ppcOcc(out12.300, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out12.300)
waicOcc(out12.300)
names(out12.300)
str(out12.300$beta.samples)
MCMCplot(out12.300$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), ci=c(50, 95))
MCMCplot(out12.300$alpha.samples, ref_ovl=TRUE, ci=c(50, 95))




####
###Linear equation: logit(ψi​)=β0​+β1​⋅Sucecional Forest^2​+β2​⋅Grassland^2​+β3​⋅Cropland^2​+β4​⋅Artificial^2​+β5​⋅infrutescence​+β6​⋅DOM_score^2​
btbw.occ.formula13.300 <- ~ I(scale(`Sucecional Forest`)^2) + I(scale(Grassland))^2 + I(scale(Cropland)^2) + I(scale(Artificial)^2) + scale(infrutescence) + I(scale(DOM_score)^2)
btbw.det.formula <- ~ scale(precip)
out13.300  <- spPGOcc(occ.formula = btbw.occ.formula13.300 ,
                      det.formula = btbw.det.formula,
                      data = data.list, 
                      n.batch = 400, 
                      batch.length = 25,
                      NNGP = TRUE, 
                      n.neighbors = 15, 
                      n.burn = 2000, 
                      n.thin = 4, 
                      n.chains = 3, 
                      verbose = FALSE, 
                      k.fold = 2)

summary(out13.300)
ppc.out13.300 <- ppcOcc(out13.300, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out13.300)
waicOcc(out13.300)
names(out13.300)
str(out13.300$beta.samples)
MCMCplot(out13.300$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), ci=c(50, 95))
MCMCplot(out13.300$alpha.samples, ref_ovl=TRUE, ci=c(50, 95))


####
###Linear equation: logit(ψi​)=β0​+β1​⋅Sucecional Forest^2​+β2​⋅Grassland^2​+β3​⋅Cropland^2​+β4​⋅Artificial^2​+β5​⋅infrutescence​+β6​⋅DOM_score^2​+β7​⋅Altitude (feet)^2​

btbw.occ.formula14.300 <- ~ I(scale(`Sucecional Forest`)^2) + I(scale(Grassland))^2 + I(scale(Cropland)^2) + I(scale(Artificial)^2) + scale(infrutescence) + I(scale(DOM_score)^2) + I(scale(`Altitude (feet)`)^2)
btbw.det.formula <- ~ scale(precip)
out14.300  <- spPGOcc(occ.formula = btbw.occ.formula14.300 ,
                      det.formula = btbw.det.formula,
                      data = data.list, 
                      n.batch = 400, 
                      batch.length = 25,
                      NNGP = TRUE, 
                      n.neighbors = 15, 
                      n.burn = 2000, 
                      n.thin = 4, 
                      n.chains = 3, 
                      verbose = FALSE, 
                      k.fold = 2)

summary(out14.300)
ppc.out14.300 <- ppcOcc(out14.300, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out14.300)
waicOcc(out14.300)
names(out14.300)
str(out14.300$beta.samples)
MCMCplot(out14.300$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), ci=c(50, 95))
MCMCplot(out14.300$alpha.samples, ref_ovl=TRUE, ci=c(50, 95))


####
####Linear equation: logit(ψi​)=β0​+β1​⋅(Successional Forest​)^2+β2​⋅(Grassland​)^2+β3​⋅(Cropland)^2+β4​⋅(Artificial Land Cover)^2+β5​⋅Infrutescence​+β6​⋅(DOM Score​)^2+β7​⋅(Altitude​)^2

btbw.occ.formula15.300 <- ~ I(scale(`Sucecional Forest`)^2) + I(scale(Grassland))^2 + I(scale(Cropland)^2) + I(scale(Artificial)^2) + scale(infrutescence) + I(scale(DOM_score)^2) + I(scale(`Altitude (feet)`)^2)
btbw.det.formula <- ~ scale(precip)
out14.300  <- spPGOcc(occ.formula = btbw.occ.formula14.300 ,
                      det.formula = btbw.det.formula,
                      data = data.list, 
                      n.batch = 400, 
                      batch.length = 25,
                      NNGP = TRUE, 
                      n.neighbors = 15, 
                      n.burn = 2000, 
                      n.thin = 4, 
                      n.chains = 3, 
                      verbose = FALSE, 
                      k.fold = 2)

summary(out14.300)
ppc.out14.300 <- ppcOcc(out14.300, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out14.300)
waicOcc(out14.300)
names(out14.300)
str(out14.300$beta.samples)
MCMCplot(out14.300$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), ci=c(50, 95))
MCMCplot(out14.300$alpha.samples, ref_ovl=TRUE, ci=c(50, 95))



####
####Linear equation: logit(ψi​)=β0​+β1​⋅(Successional Forest)^2+β2​⋅(Grassland)^2+β3​⋅Cropland+β4​⋅Artificial Land Cover+β5​⋅(DOM Score)^2+β6​⋅Infrutescence
btbw.occ.formula4.2.1 <- ~ I(scale(`Sucecional Forest`)^2) + I(scale(Grassland)^2) + scale (Cropland) + scale (Artificial) +  I(scale(DOM_score)^2) + scale(infrutescence)
btbw.det.formula <- ~ scale(precip)
out4.2.1 <- spPGOcc(occ.formula = btbw.occ.formula4.2.1,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 4000, 
                    batch.length = 25,
                    NNGP = TRUE, 
                    cov.model = "exponential", 
                    n.neighbors = 15, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out4.2.1)
ppc.out4.2.1 <- ppcOcc(out4.2.1, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out4.2.1)
waicOcc(out4.2.1)
names(out4.2.1)
str(out4.2.1$beta.samples)
MCMCplot(out4.2.1$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1),labels= c("Intercept"," % Sucecional forest","% Grassland", "% Cropland", "% Artificial", "DOM score", "Infrutescence"), ci=c(50, 89))

MCMCplot(out4.2.1$alpha.samples, ref_ovl=TRUE, ci=c(50, 89))

####The best so so far! 2024 with Precip

###Linear equation: logit(ψi​)=β0​+β1​⋅(Successional Forest)^2+β2​⋅(Native Forest)^2+β3​⋅(Grassland)^2+β4​⋅(Cropland)^2+β5​⋅(DOM Score)^2+β6​⋅Infrutescence+β7​⋅(Altitude)^2+β8​⋅Connected Trees
btbw.occ.formula4.2.2 <- ~ I(scale(`Sucecional Forest`)^2) + I(scale(`Native forest`)^2) + I(scale(`Grassland`)^2) +I(scale(Cropland)^2)+ I(scale(DOM_score)^2) + scale(infrutescence) + I(scale(`Altitude (feet)`)^2) + scale(connected_trees)
btbw.det.formula <- ~ scale(precip)
out4.2.2 <- spPGOcc(occ.formula = btbw.occ.formula4.2.2,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 5000, 
                    batch.length = 25,
                    NNGP = FALSE, 
                    cov.model = "spherical", 
                    n.neighbors = 5, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out4.2.2)
ppc.out4.2.2 <- ppcOcc(out4.2.2, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out4.2.2)
waicOcc(out4.2.2)
names(out4.2.2)
str(out4.2.2$beta.samples)

MCMCplot(object= out4.2.2$beta.samples, guide_axis = TRUE ,rank = TRUE, ref_ovl=TRUE, mar = c(6.1, 5.1, 4.1, 4.1),labels= c("Intercept", expression("Successional Forest" ^ 2), expression("Native Forest" ^ 2),
                                                                                            expression("Grassland" ^ 2), expression("Cropland" ^ 2), expression("DOM score" ^ 2), "Infructescence", expression("Elevation" ^ 2), "Number surrounding palms"), ci=c(50, 93))
rMCMCplot(out4.2.2$alpha.samples, ref_ovl=TRUE, rank = TRUE, ci=c(50, 95), mar = c(6.1, 8.1, 4.1, 9.1),labels= c("Intercept", "Precipitation"))

#________________
####The best so so far! 2024 with Max temp

###Linear equation: logit(ψi​)=β0​+β1​⋅(Successional Forest)^2+β2​⋅(Native Forest)^2+β3​⋅(DOM Score)^2+β4​⋅(Cropland)^2+β5​⋅Infrutescence+β6​⋅Connected Trees

btbw.occ.formula4.2.2 <- ~ I(scale(`Sucecional Forest`)^2) + I(scale(`Native forest`)^2)+ I(scale(DOM_score)^2)+ I(scale(Cropland)^2) + scale(infrutescence) + scale(connected_trees)
btbw.det.formula <- ~ scale(tempmax)
out4.2.2 <- spPGOcc(occ.formula = btbw.occ.formula4.2.2,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 5000, 
                    batch.length = 25,
                    NNGP = TRUE, 
                    cov.model = "spherical", 
                    n.neighbors = 8, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out4.2.2)
ppc.out4.2.2 <- ppcOcc(out4.2.2, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out4.2.2)
waicOcc(out4.2.2)
names(out4.2.2)
str(out4.2.2$beta.samples)

MCMCplot(object= out4.2.2$beta.samples,  sz_ax_txt = 1.2, horiz = FALSE, rank = TRUE, guide_axis = TRUE, mar = c(6.1, 5.1, 4.1, 4.1),labels= c("Intercept", expression("Successional Forest" ^ 2), expression("Native Forest" ^ 2), expression("DOM score" ^ 2), expression("Cropland" ^ 2), "infrutescence", "Connected Trees"), ci=c(50, 92))

MCMCplot(out4.2.2$alpha.samples, sz_ax_txt = 1.2, horiz = FALSE, guide_axis = TRUE, ref_ovl=TRUE, ci=c(50, 92), mar = c(6.1, 8.1, 4.1, 9.1),labels= c("Intercept", "Temp. Max"))


###Best Model, Temp Max, Non spatial

####Linear equation: logit(ψi​)=β0​+β1​⋅(Successional Forest)^2+β2​⋅(Native Forest)^2+β3​⋅(DOM Score)^2+β4​⋅(Cropland)^2+β5​⋅Infrutescence+β6​⋅Connected Trees

btbw.occ.formula4b <- ~ I(scale(`Sucecional Forest`)^2) + I(scale(`Native forest`)^2)+ I(scale(DOM_score)^2)+ I(scale(Cropland)^2) + scale(infrutescence) + scale(connected_trees)
btbw.det.formula <- ~ scale(tempmax)
n.samples <- 5000

prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72))

inits.list <- list(alpha = 0, beta = 0,
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
n.report <- 1000

out4b <- PGOcc(occ.formula = btbw.occ.formula4b, 
               det.formula = btbw.det.formula, 
               data = data.list, 
               inits = inits.list,
               n.samples = n.samples,
               priors = prior.list,
               n.omp.threads = 1,
               verbose = TRUE,
               n.report = n.report, 
               n.burn = 1000, 
               n.thin = 1, 
               n.chains = 40)

summary(out4b)
ppc.out4b <- ppcOcc(out4b, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out4b)
waicOcc(out4b)
MCMCplot(object= out4b$beta.samples, guide_axis = TRUE , ref_ovl=TRUE, mar = c(6.1, 5.1, 4.1, 4.1),labels= c("Intercept", expression("Successional Forest" ^ 2), expression("Native Forest" ^ 2), expression("DOM score" ^ 2),
                                                                                                                expression("Cropland" ^ 2),"Artificial"), ci=c(50, 91))

###Null model
### Linear equation, Occupancy model (ψ): logit(ψi​)=β0​
### Linear equation, Detection model (p): logit(pij​)=α0​+α1​⋅Precipitation+α2​⋅Maximum Temperature

btbw.occ.formula4a <- ~ 1
btbw.det.formula <- ~ scale(precip)+ scale(tempmax)
out4a <- spPGOcc(occ.formula = btbw.occ.formula4a,
                 det.formula = btbw.det.formula,
                 data = data.list, 
                 n.batch = 4000, 
                 batch.length = 25,
                 NNGP = TRUE, 
                 cov.model = "exponential", 
                 n.neighbors = 15, 
                 n.burn = 2000, 
                 n.thin = 4, 
                 n.chains = 3, 
                 verbose = FALSE, 
                 k.fold = 2)

summary(out4a)
ppc.out4a <- ppcOcc(out4a, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out4a)
waicOcc(out4a)
names(out4a)
str(out4a$beta.samples)
MCMCplot(out4a$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), ci=c(50, 90))

MCMCplot(out4a$alpha.samples, ref_ovl=TRUE, ci=c(50, 90))


###
### Global model, Quadratic effect, spatial 

##Linear Equation: logit(ψi​)=β0​+β1​⋅(Native Forest)^2+β2​⋅(Successional Forest)^2+β3​⋅Grassland+β4​⋅(Artificial Land Cover)^2+β5​⋅(Cropland)^2+β6​⋅(DOM Score)^2+β7​⋅Infrutescence+β8​⋅(Altitude)^2+β9​⋅Stem Height+β10​⋅Connected Trees
btbw.occ.formula4.2.3 <- ~ I(scale(`Native forest`)^2) + I(scale(`Sucecional Forest`)^2) + scale(Grassland) + I(scale(Artificial)^2) + I(scale(Cropland)^2) +  I(scale(DOM_score)^2) + scale(infrutescence) + I(scale(`Altitude (feet)`)^2)+ scale(Stem_height) + scale(connected_trees)
btbw.det.formula <- ~ scale(precip)
out4.2.3 <- spPGOcc(occ.formula = btbw.occ.formula4.2.3,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 4000, 
                    batch.length = 15,
                    NNGP = TRUE, 
                    cov.model = "exponential", 
                    n.neighbors = 8, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out4.2.3)
ppc.out4.2.3 <- ppcOcc(out4.2.3, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out4.2.3)
waicOcc(out4.2.3)
str(out4.2.3$beta.samples)
MCMCplot(out4.2.3$beta.samples, guide_axis = TRUE, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1),labels= c("Intercept", "Native forest*", "Successional Forest*",
                                                                                                       "Grassland*",  "Artificial","Cropland*", "DOM score*", "Infructescence", "Elevation*", "Stem height"), ci=c(50, 95))
MCMCplot(out4.2.3$alpha.samples, ref_ovl=TRUE, labels= c("Intercept", "Precipitation"), ci=c(50, 95))

###
### Global model, Linear effect, spatial 

##Linear Equation: logit(ψi​)=β0​+β1​⋅(Native Forest)^2+β2​⋅(Successional Forest)^2+β3​⋅Grassland+β4​⋅(Artificial Land Cover)^2+β5​⋅(Cropland)^2+β6​⋅(DOM Score)^2+β7​⋅Infrutescence+β8​⋅(Altitude)^2+β9​⋅Stem Height+β10​⋅Connected Trees

btbw.occ.formula4.2.6 <- ~ scale(`Native forest`) + scale(`Sucecional Forest`) + scale(Grassland) + scale(Artificial) + scale(Cropland) +  scale(DOM_score) + scale(infrutescence) + scale(`Altitude (feet)`)+ scale(Stem_height) 
btbw.det.formula <- ~ scale(precip)
out4.2.6 <- spPGOcc(occ.formula = btbw.occ.formula4.2.6,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 4000, 
                    batch.length = 15,
                    NNGP = TRUE, 
                    cov.model = "exponential", 
                    n.neighbors = 8, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out4.2.6)
ppc.out4.2.6 <- ppcOcc(out4.2.6, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out4.2.6)
waicOcc(out4.2.6)

##No spatial, quadractic effect, best additive 


##Global model, Quadratic effect, No-spatial 
### logit(ψi​)=β0​+β1​⋅(Native Forest)^2+β2​⋅(Successional Forest)^2+β3​⋅Grassland+β4​⋅(Artificial Land Cover)^2+β5​⋅(Cropland)^2+β6​⋅(DOM Score)^2+β7​⋅Infrutescence+β8​⋅(Altitude)^2+β9​⋅Stem Height+β10​⋅Connected Trees

btbw.occ.formula4g <- ~ I(scale(`Native forest`)^2) + I(scale(`Sucecional Forest`)^2) + scale(Grassland) + I(scale(Artificial)^2) + I(scale(Cropland)^2) +  I(scale(DOM_score)^2) + scale(infrutescence) + I(scale(`Altitude (feet)`)^2)+ scale(Stem_height) + scale(connected_trees)
btbw.det.formula <- ~ scale(precip)
n.samples <- 5000

prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72))

inits.list <- list(alpha = 0, beta = 0,
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
n.report <- 1000

out4g <- PGOcc(occ.formula = btbw.occ.formula4g, 
               det.formula = btbw.det.formula, 
               data = data.list, 
               inits = inits.list,
               n.samples = n.samples,
               priors = prior.list,
               n.omp.threads = 1,
               verbose = TRUE,
               n.report = n.report, 
               n.burn = 1000, 
               n.thin = 1, 
               n.chains = 40)

summary(out4g)
ppc.out4g <- ppcOcc(out4g, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out4g)
waicOcc(out4g)

###
##Global model, Linear effect, No-spatial 

### Linear Equation: logit(ψi​)=β0​+β1​⋅(Native Forest)^2+β2​⋅(Successional Forest)^2+β3​⋅Grassland+β4​⋅(Artificial Land Cover)^2+β5​⋅(Cropland)^2+β6​⋅(DOM Score)^2+β7​⋅Infrutescence+β8​⋅(Altitude)^2+β9​⋅Stem Height+β10​⋅Connected Trees

btbw.occ.formula4h <- ~ scale(`Native forest`) + (`Sucecional Forest`) + scale(Grassland) + scale(Artificial) + scale(Cropland) +  scale(DOM_score) + scale(infrutescence) + scale(`Altitude (feet)`) + scale(Stem_height)
btbw.det.formula <- ~ scale(precip)
n.samples <- 5000

prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72))

inits.list <- list(alpha = 0, beta = 0,
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
n.report <- 1000

out4h <- PGOcc(occ.formula = btbw.occ.formula4h, 
               det.formula = btbw.det.formula, 
               data = data.list, 
               inits = inits.list,
               n.samples = n.samples,
               priors = prior.list,
               n.omp.threads = 1,
               verbose = TRUE,
               n.report = n.report, 
               n.burn = 1000, 
               n.thin = 1, 
               n.chains = 40)

summary(out4h)
ppc.out4h <- ppcOcc(out4h, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out4h)
waicOcc(out4h)


##No spatial, Linear effect 

##Linear Equation: logit(ψi​)=β0​+β1​⋅Infrutescence+β2​⋅Stem Height+β3​⋅Successional Forest+β4​⋅Grassland

btbw.occ.formula4c <- ~ scale(infrutescence) + scale(Stem_height) + scale(`Sucecional Forest`) + scale(Grassland) +
  btbw.det.formula <- ~ scale(precip)
n.samples <- 5000

prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72))

inits.list <- list(alpha = 0, beta = 0,
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
n.report <- 1000

out4c <- PGOcc(occ.formula = btbw.occ.formula4c, 
               det.formula = btbw.det.formula, 
               data = data.list, 
               inits = inits.list,
               n.samples = n.samples,
               priors = prior.list,
               n.omp.threads = 1,
               verbose = TRUE,
               n.report = n.report, 
               n.burn = 1000, 
               n.thin = 1, 
               n.chains = 40)

summary(out4c)
ppc.out4c <- ppcOcc(out4c, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out4c)
waicOcc(out4c)

##No spatial, Null model 

###Linear Equation: logit(ψi​)=β0​

btbw.occ.formula4d <- ~ 1
btbw.det.formula <- ~ scale(precip)
n.samples <- 5000

prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72))

inits.list <- list(alpha = 0, beta = 0,
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
n.report <- 1000

out4d <- PGOcc(occ.formula = btbw.occ.formula4d, 
               det.formula = btbw.det.formula, 
               data = data.list, 
               inits = inits.list,
               n.samples = n.samples,
               priors = prior.list,
               n.omp.threads = 1,
               verbose = TRUE,
               n.report = n.report, 
               n.burn = 1000, 
               n.thin = 1, 
               n.chains = 40)

summary(out4d)
ppc.out4d <- ppcOcc(out4d, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out4d)
waicOcc(out4d)

waicOcc(out4a)
waicOcc(out4b)
waicOcc(out4.2.3)
waicOcc(out4c)
waicOcc(out4.2.4)

names(out4b)

## 

#### Best model, Linear effect, spatial (reduced)

###Linear Equation: logit(ψi​)=β0​+β1​⋅Infrutescence+β2​⋅Stem Height+β3​⋅Successional Forest+β4​⋅Grassland+β5​⋅Altitude

btbw.occ.formula4.2.5 <- ~ scale(infrutescence) + scale(Stem_height) + scale(`Sucecional Forest`) + scale(Grassland) + scale(`Altitude (feet)`)
btbw.det.formula <- ~ scale(precip)
out4.2.5 <- spPGOcc(occ.formula = btbw.occ.formula4.2.5,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 4000, 
                    batch.length = 15,
                    NNGP = TRUE, 
                    cov.model = "exponential", 
                    n.neighbors = 8, 
                    n.burn = 2000, 
                    n.thin = 4, 
                    n.chains = 3, 
                    verbose = FALSE, 
                    k.fold = 2)

summary(out4.2.5)
ppc.out4.2.5 <- ppcOcc(out4.2.5, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out4.2.5)
waicOcc(out4.2.5)


##Another one (scaling)

###Linear Equation: logit(ψi​)=β0​+β1​⋅(Successional Forest)2+β2​⋅(Native Forest)2+β3​⋅(Grassland)2+β4​⋅DOM Score+β5​⋅Connected Trees

btbw.occ.formula5 <- ~ I(scale(`Sucecional Forest`)^2) + I(scale(`Native forest`)^2) + I(scale(Grassland)^2) + scale(DOM_score) + scale(connected_trees) 
btbw.det.formula <- ~ scale(precip) + precip
out5 <- spPGOcc(occ.formula = btbw.occ.formula5,
                det.formula = btbw.det.formula,
                data = data.list, 
                n.batch = 4000, 
                batch.length = 25,
                NNGP = FALSE, 
                cov.model = "exponential", 
                search.type = 'cb', 
                n.neighbors = 15, 
                n.burn = 2000, 
                n.thin = 4, 
                n.chains = 3, 
                verbose = FALSE, 
                k.fold = 2)

summary(out5)
ppc.out5 <- ppcOcc(out5, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out5)
waicOcc(out5)
names(out5)
str(out5$beta.samples)
MCMCplot(out5$beta.samples, ref_ovl=TRUE, ci=c(50, 90))

MCMCplot(out5$alpha.samples, ref_ovl=TRUE, ci=c(50, 90))

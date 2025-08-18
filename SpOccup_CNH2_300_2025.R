
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

###Linear model: logit(ψi​)=β0​+β1​⋅Sucecional Forest^2​+β2​⋅Grassland^2​+β3​⋅infrutescence

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

###Linear Model: logit(ψi​)=β0​+β1​⋅Sucecional Forest^2​+β2​⋅Grassland^2​+β3​⋅Cropland^2​+β4​⋅Artificial^2​
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
####Linear Model: logit(ψi​)=β0​+β1​⋅Sucecional Forest^2​+β2​⋅Grassland^2​+β3​⋅Cropland^2​+β4​⋅Artificial^2​+β5​⋅infrutescence​
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


#The best so so far_Non spatial with precipitation 
btbw.occ.formula4b <- ~ I(scale(`Sucecional Forest`)^2) + I(scale(`Native forest`)^2) + I(scale(`Grassland`)^2) +I(scale(Cropland)^2)+ I(scale(DOM_score)^2) + scale(infrutescence) + I(scale(`Altitude (feet)`)^2) + scale(connected_trees)
btbw.det.formula <- ~ scale(precip)
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

MCMCplot(out4b$beta.samples, guide_axis = TRUE ,rank = TRUE, ref_ovl=TRUE, mar = c(6.1, 5.1, 4.1, 4.1),labels= c("Intercept", expression("Successional Forest" ^ 2), expression("Native Forest" ^ 2),
                                                                                                                 expression("Grassland" ^ 2), expression("Cropland" ^ 2), expression("DOM score" ^ 2), "Infructescence", expression("Elevation" ^ 2), "Number surrounding palms"), ci=c(50, 93))
###




## Best model, Quadratic effect, spatial 

btbw.occ.formula4.2.3 <- ~ I(scale(`Sucecional.Forest`)^2) * I(scale(`Native.Forest`)^2)+ I(scale(Grassland)^2) + I(scale(Artificial)^2) + I(scale(Cropland)^2) +  I(scale(DOM_score)^2) + scale(infrutescence) + I(scale(`Altitude..feet.`)^2)+ scale(Stem_height) 
btbw.det.formula <- ~ scale(precip) 
out4.2.3 <- spPGOcc(occ.formula = btbw.occ.formula4.2.3,
                    det.formula = btbw.det.formula,
                    data = data.list, 
                    n.batch = 5000, 
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
names(out4.2.3)
str(out4.2.3$beta.samples)
MCMCplot(out4.2.3$beta.samples, guide_lines = TRUE, guide_axis = TRUE, ref_ovl=TRUE, 
         mar = c(5.1, 4.1, 4.1, 2.1),labels= c("Intercept", "Successional Forest*", "Native Forest",
                                               "Grassland*", "Artifical*", "Cropland*", "DOM score*", "Infructescence", "Elevation*", "Stem height"), ci=c(75, 95), rank = TRUE)
MCMCplot(out4.2.3$alpha.samples, ref_ovl=TRUE, ci=c(85, 95))
fitted(ppc.out4.2.3)
fitted(out4.2.3)

#### Best model, Linear effect, spatial 

btbw.occ.formula4.2.4 <- ~  scale(`Sucecional Forest`) + scale(Grassland) + scale(Artificial) + scale(Cropland) +  scale(DOM_score) + scale(infrutescence) + scale(`Altitude (feet)`)+ scale(Stem_height) 
btbw.det.formula <- ~ scale(precip)
out4.2.4 <- spPGOcc(occ.formula = btbw.occ.formula4.2.4,
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

summary(out4.2.4)
ppc.out4.2.4 <- ppcOcc(out4.2.4, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out4.2.4)
waicOcc(out4.2.4)
names(out4.2.4)
str(out4.2.4$beta.samples)
MCMCplot(out4.2.4$beta.samples, guide_lines = TRUE, guide_axis = TRUE, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1),labels= c("Intercept", "Successional Forest*",
                                                                                                                           "Grassland*", "Artifical*", "Cropland*", "DOM score*", "Infructescence", "Elevation*", "Stem height"), ci=c(50, 95), rank = TRUE)
MCMCplot(out4.2.4$alpha.samples, ref_ovl=TRUE, ci=c(50, 95))
fitted(ppc.out4.2.4)
fitted(out4.2.4)


###Null model

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


##No spatial, Linear effect, best additive 

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



# Predictions -----------------------------------------------------------
# Predict occupancy along a gradient of Successional forest cover.  
# Create a set of values across the range of observed Successional forest values
Sucecional.pred.vals <- seq(min(data.list$occ.covs$`Sucecional Forest`), 
                            max(data.list$occ.covs$`Sucecional Forest`), 
                            length.out = 100)

# Scale predicted values by mean and standard deviation used to fit the model

Sucecional.pred.vals.scale <- ((Sucecional.pred.vals - mean(data.list$occ.covs$`Sucecional Forest`)) / sd(data.list$occ.covs$`Sucecional Forest`))

# Predict occupancy across forest values at mean values of all other variables
pred.df <- as.matrix(data.frame(intercept = 1, Sucecional = Sucecional.pred.vals.scale, 
                                Grassland = 0, Cropland = 0, 
                                Artificial = 0, DOM_score = 0, infrutescence=0,`Altitude (feet)`= 0, Stem_height=0))

out.pred <- predict(out4b, pred.df)
str(out.pred)
psi.0.quants <- apply(out.pred$psi.0.samples, 2, quantile, 
                      prob = c(0.05, 0.5, 0.95))
psi.plot.dat <- data.frame(psi.med = psi.0.quants[2, ], 
                           psi.low = psi.0.quants[1, ], 
                           psi.high = psi.0.quants[3, ], 
                           Sucecional = Sucecional.pred.vals)
ggplot(psi.plot.dat, aes(x = Sucecional, y = psi.med)) + 
  geom_ribbon(aes(ymin = psi.low, ymax = psi.high), alpha= 0.4, fill = 'olivedrab4') +
  geom_line() + 
  theme_bw() + 
  scale_y_continuous(limits = c(0, 1)) + 
  labs(x = 'Successional (% cover)', y = 'Occupancy Probability')

# Predict occupancy along a gradient of Elevation cover.  
# Create a set of values across the range of observed Successional forest values
Elevation.pred.vals <- seq(min(data.list$occ.covs$`Altitude (feet)`), 
                           max(data.list$occ.covs$`Altitude (feet)`), 
                           length.out = 100)

# Scale predicted values by mean and standard deviation used to fit the model

Elevation.pred.vals.scale <- ((Elevation.pred.vals - mean(data.list$occ.covs$`Altitude (feet)`)) / sd(data.list$occ.covs$`Altitude (feet)`))

# Predict occupancy across forest values at mean values of all other variables
pred.df1 <- as.matrix(data.frame(intercept = 1, Elevation = scale(Elevation.pred.vals.scale), 
                                 Grassland = 0, Cropland = 0, 
                                 Artificial = 0, DOM_score = 0, infrutescence=0,`Sucecional Forest`= 0, Stem_height=0))

out.pred1 <- predict(out4b, pred.df1)
str(out.pred1)
psi.0.quants <- apply(out.pred1$psi.0.samples, 2, quantile, 
                      prob = c(0.1, 0.5, 0.9))
psi.plot.dat <- data.frame(psi.med = psi.0.quants[2, ], 
                           psi.low = psi.0.quants[1, ], 
                           psi.high = psi.0.quants[3, ], 
                           Elevation = Elevation.pred.vals)
ggplot(psi.plot.dat, aes(x = Elevation, y = psi.med)) + 
  geom_ribbon(aes(ymin = psi.low, ymax = psi.high), alpha= 0.4, fill = 'cyan4') +
  geom_line() + 
  theme_bw() + 
  scale_y_continuous(limits = c(0, 1)) + 
  labs(x = 'Elevation', y = 'Occupancy Probability')


# Predict occupancy along a gradient of Grassland cover.  
# Create a set of values across the range of observed Successional forest values
Grassland.pred.vals <- seq(min(data.list$occ.covs$Grassland), 
                           max(data.list$occ.covs$Grassland), 
                           length.out = 100)

# Scale predicted values by mean and standard deviation used to fit the model

Grassland.pred.vals.scale <- ((Grassland.pred.vals - mean(data.list$occ.covs$Grassland)) / sd(data.list$occ.covs$Grassland))

# Predict occupancy across forest values at mean values of all other variables
pred.df2 <- as.matrix(data.frame(intercept = 1, Grassland = Grassland.pred.vals.scale, 
                                 `Altitude (feet)`=0, Cropland = 0, 
                                 Artificial = 0, DOM_score = 0, infrutescence=0,`Sucecional Forest`= 0, Stem_height=0))

out.pred2 <- predict(out4b, pred.df2)
str(out.pred2)
psi.0.quants <- apply(out.pred2$psi.0.samples, 2, quantile, 
                      prob = c(0.1, 0.5, 0.9))
psi.plot.dat <- data.frame(psi.med = psi.0.quants[2, ], 
                           psi.low = psi.0.quants[1, ], 
                           psi.high = psi.0.quants[3, ], 
                           Grassland = Grassland.pred.vals)
ggplot(psi.plot.dat, aes(x =Grassland, y = psi.med)) + 
  geom_ribbon(aes(ymin = psi.low, ymax = psi.high), alpha= 0.4, fill = 'tomato3') +
  geom_line() + 
  theme_bw() + 
  scale_y_continuous(limits = c(0, 1)) + 
  labs(x = 'Grassland', y = 'Occupancy Probability')


###################################


##More models


##Another one (scaling)

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


###Just land use

btbw.occ.formula6 <- ~ scale(`Sucecional Forest`) + scale(Grassland) + scale(`Native forest`) + scale (Cropland) + scale (Artificial)
btbw.det.formula <- ~ scale(precip)
out6 <- spPGOcc(occ.formula = btbw.occ.formula6,
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

summary(out6)
ppc.out6 <- ppcOcc(out6, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out6)
waicOcc(out6)
names(out6)
str(out6$beta.samples)
MCMCplot(out6$beta.samples, ref_ovl=TRUE, ci=c(50, 89))

MCMCplot(out6$alpha.samples, ref_ovl=TRUE, ci=c(50, 89))


#stPGOcc

# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = 0, var = 2.72), 
                   sigma.sq.ig = c(2, 2), 
                   phi.unif = c(3 / 1, 3 / 0.1), 
                   rho.unif = c(-1, 1),
                   sigma.sq.t.ig = c(2, 1))

# Initial values
z.init <- apply(data.list$y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(beta = 0, alpha = 0, z = z.init, phi = 3 / .5, sigma.sq = 2, rho = 0, sigma.sq.t = 0.5)


# Detection-nondetection data
y <- data.list$y
# Occupancy covariates
X <- data.list$occ.covs
# Detection covarites
X.p <- data.list$det.covs
# Spatial coordinates
coords <- data.list$coords

# Package all data into a list
occ.covs <- X[, -1, drop = FALSE]
colnames(occ.covs) <- c('occ.cov')
det.covs <- list(det.cov.1 = X.p[, , 2])
data.list <- list(y = y, 
                  occ.covs = occ.covs, 
                  det.covs = det.covs, 
                  coords = coords)

# Tuning
tuning.list <- list(phi = 1, rho = 1)
# Number of batches
n.batch <- 10
# Batch length
batch.length <- 25
n.iter <- n.batch * batch.length

trend<-TRUE

occ.covs <- list(int = data.list$X[, , 1],
                 trend = data.list$X[, , 2],
                 occ.cov.1 = data.list$X[, , 3])

out7 <- stPGOcc(btbw.occ.formula4,
                btbw.det.formula, 
                data = data.list, 
                inits = inits.list, 
                n.batch = n.batch, 
                batch.length = batch.length, 
                priors = prior.list,
                cov.model = "exponential", 
                tuning = tuning.list, 
                NNGP = TRUE, 
                ar1 = TRUE,
                n.neighbors = 5, 
                search.type = 'cb', 
                n.report = 10, 
                n.burn = 50, 
                n.chains = 1)


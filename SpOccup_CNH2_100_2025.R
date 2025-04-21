####Code for landscape compositions withing buffers of 100 m radius

library(spOccupancy)
library(datasets)
library(base)
library(MCMCvis)


# 1. Defining matrices Preparing the data -----------------------------------------------------------

y<-as.matrix(occurrance_300c)
y<-data.matrix(occurrance_300c)
str(y)
occ.covs<- as.matrix(cov_occ_100)
View(occ.covs)
#occ.covs<- as.data.frame(cov_occ_100)
str(occ.covs)
precip<- as.matrix(precip_300_c)
tempmax<- as.matrix(temp_300_b)
tempmean<- as.matrix(MeanTemp)
tempmin<- as.matrix(MinTemp)

det.covs<-list(precip = precip, 
               tempmax = tempmax,
               tempmean= tempmean,
               tempmin=tempmin)
str(det.covs)

coords<- as.matrix(cord_UTM_300)
str(coords)

data.list_2 <- list(y = y, 
                  occ.covs = occ.covs, 
                  det.covs = det.covs, 
                  coords = coords)

plot(data.list$coords, pch=19)

str(data.list_2)

# 2. Models -----------------------------------------------------------

####The best so so far!
btbw.occ.formula4.2.2 <- ~ I(scale(`Sucecional Forest`)^2) + I(scale(`Native forest`)^2) + I(scale(`Grassland`)^2) +I(scale(Cropland)^2)+ I(scale(DOM_score)^2) + scale(infrutescence) + I(scale(`Altitude (feet)`)^2) + scale(connected_trees)
btbw.det.formula <- ~ scale(precip)
out4.2.2 <- spPGOcc(occ.formula = btbw.occ.formula4.2.2,
                    det.formula = btbw.det.formula,
                    data = data.list_2, 
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

summary(out4.2.2)
ppc.out4.2.2 <- ppcOcc(out4.2.2, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out4.2.2)
waicOcc(out4.2.2)

names(out4.2.2)
str(out4.2.2$beta.samples)
MCMCplot(out4.2.2$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", expression("Successional Forest" ^ 2), expression("Native Forest" ^ 2),
                                                                                     expression("Grassland" ^ 2), expression("Cropland" ^ 2), expression("DOM score" ^ 2), "Infructescence", expression("Elevation" ^ 2), "Number surrounding palms"), ci=c(50, 90))

MCMCplot(out4.2.2$alpha.samples, ref_ovl=TRUE, ci=c(50, 89))



####Single variables
btbw.occ.formula1<- ~ scale(`Sucecional Forest`)
btbw.det.formula <- ~ scale(tempmax) 
out1 <- spPGOcc(occ.formula = btbw.occ.formula1,
                    det.formula = btbw.det.formula,
                    data = data.list_2, 
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

summary(out1)
ppc.out1 <- ppcOcc(out1, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out1)
waicOcc(out1)
names

str(out1$beta.samples)
MCMCplot(out1$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "successional forest"),ci=c(50, 95))

#Native Forest

btbw.occ.formula2<- ~ I(scale(`Native forest`)^2)
btbw.det.formula <- ~ scale(tempmax)
out2 <- spPGOcc(occ.formula = btbw.occ.formula2,
                det.formula = btbw.det.formula,
                data = data.list_2, 
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

summary(out2)
ppc.out2 <- ppcOcc(out2, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out2)
waicOcc(out2)
names(out2)
str(out2$beta.samples)
MCMCplot(out2$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "Native forest"),ci=c(50, 95))


###Grassland

btbw.occ.formula3<- ~ I(scale(Grassland)^2)
btbw.det.formula <- ~ scale(tempmax)
out3 <- spPGOcc(occ.formula = btbw.occ.formula3,
                det.formula = btbw.det.formula,
                data = data.list_2, 
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

summary(out3)
ppc.out3 <- ppcOcc(out3, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out3)
waicOcc(out3)
names(out3)
str(out3$beta.samples)
MCMCplot(out3$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "Grassland"),ci=c(50, 95))


###Cropland

btbw.occ.formula4<- ~ I(scale(Cropland)^2)
btbw.det.formula <- ~ scale(precip)
out4 <- spPGOcc(occ.formula = btbw.occ.formula4,
                det.formula = btbw.det.formula,
                data = data.list_2, 
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

summary(out4)
ppc.out4 <- ppcOcc(out4, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out4)
waicOcc(out4)
names(out4)
str(out4$beta.samples)
MCMCplot(out4$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "Cropland"),ci=c(50, 95))

### Elevation

btbw.occ.formula5<- ~ scale(`Altitude (feet)`)
btbw.det.formula <- ~ scale(precip)
out5 <- spPGOcc(occ.formula = btbw.occ.formula5,
                det.formula = btbw.det.formula,
                data = data.list_2, 
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

summary(out5)
ppc.out5 <- ppcOcc(out5, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out5)
waicOcc(out5)
names(out5)
str(out5$beta.samples)
MCMCplot(out5$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "Elevation"),ci=c(50, 90))

###

btbw.occ.formula6<- ~ I(scale(Stem_height)^2)
btbw.det.formula <- ~ scale(tempmax)
out6 <- spPGOcc(occ.formula = btbw.occ.formula6,
                det.formula = btbw.det.formula,
                data = data.list_2, 
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
MCMCplot(out6$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "Stem height"),ci=c(50, 89))

###DOM score

btbw.occ.formula7<- ~ I(scale(DOM_score)^2)
btbw.det.formula <- ~ scale(tempmax)
out7 <- spPGOcc(occ.formula = btbw.occ.formula7,
                det.formula = btbw.det.formula,
                data = data.list_2, 
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

summary(out7)
ppc.out7 <- ppcOcc(out7, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out7)
waicOcc(out7)
names(out7)
str(out7$beta.samples)
MCMCplot(out7$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "DOM score"),ci=c(50, 89))

### Infructescence 

btbw.occ.formula8<- ~ scale(infrutescence)
btbw.det.formula <- ~ scale(tempmax)
out8 <- spPGOcc(occ.formula = btbw.occ.formula8,
                det.formula = btbw.det.formula,
                data = data.list_2, 
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

summary(out8)
ppc.out8 <- ppcOcc(out8, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out8)
waicOcc(out8)
names(out8)
str(out8$beta.samples)
MCMCplot(out8$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "infrutescence"),ci=c(50, 90))

###

btbw.occ.formula9<- ~ scale(connected_trees)
btbw.det.formula <- ~ scale(tempmax)
out9 <- spPGOcc(occ.formula = btbw.occ.formula9,
                det.formula = btbw.det.formula,
                data = data.list_2, 
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

summary(out9)
ppc.out9 <- ppcOcc(out9, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out9)
waicOcc(out9)
names(out9)
str(out9$beta.samples)
MCMCplot(out9$beta.samples, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1), labels= c("Intercept", "Connected Trees"),ci=c(50, 90))


### Null model

btbw.occ.formula0<- ~ 1
btbw.det.formula <- ~ scale(tempmean)
out0 <- spPGOcc(occ.formula = btbw.occ.formula0,
                det.formula = btbw.det.formula,
                data = data.list_2, 
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

summary(out0)
ppc.out0 <- ppcOcc(out0, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out0)
waicOcc(out0)
names(out0)
str(out0$beta.samples)
MCMCplot(out0$alpha.samples, ref_ovl=TRUE, ci=c(50, 89))

# Model comparisons __________ 

waicOcc(out1)
waicOcc(out2)
waicOcc(out3)
waicOcc(out4)
waicOcc(out5)
waicOcc(out6)
waicOcc(out7)
waicOcc(out8)
waicOcc(out9)
waicOcc(out0)

### Additive

#Best Additive model
btbw.occ.formula10 <- ~ scale(`Sucecional Forest`) + scale(Cropland) + I(scale(DOM_score)^2)+ scale(infrutescence)
btbw.det.formula <- ~ scale(tempmax)
out.10 <- spPGOcc(occ.formula = btbw.occ.formula10,
                    det.formula = btbw.det.formula,
                    data = data.list_2, 
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

summary(out.10)
ppc.out.10<- ppcOcc(out.10, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out.10)
waicOcc(out.10)
names(out.10)
str(out.10$beta.samples)

MCMCplot(out.10$beta.samples, horiz = FALSE, sz_ax_txt = 1.2, guide_axis = TRUE, mar = c(5.1, 4.1, 4.1, 2.1),labels= c("Intercept", "Successional forest", "Cropland", expression("DOM score" ^ 2), "infructescence"), ci=c(50, 92), rank = TRUE)


MCMCplot(out.10$alpha.samples, guide_axis = TRUE, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1),labels= c("Intercept", "Daily Maximum Temperature "), ci=c(50, 95), rank = TRUE)

##Best model_nonspatial
btbw.occ.formula4b <-  ~ scale(`Sucecional Forest`) + scale(Cropland) + I(scale(DOM_score)^2) + scale(infrutescence)
btbw.det.formula <- ~ scale(tempmax)
n.samples <- 5000

prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72))

inits.list <- list(alpha = 0, beta = 0,
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
n.report <- 1000

out4b <- PGOcc(occ.formula = btbw.occ.formula4b, 
               det.formula = btbw.det.formula, 
               data = data.list_2, 
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

MCMCplot(out.10$beta.samples, guide_lines = TRUE, guide_axis = TRUE, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1),labels= c("Intercept", "Successional", "Cropland"), ci=c(50, 95), rank = TRUE)
MCMCplot(out.10$alpha.samples, ref_ovl=TRUE, ci=c(50, 95))
fitted(ppc.out.10)
fitted(out.10)

###

btbw.occ.formula11 <- ~  I(scale(Stem_height)^2) + I(scale(DOM_score)^2) + scale(infrutescence) + I(scale(`Sucecional Forest`)^2) 
btbw.det.formula <- ~ scale(precip)
out.11 <- spPGOcc(occ.formula = btbw.occ.formula11,
                  det.formula = btbw.det.formula,
                  data = data.list_2, 
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

summary(out.11)
ppc.out.11<- ppcOcc(out.11, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out.11)
waicOcc(out.11)
names(out.11)
str(out.11$beta.samples)
MCMCplot(out.11$beta.samples, guide_lines = TRUE, guide_axis = TRUE, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1),labels= c("Intercept", "Stem height", "DOM score", "Infrutescence", "sucessional forest"), ci=c(50, 95), rank = TRUE)
MCMCplot(out.11$alpha.samples, ref_ovl=TRUE, ci=c(50, 95))
fitted(ppc.out.11)
fitted(out.11)

###

btbw.occ.formula12 <- ~  I(scale(Stem_height)^2) + I(scale(DOM_score)^2) + scale(infrutescence) + I(scale(`Sucecional Forest`)^2) + I(scale(`Altitude (feet)`)^2)
btbw.det.formula <- ~ scale(precip)
out.12 <- spPGOcc(occ.formula = btbw.occ.formula12,
                  det.formula = btbw.det.formula,
                  data = data.list_2, 
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

summary(out.12)
ppc.out.12<- ppcOcc(out.12, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out.12)
waicOcc(out.12)
names(out.12)
str(out.12$beta.samples)
MCMCplot(out.12$beta.samples, guide_lines = TRUE, guide_axis = TRUE, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1),labels= c("Intercept", "Stem height", "DOM score", "Infrutescence", "sucessional forest", "Elevation"), ci=c(50, 95), rank = TRUE)
MCMCplot(out.12$alpha.samples, ref_ovl=TRUE, ci=c(50, 95))
fitted(ppc.out.12)
fitted(out.12)

###

btbw.occ.formula13 <- ~  I(scale(Stem_height)^2) + I(scale(DOM_score)^2) + scale(infrutescence) + I(scale(`Altitude (feet)`)^2)
btbw.det.formula <- ~ scale(precip)
out.13 <- spPGOcc(occ.formula = btbw.occ.formula13,
                  det.formula = btbw.det.formula,
                  data = data.list_2, 
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

summary(out.13)
ppc.out.13<- ppcOcc(out.13, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out.13)
waicOcc(out.13)
names(out.13)
str(out.13$beta.samples)
MCMCplot(out.13$beta.samples, guide_lines = TRUE, guide_axis = TRUE, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1),labels= c("Intercept", "Stem height", "DOM score", "Infrutescence", "Elevation"), ci=c(50, 95), rank = TRUE)
MCMCplot(out.13$alpha.samples, ref_ovl=TRUE, ci=c(50, 95))
fitted(ppc.out.13)
fitted(out.13)

###
btbw.occ.formula14 <- ~  I(scale(Stem_height)^2) + I(scale(DOM_score)^2) + scale(infrutescence) + I(scale(`Sucecional Forest`)^2) + I(scale(Cropland)^2)
btbw.det.formula <- ~ scale(precip)
out.14 <- spPGOcc(occ.formula = btbw.occ.formula14,
                  det.formula = btbw.det.formula,
                  data = data.list_2, 
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

summary(out.14)
ppc.out.14<- ppcOcc(out.14, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out.14)
waicOcc(out.14)
names(out.14)
str(out.14$beta.samples)
MCMCplot(out.14$beta.samples, guide_lines = TRUE, guide_axis = TRUE, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1),labels= c("Intercept", "Stem height", "DOM score", "Infrutescence", "sucessional forest", "Grassland"), ci=c(50, 95), rank = TRUE)
MCMCplot(out.14$alpha.samples, ref_ovl=TRUE, ci=c(50, 95))
fitted(ppc.out.14)
fitted(out.14)


###
btbw.occ.formula15 <- ~  I(scale(Stem_height)^2) + I(scale(DOM_score)^2) + scale(infrutescence) + I(scale(`Sucecional Forest`)^2) + scale(connected_trees)
btbw.det.formula <- ~ scale(precip)
out.15 <- spPGOcc(occ.formula = btbw.occ.formula15,
                  det.formula = btbw.det.formula,
                  data = data.list_2, 
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

summary(out.15)
ppc.out.15<- ppcOcc(out.15, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out.15)
waicOcc(out.15)
names(out.15)
str(out.15$beta.samples)
MCMCplot(out.15$beta.samples, guide_lines = TRUE, guide_axis = TRUE, ref_ovl=TRUE, mar = c(5.1, 4.1, 4.1, 2.1),labels= c("Intercept", "Stem height", "DOM score", "Infrutescence", "sucessional forest", "Connected trees"), ci=c(50, 95), rank = TRUE)
MCMCplot(out.15$alpha.samples, ref_ovl=TRUE, ci=c(50, 95))
fitted(ppc.out.15)
fitted(out.15)


waicOcc(out1)
waicOcc(out2)
waicOcc(out3)
waicOcc(out4)
waicOcc(out5)
waicOcc(out6)
waicOcc(out7)
waicOcc(out8)
waicOcc(out9)
waicOcc(out0)
waicOcc(out.10)

###The best
waicOcc(out.11)
waicOcc(out.12)
waicOcc(out.13)
waicOcc(out.14)
waicOcc(out.15)

#####

btbw.occ.formula4b <- ~ I(scale(`Native forest`)^2) +I(scale(`Sucecional Forest`)^2) + I(scale(Grassland)^2) + scale(Cropland) +  I(scale(DOM_score)^2) + scale(infrutescence) + I(scale(`Altitude (feet)`)^2)+ scale(Stem_height)
btbw.det.formula <- ~ scale(precip)
n.samples <- 5000

prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72))

inits.list <- list(alpha = 0, beta = 0,
                   z = apply(data.list$y, 1, max, na.rm = TRUE))
n.report <- 1000

out4b <- PGOcc(occ.formula = btbw.occ.formula4b, 
               det.formula = btbw.det.formula, 
               data = data.list_2, 
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


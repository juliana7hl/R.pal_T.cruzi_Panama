library(mgcv)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(tidymv)
library(lme4)
library(glmmTMB)
library("ggpubr")
library(pROC)
library(MuMIn)
library(lmtest)
library(MASS)
library(effects)

#####
Datasetm<-data.frame(Dataset_300mtemp)

#GLMs
glm_null <- glm(T_cruzi ~ 1, family = binomial, data = Datasetm)
summary(glm_null)

##Step_function
#Excluded Native Forest since showed High corraltion (spearman) with Grassland, Artificial surfaces and Successional Forest
glm_model <- glm(T_cruzi ~ Artificial_surfaces + factor(Sampling.Season) + Successional_forest + Grassland*connected_trees + Daily_Precipitacion + Mean_Temp + X*Y + Native_forest*connected_trees + factor(Stage) + infrutescence + DOM_score, family = binomial, data = Datasetm)
summary(glm_model)
AIC(glm_model)
# Get the null deviance and residual deviance
null_deviance <- GAM_Step$null.deviance
residual_deviance <- GAM_Step$deviance
# Calculate the deviance explained
deviance_explained <- 1 - (residual_deviance / null_deviance)
# Print the result
deviance_explained
#Step model
step_model <- step(glm_model)

#After Step selection
glm_model_1 <- glm(T_cruzi ~  Native_forest + Successional_forest + Grassland + Daily_Precipitacion  + X*Y + Stage + infrutescence + connected_trees, family = binomial, data = Datasetm)
summary(glm_model_1)
AIC(glm_model_1)
cor(predict(glm_model_1, newdata=Datasetm),Datasetm$T_cruzi)^2



#GAMs __________________________________
##Null model
Nullmodel <- gam(T_cruzi ~ 1, family = binomial, data = Datasetm)
summary(Nullmodel)
#AIC=182.8356
AIC(Nullmodel)

#GAMs, landscape composition
#Global model (Deviance explained = 49.4%), AIC:149.8632
GAM_Global <- gam(T_cruzi ~ s(Artificial_surfaces, bs = "cc")+factor(Sampling.Season) +s(Successional_forest, bs = "cs") +s(Native_forest, bs = "cr") +s(Grassland, bs = "cc") + s(Daily_Precipitacion, bs = "cc") + te(X, Y) + factor(Stage) + infrutescence + DOM_score + connected_trees + Altitude..feet., family = binomial, data = Datasetm)
summary(GAM_Global)
AIC(GAM_Global)
plot(GAM_Global, pages=1, residuals=TRUE, cex = 0.7, pch = 20, shade=TRUE, scheme = 2, shade.col="khaki1")

#Only significant predictor from GLobal Model
GAM_In <- gam(T_cruzi ~ Sampling.Season + infrutescence  + connected_trees + te(X, Y), family = binomial, data = Datasetm)
summary(GAM_In)

###################################
#Additive models with Cubic splines
###################################

#Defining factors
Datasetm$Palm_ID <- as.factor(Datasetm$Palm_ID)
Datasetm$Stage <- as.factor(Datasetm$Stage)

#AIC=143.4866, Deviance explained = 45.2%
GAM_Step <- gam(T_cruzi ~  s(Palm_ID, bs = "re") + s(Daily_Precipitacion, bs = "cc") + s(Grassland, bs = "cr") +s(Native_forest, bs = "cs") + connected_trees + ti(X, Y) + factor(Stage), family = binomial, data = Datasetm)
summary(GAM_Step)
AIC(GAM_Step)
plot(GAM_Step, residuals=TRUE, lwd= 2, cex = 0.6, pch = 19, shade=TRUE, scheme = 2, shade.col="khaki1")
cor(predict(GAM_Step, newdata=Datasetm),Datasetm$T_cruzi)^2

#AIC= 151.4037
GAM_Step2 <- gam(T_cruzi ~s(Daily_Precipitacion, bs = "cc")+s(Native_forest, bs = "cs") + connected_trees + ti(X, Y) + factor(Stage)  + infrutescence + connected_trees , family = binomial, data = Datasetm)
summary(GAM_Step2)
AIC(GAM_Step2)
plot(GAM_Step2, page=1, residuals=TRUE, cex = 0.7, pch = 20, shade=TRUE, scheme = 2, shade.col="khaki1")

vis.gam(GAM_Step2, view = c("Grassland", "Native_forest"), plot.type = "contour", color = "terrain",
        main = "Partial Dependence", xlab = "Grassland", ylab = "Native_forest", 
        n.grid = 100)  # n.grid sets the resolution of the plot

##GLM alternative
GLM_Step2 <- glm(T_cruzi ~ s(Palm_ID, bs = "re") +Daily_Precipitacion + Grassland + Native_forest + Grassland*Native_forest + Grassland*Successional_forest + Native_forest*Successional_forest + Native_forest*Artificial_surfaces +connected_trees + X*Y + factor(Stage)  + infrutescence + connected_trees , family = binomial, data = Datasetm)
summary(GLM_Step2)
AIC(GLM_Step2)
plot(GLM_Step2)

#AIC=147.8935, Deviance explained = 40%
GAM_additive <- gam(T_cruzi ~ s(Artificial_surfaces, bs = "cc") +s(Daily_Precipitacion, bs = "cc") + s(Successional_forest, bs = "cr")+ s(Grassland, bs = "cr") + connected_trees + ti(X, Y) + factor(Stage) + infrutescence + DOM_score, family = binomial, data = Datasetm)
summary(GAM_additive)
AIC(GAM_additive)
summary(GAM_additive)$p.coeff  # For parametric coefficients
summary(GAM_additive)$s.table
coef(GAM_additive)
gam.check(GAM_additive)
plot(GAM_additive, pages=1, residuals=TRUE, cex = 0.7, pch = 20, shade=TRUE, scheme = 2, shade.col="lightgrey")

GLM_additive <- glm(T_cruzi ~ Artificial_surfaces + Daily_Precipitacion + Successional_forest + Grassland + connected_trees + Y*X + factor(Stage) + infrutescence + DOM_score, family = binomial, data = Datasetm)
AIC(GLM_additive)
summary(GLM_additive)
plot(GLM_additive)

#Deviance explained (31%) AIC(155.1057)
GAM_additive_2 <- gam(T_cruzi ~ s(Successional_forest, bs = "cs") + s(Daily_Precipitacion, bs = "cc") + ti(X, Y) + factor(Stage) + infrutescence + DOM_score, family = binomial, data = Datasetm)
summary(GAM_additive_2)
AIC(GAM_additive_2)
plot(GAM_additive_2, pages=1, residuals=TRUE, cex = 0.7, pch = 20, shade=TRUE, scheme = 2, shade.col="khaki1")

#Deviance explained (37.8%), AIC 146.9763
GAM_additive_3 <- gam(T_cruzi ~ s(Artificial_surfaces, bs = "cc") +s(Successional_forest, bs = "cr") +s(Native_forest, bs = "cr") + s(Daily_Precipitacion, bs = "cc") + ti(X, Y) + factor(Stage) + infrutescence + DOM_score + connected_trees, family = binomial, data = Datasetm)
summary(GAM_additive_3)
AIC(GAM_additive_3)
logLik(GAM_additive_3)
plot(GAM_additive_3, pages=1, residuals=TRUE, cex = 0.7, pch = 20, shade=TRUE, scheme = 2, shade.col="khaki1")

#Model adding Mean Temp, (42.2%) AIC 145.8177
GAM_Global_Clim <- gam(T_cruzi ~ s(Grassland, bs = "cr")+s(Native_forest, bs = "cs")  + s(Daily_Precipitacion, bs = "cc")+ ti(X, Y) + Stage + connected_trees, family = binomial, data = Datasetm)
summary(GAM_Global_Clim)
AIC(GAM_Global_Clim)
plot(GAM_Global_Clim, residuals=TRUE, cex = 0.7, pch = 20, shade=TRUE, scheme = 2, shade.col="khaki1")

##GAMs without cubic splines
#Deviance explained = 36.7%, AIC=151.323
GAM_wosplines <- gam(T_cruzi ~ s(Artificial_surfaces) +s(Daily_Precipitacion) + s(Successional_forest)+ s(Grassland) + connected_trees + ti(X, Y) + factor(Stage) + infrutescence + DOM_score, family = binomial, data = Datasetm)
summary(GAM_wosplines)
AIC(GAM_wosplines)
plot(GAM_wosplines, pages=1, residuals=TRUE, cex = 0.7, pch = 20, shade=TRUE, scheme = 2, shade.col="khaki1")

# Compare models using Likelihood Ratio Test

#Likelihood Ratio Test (LRT)
#Higher log-likelihood values indicate a better fit to the data.
#Low p-value (< 0.05): The complex model significantly improves the fit, so it’s preferable.
#High p-value (≥ 0.05): The complex model does not significantly improve the fit, so the simpler model is adequate

lrtest(Nullmodel, GAM_additive, GAM_Step, GAM_Step2, GAM_wosplines, GAM_additive_3)

# Calculate AIC for each model
aic_GAM_Global <- AIC(GAM_Global)
aic_Nullmodel <- AIC(Nullmodel)
aic_GAM_additive<- AIC(GAM_additive)
aic_GAM_Step<- AIC(GAM_Step)
aic_GAM_additive_3<- AIC(GAM_additive_3)
aic_GAM_Step2<- AIC(GAM_Step2)

# Create a data frame to store the AIC values
aic_values <- data.frame(
  model = c("GAM_Global", "Nullmodel","GAM_additive","GAM_Step","GAM_additive_3","GAM_Step2"),
  AIC = c(aic_GAM_Global, aic_Nullmodel, aic_GAM_additive, aic_GAM_Step, aic_GAM_additive_3, aic_GAM_Step2)
)

# Calculate ΔAIC (difference from the minimum AIC)
aic_values$delta_AIC <- aic_values$AIC - min(aic_values$AIC)

# Calculate AIC weights
aic_values$AIC_weight <- exp(-0.5 * aic_values$delta_AIC) / sum(exp(-0.5 * aic_values$delta_AIC))

# Print AIC values and weights
print(aic_values)

########PLOTING RESULTS 

gam.check(GAM_Step)

# Extract fitted values (on the response scale)
fitted_values <- fitted(GAM_Step)

# Extract residuals (e.g., deviance residuals)
residuals <- residuals(GAM_Step, type = "deviance")

# Plot residuals vs. fitted values
plot(fitted_values, residuals,
     xlab = "Fitted Values",
     ylab = "Deviance Residuals",
     main = "Residuals vs. Fitted Values",
     pch = 20, col = "blue")

# Add a horizontal line at zero to highlight the reference
abline(h = 0, col = "red", lty = 2)


##QQ PLOT
deviance_residuals <- residuals(GAM_Step, type = "deviance")

# Generate a QQ-plot for the deviance residuals
qqnorm(deviance_residuals, main = "QQ-Plot of Deviance Residuals")
qqline(deviance_residuals, col = "red", lwd = 2)

qqnorm(deviance_residuals, main = "QQ-Plot of Deviance Residuals for GAM", 
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(deviance_residuals, col = "blue", lwd = 2, lty = 2)

#AUC

# Obtain predicted probabilities from the model
predicted_probs <- predict(GAM_Step, type = "response")
Tcruzi<-as.vector(Datasetm$T_cruzi)

# Calculate the AUC using pROC
roc_curve <- roc(Tcruzi, predicted_probs)
auc_value <- auc(roc_curve)

# Print the AUC
print(paste("AUC:", auc_value))

# Plot the ROC curve
plot(roc_curve, main = "ROC Curve for GAM")

# -------------------------------------------------------------------------
# Housekeeping

rm(list=ls()) # remove everything currently held in the R memory

graphics.off() # close all open graphic windows

# -------------------------------------------------------------------------

# code developed in January 2024 following Guisan et al 2017 book
# tutorial GLM for community ecology https://bios2.github.io/posts/2021-07-19-glm-community-ecology/index.html


library(ggplot2)
library(dplyr)
library(tibble)

# source files

setwd("~/Fennoscandia")

df <- read.table('Fennoscandia_1081lakes.txt',header=TRUE)
nrow(df) # 1081 lakes
head(df) 

str(df)
df$Country <- as.factor(df$Country)

#########################################
####       Models and graphs         ####
#########################################

m.full <- glm(ISOE_LAC ~ 
                poly(log1p(Alk),degree=2, raw=TRUE) + # raw=FALSE -> orthogonal polynomial
                poly(log(TP),degree=2, raw=TRUE) +
                poly(log(TN.TP),degree=2, raw=TRUE) +
                poly(log(Colour),degree=2, raw=TRUE) +
                log(Area_km2) +
                poly(log(deg.sum),degree=2, raw=TRUE) ,
                data = df, family = "binomial"
              )

m.full
summary(m.full) 

# stepwise selection based on AIC 
m.step <- step(m.full, direction="both", trace=F, k=2) # default k=2 (see below)
m.step
summary(m.step)

# select manually co-variates with negative b2 to calculate optimum and tolerance

ISOE_LAC <- as.data.frame(coefficients(m.step))
ISOE_LAC <- ISOE_LAC %>% rename_at('coefficients(m.step)', ~'ISOE_LAC.coeff')
row.names(ISOE_LAC)[row.names(ISOE_LAC) == "(Intercept)"] <- "b0"
row.names(ISOE_LAC)[row.names(ISOE_LAC) == "poly(log1p(Alk), degree = 2, raw = TRUE)1"] <- "b1.Alk"
row.names(ISOE_LAC)[row.names(ISOE_LAC) == "poly(log1p(Alk), degree = 2, raw = TRUE)2"] <- "b2.Alk"
row.names(ISOE_LAC)[row.names(ISOE_LAC) == "poly(log(TP), degree = 2, raw = TRUE)1"] <- "b1.TP"
row.names(ISOE_LAC)[row.names(ISOE_LAC) == "poly(log(TP), degree = 2, raw = TRUE)2"] <- "b2.TP"
row.names(ISOE_LAC)[row.names(ISOE_LAC) == "poly(log(Colour), degree = 2, raw = TRUE)1"] <- "b1.Colour"
row.names(ISOE_LAC)[row.names(ISOE_LAC) == "poly(log(Colour), degree = 2, raw = TRUE)2"] <- "b2.Colour"
row.names(ISOE_LAC)[row.names(ISOE_LAC) == "poly(log(deg.sum), degree = 2, raw = TRUE)1"] <- "b1.deg.sum"
row.names(ISOE_LAC)[row.names(ISOE_LAC) == "poly(log(deg.sum), degree = 2, raw = TRUE)2"] <- "b2.deg.sum"
row.names(ISOE_LAC)[row.names(ISOE_LAC) == "log(Area_km2)"] <- "b1.Area_km2"

ISOE_LAC <- as.data.frame(t(ISOE_LAC))

# using tibble
ISOE_LAC <- tibble::rownames_to_column(ISOE_LAC, "coefficients")
ISOE_LAC <- as.data.frame(ISOE_LAC)


# Select manually relevant co-variates defined above to calculate optimum and tolerance
# Alkalinity
ISOE_LAC$Alk.optimum <- exp(-(ISOE_LAC$b1.Alk)/(2*ISOE_LAC$b2.Alk))-1
ISOE_LAC$Alk.tolerance <- 1/sqrt(-2*ISOE_LAC$b2.Alk) # on log1p scale
ISOE_LAC$Alk.lwr.tol <- exp((-(ISOE_LAC$b1.Alk)/(2*ISOE_LAC$b2.Alk))-(1/sqrt(-2*ISOE_LAC$b2.Alk)))-1
ISOE_LAC$Alk.upr.tol <- exp((-(ISOE_LAC$b1.Alk)/(2*ISOE_LAC$b2.Alk))+(1/sqrt(-2*ISOE_LAC$b2.Alk)))-1
# TP
ISOE_LAC$TP.optimum <- exp(-(ISOE_LAC$b1.TP)/(2*ISOE_LAC$b2.TP))
ISOE_LAC$TP.tolerance <- 1/sqrt(-2*ISOE_LAC$b2.TP) # on log scale
ISOE_LAC$TP.lwr.tol <- exp((-(ISOE_LAC$b1.TP)/(2*ISOE_LAC$b2.TP))-(1/sqrt(-2*ISOE_LAC$b2.TP)))
ISOE_LAC$TP.upr.tol <- exp((-(ISOE_LAC$b1.TP)/(2*ISOE_LAC$b2.TP))+(1/sqrt(-2*ISOE_LAC$b2.TP)))
# TN.TP
ISOE_LAC$TN.TP.optimum <- exp(-(ISOE_LAC$b1.TN.TP)/(2*ISOE_LAC$b2.TN.TP))
ISOE_LAC$TN.TP.tolerance <- 1/sqrt(-2*ISOE_LAC$b2.TN.TP) # on log scale
ISOE_LAC$TN.TP.lwr.tol <- exp((-(ISOE_LAC$b1.TN.TP)/(2*ISOE_LAC$b2.TN.TP))-(1/sqrt(-2*ISOE_LAC$b2.TN.TP)))
ISOE_LAC$TN.TP.upr.tol <- exp((-(ISOE_LAC$b1.TN.TP)/(2*ISOE_LAC$b2.TN.TP))+(1/sqrt(-2*ISOE_LAC$b2.TN.TP)))
# Colour
ISOE_LAC$Colour.optimum <- exp(-(ISOE_LAC$b1.Colour)/(2*ISOE_LAC$b2.Colour))
ISOE_LAC$Colour.tolerance <- 1/sqrt(-2*ISOE_LAC$b2.Colour) # on log scale
ISOE_LAC$Colour.lwr.tol <- exp((-(ISOE_LAC$b1.Colour)/(2*ISOE_LAC$b2.Colour))-(1/sqrt(-2*ISOE_LAC$b2.Colour)))
ISOE_LAC$Colour.upr.tol <- exp((-(ISOE_LAC$b1.Colour)/(2*ISOE_LAC$b2.Colour))+(1/sqrt(-2*ISOE_LAC$b2.Colour)))
# deg.sum
ISOE_LAC$deg.sum.optimum <- exp(-(ISOE_LAC$b1.deg.sum)/(2*ISOE_LAC$b2.deg.sum))
ISOE_LAC$deg.sum.tolerance <- 1/sqrt(-2*ISOE_LAC$b2.deg.sum) # on log scale
ISOE_LAC$deg.sum.lwr.tol <- exp((-(ISOE_LAC$b1.deg.sum)/(2*ISOE_LAC$b2.deg.sum))-(1/sqrt(-2*ISOE_LAC$b2.deg.sum)))
ISOE_LAC$deg.sum.upr.tol <- exp((-(ISOE_LAC$b1.deg.sum)/(2*ISOE_LAC$b2.deg.sum))+(1/sqrt(-2*ISOE_LAC$b2.deg.sum)))


head(ISOE_LAC)
write.table(ISOE_LAC,"Data analyses/glm_coefficient/ISOE_LAC_glm_coeff.txt",sep="\t",row.names=FALSE)



# Adjusted D2 calibration strength for GLM based on Weisberg's (1980) formula
library(ecospat) 
ecospat.adj.D2.glm(m.full) # 0.42
ecospat.adj.D2.glm(m.step) # 0.42

######################################
#####  checking model accuracy  ######
######################################

library(PresenceAbsence)
# calibration plots 
df %>% select(EuropaBon_ID, ISOE_LAC) -> EvalData
EvalData$GLM <- predict(m.step, newdata = df, type = "response")
head(EvalData)
calibration.plot(EvalData)
calibration.plot(EvalData, which.model = 1, na.rm = FALSE, alpha = 0.05, N.bins = 10, 
                 xlab = "Probability of Occurrence", # Predicted Probability of Occurrence
                 ylab = "Observed Occurrence", # Observed Occurrence as Proportion of Sites Surveyed
                 main = NULL, color= NULL, model.names= NULL)


###########################
####  Thresholds   ########
###########################
# determination of a threshold to change continuous probability into binary data

library(ecospat)
# limited to 0.01 increment or 99 thresholds
kappa100 <- ecospat.max.kappa(EvalData$GLM, EvalData$ISOE_LAC) 
kappa100
head(kappa100)
kappa100$max.Kappa # 0.641
Kappa = kappa100$max.threshold
Kappa # 0.42

TSS100 <- ecospat.max.tss(EvalData$GLM, EvalData$ISOE_LAC)
TSS100
TSS100$max.TSS # 0.65
TSS = TSS100$max.threshold 
TSS # 0.42



###################################
#### spatial autocorrelation ######
###################################

# check for spatial autocorrelation in model residuals

# residuals
devi.res <- resid(m.step, type="deviance") # deviance residuals
boxplot(devi.res, horizontal = TRUE) # more normally distributed, good for Moran I test
resp.res.m <- resid(m.step, type="response") # response residuals (observed - fitted)
boxplot(resp.res.m, horizontal = TRUE) # same results as my calculations in dfx below

# response residuals are similarly distributed to deviance residuals 
# use either deviance or response residuals for autocorrelation tests

# save model results
df %>% select(EuropaBon_ID, X_3035, Y_3035, ISOE_LAC) -> dfx
head(dfx)
dfx$fit <- predict(m.step, newdata = df, type = "response")
pp <- predict(m.step, se.fit = TRUE) # since newdata is omitted the predictions are 
# based on the data used for the fit.
dfx$se_lwr <- with(pp, plogis(fit + qnorm(0.025)*se.fit)) # qnorm(0.025) = -1.96
dfx$se_upr <- with(pp, plogis(fit + qnorm(0.975)*se.fit)) # qnorm(0.975) = +1.96
# 95% confidence interval (range or model response uncertainty)
dfx$se_range <- dfx$se_upr - dfx$se_lwr
# glm response residuals
dfx$resp.res <- dfx$ISOE_LAC - dfx$fit

# Fitted probabilities on binary scale using TSS results
dfx$fit_TSS <-  ifelse(dfx$fit <= TSS, "0", "1")
dfx$fit_TSS <- as.numeric(dfx$fit_TSS)
# Fitted probabilities on binary scale using Kappa results 
dfx$fit_Kappa <-  ifelse(dfx$fit <= Kappa, "0", "1")
dfx$fit_Kappa <- as.numeric(dfx$fit_Kappa)

head(dfx)
str(dfx)
sum(dfx$fit_TSS) # 580
sum(dfx$fit_Kappa) # 580

# write file 
write.table(dfx,"Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_predictions.txt",sep="\t",row.names=FALSE)

 



# Transform coordinate units to km
xy <- df[,5:6]
xy$X_km <- xy$X_3035/1000 # longitude in km
xy$Y_km <- xy$Y_3035/1000 # latitude in km

# isotropic spline correlogram (time consuming so limited to 99 re-sampling)
library(ncf)

sp.corr <- spline.correlog(
  xy$X_km, # units are km
  xy$Y_km, # units are km
  dfx$resp.res,
  w = NULL,
  df = NULL,
  type = "boot",
  resamp = 99,
  npoints = 50, # every 10 km with xmax = 500 km
  save = FALSE,
  filter = FALSE,
  fw = 0,
  max.it = 25,
  xmax = 500, # units are km
  latlon = FALSE,
  na.rm = FALSE,
  quiet = FALSE
)
sp.corr
summary(sp.corr)

autocorr.plot <- plot(sp.corr) # export image .tiff 400 x 400 


#################################
######   PARTIAL PLOTS    #######
#################################
 
# Partial plots of model prediction 
library(psych) # help to calculate geometric.mean below
library(scales) # help with x axis log transformation

# Alk
sort(df$Alk)
p.Alk <- seq(0.01, max(df$Alk), by = 0.01) # minimum value set just above zero to avoid log issues

newdata = data.frame(Alk = p.Alk,
                     Area_km2 = geometric.mean(df$Area_km2),
                     Colour = geometric.mean(df$Colour),
                     TP = geometric.mean(df$TP),
                     deg.sum = geometric.mean(df$deg.sum))
head(newdata)
colSums(is.na(newdata)) # no NAs

pglm <- predict(m.step, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pglm$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- with(pglm, plogis(fit + qnorm(0.025)*se.fit)) # qnorm(0.025) = -1.96
newdata$se_upr <- with(pglm, plogis(fit + qnorm(0.975)*se.fit)) # qnorm(0.975) = +1.96

# Plot the fit and se using ggplot2
p.Alk.plot <- ggplot (newdata, aes (x = Alk, y = fit)) +
  geom_line (color = "black") +
  geom_ribbon (aes (ymin = se_lwr, ymax = se_upr), alpha = 0.2, fill = "#333333") +
  labs (x = "Alkalinity", y = "predicted response") +
  theme_bw () +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(color = "black"))+ 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "tb")  # tb top and bottom only

p.Alk.plot 

ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_p-ALK_1081lakes.png", 
       plot = p.Alk.plot ,
       width = 1539, height = 1300, 
       units = "px"
) # save image (600 x 443)


# TP
sort(df$TP)
p.TP <- seq(min(df$TP), max(df$TP), by = 0.01)

newdata = data.frame(TP = p.TP,
                     Area_km2 = geometric.mean(df$Area_km2),
                     Colour = geometric.mean(df$Colour),
                     Alk = 0.1, # use near optimum, more relevant
                     deg.sum = geometric.mean(df$deg.sum))
head(newdata)
colSums(is.na(newdata)) # no NAs

pglm <- predict(m.step, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pglm$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- with(pglm, plogis(fit + qnorm(0.025)*se.fit)) # qnorm(0.025) = -1.96
newdata$se_upr <- with(pglm, plogis(fit + qnorm(0.975)*se.fit)) # qnorm(0.975) = +1.96

# Plot the fit and se using ggplot2
p.TP.plot <- ggplot (newdata, aes (x = TP, y = fit)) +
  geom_line (color = "black") +
  geom_ribbon (aes (ymin = se_lwr, ymax = se_upr), alpha = 0.2, fill = "#333333") +
  labs (x = "Total phosphorus", y = "predicted response") +
  theme_bw () +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(color = "black")) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "tb") 

p.TP.plot   

ggsave(filename = "Data analyses/ISOE_LAC/glm//ISOE_LAC_GLM_p-TP_1081lakes.png", 
       plot = p.TP.plot ,
       width = 1539, height = 1300, 
       units = "px"
)



# Colour 
sort(df$Colour)
p.Colour <- seq(min(df$Colour), max(df$Colour), by = 0.01)

newdata = data.frame(Colour = p.Colour,
                     Area_km2 = geometric.mean(df$Area_km2),
                     TP = geometric.mean(df$TP),
                     Alk = 0.1, # use near optimum, more relevant
                     deg.sum = geometric.mean(df$deg.sum))
head(newdata)
colSums(is.na(newdata)) # no NAs

pglm <- predict(m.step, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pglm$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- with(pglm, plogis(fit + qnorm(0.025)*se.fit)) # qnorm(0.025) = -1.96
newdata$se_upr <- with(pglm, plogis(fit + qnorm(0.975)*se.fit)) # qnorm(0.975) = +1.96


# Plot the fit and se using ggplot2
p.Colour.plot <- ggplot (newdata, aes (x = Colour, y = fit)) +
  geom_line (color = "black") +
  geom_ribbon (aes (ymin = se_lwr, ymax = se_upr), alpha = 0.2, fill = "#333333") +
  labs (x = "Colour", y = "predicted response") +
  theme_bw () +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(color = "black"))+ 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "tb")  

p.Colour.plot 

ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_p-Colour_1081lakes.png", 
       plot = p.Colour.plot ,
       width = 1539, height = 1300, 
       units = "px"
)


# annual sum of degree days
sort(df$deg.sum)
p.deg.sum <- seq(min(df$deg.sum), max(df$deg.sum), by = 0.01)

newdata = data.frame(deg.sum = p.deg.sum,
                     Area_km2 = geometric.mean(df$Area_km2),
                     TP = geometric.mean(df$TP),
                     Alk = 0.1, # use near optimum, more relevant
                     Colour = geometric.mean(df$Colour))
head(newdata)
colSums(is.na(newdata)) # no NAs

pglm <- predict(m.step, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pglm$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- with(pglm, plogis(fit + qnorm(0.025)*se.fit)) # qnorm(0.025) = -1.96
newdata$se_upr <- with(pglm, plogis(fit + qnorm(0.975)*se.fit)) # qnorm(0.975) = +1.96


# Plot the fit and se using ggplot2
p.deg.sum.plot <- ggplot (newdata, aes (x = deg.sum, y = fit)) +
  geom_line (color = "black") +
  geom_ribbon (aes (ymin = se_lwr, ymax = se_upr), alpha = 0.2, fill = "#333333") +
  labs (x = "Annual sum of degree days", y = "predicted response") +
  theme_bw () +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(color = "black"))+ 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "tb")  

p.deg.sum.plot 

ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_p-deg.sum_1081lakes.png", 
       plot = p.deg.sum.plot ,
       width = 1539, height = 1300, 
       units = "px"
)


# Lake area
sort(df$Area_km2)
p.Area_km2 <- seq(min(df$Area_km2), max(df$Area_km2), by = 0.01)

newdata = data.frame(Area_km2 = p.Area_km2,
                     deg.sum = geometric.mean(df$deg.sum),
                     TP = geometric.mean(df$TP),
                     Alk = 0.1, # use near optimum, more relevant
                     Colour = geometric.mean(df$Colour))
head(newdata)
colSums(is.na(newdata)) # no NAs

pglm <- predict(m.step, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pglm$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- with(pglm, plogis(fit + qnorm(0.025)*se.fit)) # qnorm(0.025) = -1.96
newdata$se_upr <- with(pglm, plogis(fit + qnorm(0.975)*se.fit)) # qnorm(0.975) = +1.96


# Plot the fit and se using ggplot2
p.Area_km2.plot <- ggplot (newdata, aes (x = Area_km2, y = fit)) +
  geom_line (color = "black") +
  geom_ribbon (aes (ymin = se_lwr, ymax = se_upr), alpha = 0.2, fill = "#333333") +
  labs (x = "Lake area", y = "predicted response") +
  theme_bw () +
  theme(text = element_text(size = 20)) +
  theme(axis.text = element_text(color = "black"))+ 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "tb")  

p.Area_km2.plot 

ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_p-Area_km2_1081lakes.png", 
       plot = p.Area_km2.plot ,
       width = 1539, height = 1300, 
       units = "px"
) 


#################################################
### Preparation for Alk maps under varying TP ###
#################################################

# TP=5 ppb

newdata = data.frame(Alk = df$Alk,
                     Area_km2 = df$Area_km2,
                     Colour = df$Colour,
                     TP = 5,
                     deg.sum = df$deg.sum)
head(newdata)

colSums(is.na(newdata)) # no NAs

pglm <- predict(m.step, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pglm$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pglm$fit + qnorm(0.025)* pglm$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pglm$fit + qnorm(0.975)* pglm$se.fit) # qnorm(0.975) = +1.96
newdata$se_range <- newdata$se_upr - newdata$se_lwr
# Fitted probabilities on binary scale using TSS results
newdata$fit_TSS <-  ifelse(newdata$fit <= TSS, "0", "1") 
newdata$fit_TSS <- as.numeric(newdata$fit_TSS)
# Fitted probabilities on binary scale using Kappa results 
newdata$fit_Kappa <-  ifelse(newdata$fit <= Kappa, "0", "1") 
newdata$fit_Kappa <- as.numeric(newdata$fit_Kappa)

sum(newdata$fit_TSS ) # 674
sum(newdata$fit_Kappa) # 674
head(newdata)

species = which( colnames(df)=="ISOE_LAC")
df.sp.pred <- cbind(df[, c(1, 5:6, species)], newdata) # change for ISOE_LAC

head(df.sp.pred)
write.table(df.sp.pred,"Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_5.txt",sep="\t",row.names=FALSE)


# TP= 15 ppb

newdata = data.frame(Alk = df$Alk,
                     Area_km2 = df$Area_km2,
                     Colour = df$Colour,
                     TP = 15,
                     deg.sum = df$deg.sum)
head(newdata)

colSums(is.na(newdata)) # no NAs

pglm <- predict(m.step, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pglm$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pglm$fit + qnorm(0.025)* pglm$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pglm$fit + qnorm(0.975)* pglm$se.fit) # qnorm(0.975) = +1.96
newdata$se_range <- newdata$se_upr - newdata$se_lwr
# Fitted probabilities on binary scale using TSS results 
newdata$fit_TSS <-  ifelse(newdata$fit <= TSS, "0", "1") 
newdata$fit_TSS <- as.numeric(newdata$fit_TSS)
# Fitted probabilities on binary scale using Kappa results
newdata$fit_Kappa <-  ifelse(newdata$fit <= Kappa, "0", "1") 
newdata$fit_Kappa <- as.numeric(newdata$fit_Kappa)

sum(newdata$fit_TSS ) # 696
sum(newdata$fit_Kappa) # 696
head(newdata)

df.sp.pred <- cbind(df[, c(1, 5:6, species)], newdata) 

head(df.sp.pred)
write.table(df.sp.pred,"Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_15.txt",sep="\t",row.names=FALSE)

# TP= 50 ppb

newdata = data.frame(Alk = df$Alk,
                     Area_km2 = df$Area_km2,
                     Colour = df$Colour,
                     TP = 50,
                     deg.sum = df$deg.sum)
head(newdata)

colSums(is.na(newdata)) # no NAs

pglm <- predict(m.step, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pglm$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pglm$fit + qnorm(0.025)* pglm$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pglm$fit + qnorm(0.975)* pglm$se.fit) # qnorm(0.975) = +1.96
newdata$se_range <- newdata$se_upr - newdata$se_lwr
# Fitted probabilities on binary scale using TSS results 
newdata$fit_TSS <-  ifelse(newdata$fit <= TSS, "0", "1") 
newdata$fit_TSS <- as.numeric(newdata$fit_TSS)
# Fitted probabilities on binary scale using Kappa results 
newdata$fit_Kappa <-  ifelse(newdata$fit <= Kappa, "0", "1")  
newdata$fit_Kappa <- as.numeric(newdata$fit_Kappa)

sum(newdata$fit_TSS ) # 316
sum(newdata$fit_Kappa) # 316
head(newdata)

df.sp.pred <- cbind(df[, c(1, 5:6, species)], newdata) 

head(df.sp.pred)
write.table(df.sp.pred,"Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_50.txt",sep="\t",row.names=FALSE)

# TP= 100 ppb

newdata = data.frame(Alk = df$Alk,
                     Area_km2 = df$Area_km2,
                     Colour = df$Colour,
                     TP = 100,
                     deg.sum = df$deg.sum)
head(newdata)

colSums(is.na(newdata)) # no NAs

pglm <- predict(m.step, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pglm$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pglm$fit + qnorm(0.025)* pglm$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pglm$fit + qnorm(0.975)* pglm$se.fit) # qnorm(0.975) = +1.96
newdata$se_range <- newdata$se_upr - newdata$se_lwr
# Fitted probabilities on binary scale using TSS results
newdata$fit_TSS <-  ifelse(newdata$fit <= TSS, "0", "1") 
newdata$fit_TSS <- as.numeric(newdata$fit_TSS)
# Fitted probabilities on binary scale using Kappa results 
newdata$fit_Kappa <-  ifelse(newdata$fit <= Kappa, "0", "1") 
newdata$fit_Kappa <- as.numeric(newdata$fit_Kappa)

sum(newdata$fit_TSS ) # 44
sum(newdata$fit_Kappa) # 44
head(newdata)

df.sp.pred <- cbind(df[, c(1, 5:6, species)], newdata) 

head(df.sp.pred)
write.table(df.sp.pred,"Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_100.txt",sep="\t",row.names=FALSE)


# TP *0.5 (restoration scenario)

newdata = data.frame(Alk = df$Alk,
                     Area_km2 = df$Area_km2,
                     Colour = df$Colour,
                     TP = df$TP * 0.5,
                     deg.sum = df$deg.sum)
head(newdata)

colSums(is.na(newdata)) # no NAs

pglm <- predict(m.step, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pglm$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pglm$fit + qnorm(0.025)* pglm$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pglm$fit + qnorm(0.975)* pglm$se.fit) # qnorm(0.975) = +1.96
newdata$se_range <- newdata$se_upr - newdata$se_lwr
# Fitted probabilities on binary scale using TSS results 
newdata$fit_TSS <-  ifelse(newdata$fit <= TSS, "0", "1") 
newdata$fit_TSS <- as.numeric(newdata$fit_TSS)
# Fitted probabilities on binary scale using Kappa results 
newdata$fit_Kappa <-  ifelse(newdata$fit <= Kappa, "0", "1") 
newdata$fit_Kappa <- as.numeric(newdata$fit_Kappa)

sum(newdata$fit_TSS ) # 624
sum(newdata$fit_Kappa) # 624
head(newdata)

df.sp.pred <- cbind(df[, c(1, 5:6, species)], newdata) 

head(df.sp.pred)
write.table(df.sp.pred,"Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_50pc.txt",sep="\t",row.names=FALSE)



# TP *2 (pollution scenario)

newdata = data.frame(Alk = df$Alk,
                     Area_km2 = df$Area_km2,
                     Colour = df$Colour,
                     TP = df$TP * 2,
                     deg.sum = df$deg.sum)
head(newdata)

colSums(is.na(newdata)) # no NAs

pglm <- predict(m.step, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pglm$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pglm$fit + qnorm(0.025)* pglm$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pglm$fit + qnorm(0.975)* pglm$se.fit) # qnorm(0.975) = +1.96
newdata$se_range <- newdata$se_upr - newdata$se_lwr
# Fitted probabilities on binary scale using TSS results
newdata$fit_TSS <-  ifelse(newdata$fit <= TSS, "0", "1") 
newdata$fit_TSS <- as.numeric(newdata$fit_TSS)
# Fitted probabilities on binary scale using Kappa results
newdata$fit_Kappa <-  ifelse(newdata$fit <= Kappa, "0", "1")  
newdata$fit_Kappa <- as.numeric(newdata$fit_Kappa)

sum(newdata$fit_TSS ) # 474
sum(newdata$fit_Kappa) # 474
head(newdata)

df.sp.pred <- cbind(df[, c(1, 5:6, species)], newdata) 

head(df.sp.pred)
write.table(df.sp.pred,"Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_200pc.txt",sep="\t",row.names=FALSE)



# TP reference conditions
df.ref <- read.table('Data analyses/TP reference/Fennoscandia_1081lakes_ref.txt',header=TRUE)
head(df.ref)

newdata = data.frame(Alk = df.ref$Alk,
                     Area_km2 = df.ref$Area_km2,
                     Colour = df.ref$Colour,
                     TP = df.ref$TP.fit,
                     deg.sum = df.ref$deg.sum)
head(newdata)

colSums(is.na(newdata)) # no NAs

pglm <- predict(m.step, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pglm$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pglm$fit + qnorm(0.025)* pglm$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pglm$fit + qnorm(0.975)* pglm$se.fit) # qnorm(0.975) = +1.96
newdata$se_range <- newdata$se_upr - newdata$se_lwr
# Fitted probabilities on binary scale using TSS results 
newdata$fit_TSS <-  ifelse(newdata$fit <= TSS, "0", "1")
newdata$fit_TSS <- as.numeric(newdata$fit_TSS)
# Fitted probabilities on binary scale using Kappa results
newdata$fit_Kappa <-  ifelse(newdata$fit <= Kappa, "0", "1")
newdata$fit_Kappa <- as.numeric(newdata$fit_Kappa)

sum(newdata$fit_TSS ) # 641
sum(newdata$fit_Kappa) # 641
head(newdata)

species = which( colnames(df)=="ISOE_LAC")
df.sp.pred <- cbind(df.ref[, c(1, 5:6, species)], newdata)

head(df.sp.pred)
write.table(df.sp.pred,"Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_ref.txt",sep="\t",row.names=FALSE)




# TP species optimum conditions

newdata = data.frame(Alk = df$Alk, # near optimum, within range of observed values
                     Area_km2 = df$Area_km2,
                     Colour = df$Colour,
                     TP = 9.42,
                     deg.sum = df$deg.sum)
head(newdata)

colSums(is.na(newdata)) # no NAs

pglm <- predict(m.step, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pglm$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pglm$fit + qnorm(0.025)* pglm$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pglm$fit + qnorm(0.975)* pglm$se.fit) # qnorm(0.975) = +1.96
newdata$se_range <- newdata$se_upr - newdata$se_lwr
# Fitted probabilities on binary scale using TSS results 
newdata$fit_TSS <-  ifelse(newdata$fit <= TSS, "0", "1")
newdata$fit_TSS <- as.numeric(newdata$fit_TSS)
# Fitted probabilities on binary scale using Kappa results
newdata$fit_Kappa <-  ifelse(newdata$fit <= Kappa, "0", "1")
newdata$fit_Kappa <- as.numeric(newdata$fit_Kappa)

sum(newdata$fit_TSS ) # 711
sum(newdata$fit_Kappa) # 711
head(newdata)

species = which( colnames(df)=="ISOE_LAC")
df.sp.pred <- cbind(df.ref[, c(1, 5:6, species)], newdata)

head(df.sp.pred)
write.table(df.sp.pred,"Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_optimum.txt",sep="\t",row.names=FALSE)





# Reference TP conditions - current projections: common, gain and loss
df.TP.ref <- read.table('Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_ref.txt',header=TRUE)
df.TP.ref %>% select(fit_TSS) -> TSS.ref
head(TSS.ref)
TSS.ref <- TSS.ref %>% rename_at('fit_TSS', ~'fit_TSS_ref')

df.TP.current <- read.table('Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_predictions.txt',header=TRUE)
df.TP.current <- cbind(df.TP.current, TSS.ref)
df.TP.current$ref_current <- (df.TP.current$fit_TSS_ref) - (df.TP.current$fit_TSS)
df.TP.current$common <- ifelse(df.TP.current$fit_TSS_ref == 1 
                                     & df.TP.current$fit_TSS == 1, "1", "0")
df.TP.current$common <- as.numeric(df.TP.current$common)
sum(df.TP.current$common) # 574
df.TP.current$gain <- ifelse(df.TP.current$ref_current < 0, "1", "0")
df.TP.current$gain <- as.numeric(df.TP.current$gain)
sum(df.TP.current$gain) # 6
df.TP.current$loss <- ifelse(df.TP.current$ref_current > 0, "1", "0")
df.TP.current$loss <- as.numeric(df.TP.current$loss)
sum(df.TP.current$loss) # 67
df.TP.current <- df.TP.current %>% select(-c(fit_TSS_ref, ref_current))
head(df.TP.current)
# write file 
write.table(df.TP.current,"Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_predictions.txt",sep="\t",row.names=FALSE)

# turnover
Sorensen <-  (sum(df.TP.current$gain) + sum(df.TP.current$loss)) / (
  2* sum(df.TP.current$common) + sum(df.TP.current$gain) + sum(df.TP.current$loss) )
Simpson <- min(sum(df.TP.current$gain),sum(df.TP.current$loss)) / (
  sum(df.TP.current$common) + min(sum(df.TP.current$gain),sum(df.TP.current$loss)) )
nestedness <- Sorensen - Simpson
Sorensen # 0.06
Simpson # 0.01
nestedness # 0.05




# Optimum TP conditions - current projections: common, gain and loss
df.TP.opt <- read.table('Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_optimum.txt',header=TRUE)
df.TP.opt %>% select(fit_TSS) -> TSS.opt
head(TSS.opt)
TSS.opt <- TSS.opt %>% rename_at('fit_TSS', ~'fit_TSS_opt')

df.TP.current <- read.table('Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_predictions.txt',header=TRUE)
df.TP.current <- cbind(df.TP.current, TSS.opt)
df.TP.current$opt_current <- (df.TP.current$fit_TSS_opt) - (df.TP.current$fit_TSS)
df.TP.current$common <- ifelse(df.TP.current$fit_TSS_opt == 1 
                               & df.TP.current$fit_TSS == 1, "1", "0")
df.TP.current$common <- as.numeric(df.TP.current$common)
sum(df.TP.current$common) # 580
df.TP.current$gain <- ifelse(df.TP.current$opt_current < 0, "1", "0")
df.TP.current$gain <- as.numeric(df.TP.current$gain)
sum(df.TP.current$gain) # 0
df.TP.current$loss <- ifelse(df.TP.current$opt_current > 0, "1", "0")
df.TP.current$loss <- as.numeric(df.TP.current$loss)
sum(df.TP.current$loss) # 131
df.TP.current <- df.TP.current %>% select(-c(fit_TSS_opt, opt_current))
head(df.TP.current)
# write file 
write.table(df.TP.current,"Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_optimum_current.txt",sep="\t",row.names=FALSE)

# turnover
Sorensen <-  (sum(df.TP.current$gain) + sum(df.TP.current$loss)) / (
  2* sum(df.TP.current$common) + sum(df.TP.current$gain) + sum(df.TP.current$loss) )
Simpson <- min(sum(df.TP.current$gain),sum(df.TP.current$loss)) / (
  sum(df.TP.current$common) + min(sum(df.TP.current$gain),sum(df.TP.current$loss)) )
nestedness <- Sorensen - Simpson
Sorensen # 0.10
Simpson # 0.00
nestedness # 0.10

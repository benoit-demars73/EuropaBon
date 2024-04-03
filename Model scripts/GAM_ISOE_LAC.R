# GAM_Fennoscandia

# -------------------------------------------------------------------------
# Housekeeping

rm(list=ls()) # remove everything currently held in the R memory

graphics.off() # close all open graphic windows

# -------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(mgcv) 

setwd("~/Fennoscandia")

df <- read.table('Fennoscandia_1081lakes.txt',header=TRUE)
nrow(df) # 1081 lakes
head(df) 

str(df)
df$Country <- as.factor(df$Country)


#########################################
####       Models and graphs         ####
#########################################

g.selec <- gam(ISOE_LAC ~  log(Area_km2) +
          s(log1p(Alk)) +
          s(log(Colour)) +
          s(log(TP)) +
      #    s(log(TN.TP)) +
          s(log(deg.sum)) ,
     #   s(X_3035,Y_3035) ,
          data = df, 
          family = binomial(link=logit), 
          method = "REML",
          trace = TRUE)
g.selec
summary(g.selec) # expl. deviance = 45.2%, full convergence after 10 iterations
AIC(g.selec)
# Choosing best model according to AIC and hand selection 
AIC (g0,g1,g2,g3)

# model checking
coef(g.selec) # provide coefficients 
gam.check(g.selec) 
concurvity(g.selec, full = TRUE) # overall concurvity 0.45-0.72
concurvity(g.selec, full = FALSE) # pair-wise generally fine, highest Colour v. TP: 0.46


######################################
#####  checking model accuracy  ######
######################################
library(PresenceAbsence)
# calibration plots 
df %>% select(EuropaBon_ID, ISOE_LAC) -> EvalData
EvalData$GAM <- predict(g.selec, newdata = df, type = "response")
head(EvalData) # only reports model fit here
calibration.plot(EvalData)
calibration.plot(EvalData, which.model = 1, na.rm = FALSE, alpha = 0.05, N.bins = 10, 
                 xlab = "Probability of Occurrence", # Predicted Probability of Occurrence
                 ylab = "Observed Occurrence", # Observed Occurrence as Proportion of Sites Surveyed
                 main = NULL, color= NULL, model.names= NULL)

# determination of a threshold to change continuous probability into binary data
# Guisan et al 2017, p. 256
accu <- presence.absence.accuracy(EvalData, which.model=1, threshold = 11)
accu [, -c(1,2)] <- signif(accu [, -c(1,2)], digits=2)
accu [c("threshold", "PCC", "sensitivity", "specificity", "Kappa", "AUC")] 

# AUC=0.95, scale of judgement Guisan et al 2017 p. 263
library(ecospat)
kappa100 <- ecospat.max.kappa(EvalData$GAM, EvalData$ISOE_LAC)
kappa100
head(kappa100)
kappa100$max.Kappa # 0.652
Kappa = kappa100$max.threshold # 0.42
Kappa

TSS100 <- ecospat.max.tss(EvalData$GAM, EvalData$ISOE_LAC)
TSS100
TSS100$max.TSS # 0.652
TSS = TSS100$max.threshold # 0.42
TSS

###################################
#### spatial autocorrelation ######
###################################

# check for spatial autocorrelation in model residuals

# Extract the fit (response), lwr and upr CI and range CI and append to dataframe
species = which( colnames(df)=="ISOE_LAC")
df.x <- cbind(df[, c(1, 5:6, species)]) # ID, X, Y, species

pp <- predict(g.selec, type = "link", se.fit = TRUE)
df.x$fit <- plogis(pp$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
df.x$se_lwr <- plogis(pp$fit + qnorm(0.025)* pp$se.fit) # qnorm(0.025) = -1.96
df.x$se_upr <- plogis(pp$fit + qnorm(0.975)* pp$se.fit) # qnorm(0.975) = +1.96
df.x$se_range <- df.x$se_upr - df.x$se_lwr

# Fitted probabilities on binary scale using TSS results 
df.x$fit_TSS <-  ifelse(df.x$fit <= TSS, "0", "1") 
df.x$fit_TSS <- as.numeric(df.x$fit_TSS)
# Fitted probabilities on binary scale using Kappa results
df.x$fit_Kappa <-  ifelse(df.x$fit <= Kappa, "0", "1") 
df.x$fit_Kappa <- as.numeric(df.x$fit_Kappa)

head(df.x)
str(df.x)
sum(df.x$fit_TSS) # 572
sum(df.x$fit_Kappa) # 572

# write file 
write.table(df.x,"Data analyses/ISOE_LAC/mgcv_gam/ISOE_LAC_GAM_predictions.txt",
            sep="\t",row.names=FALSE)


# gam residuals 
# also available: "pearson","scaled.pearson", "working" residuals
devi.res <- residuals(g.selec, type = "deviance")
boxplot(devi.res, horizontal = TRUE)
resp.res <- residuals(g.selec, type = "response")
boxplot(resp.res, horizontal = TRUE)


# Transform coordinate units to km
xy <- df[,5:6]
xy$X_km <- xy$X_3035/1000 # longitude in km
xy$Y_km <- xy$Y_3035/1000 # latitude in km


library(ncf)
sp.corr <- spline.correlog(
  xy$X_km, # units are km
  xy$Y_km, # units are km
  devi.res,
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

autocorr.plot <- plot(sp.corr) # save as .tiff


#############################
#####    PARTIAL PLOT   #####
#############################

# Partial plots of model prediction 
library(psych) # help to calculate geometric.mean below
library(scales) # help with x axis log transformation
# custom margins
# par(mar=c(4,4,1,1))

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

pgam <- predict(g.selec, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pgam$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pgam$fit + qnorm(0.025)* pgam$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pgam$fit + qnorm(0.975)* pgam$se.fit) # qnorm(0.975) = +1.96


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
ggsave(filename = "Data analyses/ISOE_LAC/mgcv_gam/ISOE_LAC_GAM_p-ALK_1081lakes.png", 
       plot = p.Alk.plot ,
       width = 1539, height = 1300, 
       units = "px"
 )


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

pgam <- predict(g.selec, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pgam$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pgam$fit + qnorm(0.025)* pgam$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pgam$fit + qnorm(0.975)* pgam$se.fit) # qnorm(0.975) = +1.96


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

ggsave(filename = "Data analyses/ISOE_LAC/mgcv_gam/ISOE_LAC_GAM_p-TP_1081lakes.png", 
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

pgam <- predict(g.selec, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pgam$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pgam$fit + qnorm(0.025)* pgam$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pgam$fit + qnorm(0.975)* pgam$se.fit) # qnorm(0.975) = +1.96


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

ggsave(filename = "Data analyses/ISOE_LAC/mgcv_gam/ISOE_LAC_GAM_p-Colour_1081lakes.png", 
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

pgam <- predict(g.selec, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pgam$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pgam$fit + qnorm(0.025)* pgam$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pgam$fit + qnorm(0.975)* pgam$se.fit) # qnorm(0.975) = +1.96


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

ggsave(filename = "Data analyses/ISOE_LAC/mgcv_gam/ISOE_LAC_GAM_p-deg.sum_1081lakes.png", 
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

pgam <- predict(g.selec, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pgam$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pgam$fit + qnorm(0.025)* pgam$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pgam$fit + qnorm(0.975)* pgam$se.fit) # qnorm(0.975) = +1.96


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

ggsave(filename = "Data analyses/ISOE_LAC/mgcv_gam/ISOE_LAC_GAM_p-Area_km2_1081lakes.png", 
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

pgam <- predict(g.selec, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pgam$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pgam$fit + qnorm(0.025)* pgam$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pgam$fit + qnorm(0.975)* pgam$se.fit) # qnorm(0.975) = +1.96
newdata$se_range <- newdata$se_upr - newdata$se_lwr
# Fitted probabilities on binary scale using TSS results 
newdata$fit_TSS <-  ifelse(newdata$fit <= TSS, "0", "1") 
newdata$fit_TSS <- as.numeric(newdata$fit_TSS)
# Fitted probabilities on binary scale using Kappa results 
newdata$fit_Kappa <-  ifelse(newdata$fit <= Kappa, "0", "1") 
newdata$fit_Kappa <- as.numeric(newdata$fit_Kappa)

sum(newdata$fit_TSS ) # 663
sum(newdata$fit_Kappa) # 663
head(newdata)

df.sp.pred <- cbind(df[, c(1, 5:6, species)], newdata)

head(df.sp.pred)
write.table(df.sp.pred,"Data analyses/ISOE_LAC/mgcv_gam/ISOE_LAC_GAM_TP_5.txt",sep="\t",row.names=FALSE)


# TP= 15 ppb

newdata = data.frame(Alk = df$Alk,
                     Area_km2 = df$Area_km2,
                     Colour = df$Colour,
                     TP = 15,
                     deg.sum = df$deg.sum)
head(newdata)

colSums(is.na(newdata)) # no NAs

pgam <- predict(g.selec, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pgam$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pgam$fit + qnorm(0.025)* pgam$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pgam$fit + qnorm(0.975)* pgam$se.fit) # qnorm(0.975) = +1.96
newdata$se_range <- newdata$se_upr - newdata$se_lwr
# Fitted probabilities on binary scale using TSS results
newdata$fit_TSS <-  ifelse(newdata$fit <= TSS, "0", "1") 
newdata$fit_TSS <- as.numeric(newdata$fit_TSS)
# Fitted probabilities on binary scale using Kappa results
newdata$fit_Kappa <-  ifelse(newdata$fit <= Kappa, "0", "1") 
newdata$fit_Kappa <- as.numeric(newdata$fit_Kappa)

sum(newdata$fit_TSS ) # 693
sum(newdata$fit_Kappa) # 693
head(newdata)

df.sp.pred <- cbind(df[, c(1, 5:6, species)], newdata)

head(df.sp.pred)
write.table(df.sp.pred,"Data analyses/ISOE_LAC/mgcv_gam/ISOE_LAC_GAM_TP_15.txt",sep="\t",row.names=FALSE)

# TP= 50 ppb

newdata = data.frame(Alk = df$Alk,
                     Area_km2 = df$Area_km2,
                     Colour = df$Colour,
                     TP = 50,
                     deg.sum = df$deg.sum)
head(newdata)

colSums(is.na(newdata)) # no NAs

pgam <- predict(g.selec, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pgam$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pgam$fit + qnorm(0.025)* pgam$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pgam$fit + qnorm(0.975)* pgam$se.fit) # qnorm(0.975) = +1.96
newdata$se_range <- newdata$se_upr - newdata$se_lwr
# Fitted probabilities on binary scale using TSS results 
newdata$fit_TSS <-  ifelse(newdata$fit <= TSS, "0", "1") 
newdata$fit_TSS <- as.numeric(newdata$fit_TSS)
# Fitted probabilities on binary scale using Kappa results 
newdata$fit_Kappa <-  ifelse(newdata$fit <= Kappa, "0", "1") 
newdata$fit_Kappa <- as.numeric(newdata$fit_Kappa)

sum(newdata$fit_TSS ) # 325
sum(newdata$fit_Kappa) # 325
head(newdata)

df.sp.pred <- cbind(df[, c(1, 5:6, species)], newdata)

head(df.sp.pred)
write.table(df.sp.pred,"Data analyses/ISOE_LAC/mgcv_gam/ISOE_LAC_GAM_TP_50.txt",sep="\t",row.names=FALSE)

# TP= 100 ppb

newdata = data.frame(Alk = df$Alk,
                     Area_km2 = df$Area_km2,
                     Colour = df$Colour,
                     TP = 100,
                     deg.sum = df$deg.sum)
head(newdata)

colSums(is.na(newdata)) # no NAs

pgam <- predict(g.selec, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pgam$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pgam$fit + qnorm(0.025)* pgam$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pgam$fit + qnorm(0.975)* pgam$se.fit) # qnorm(0.975) = +1.96
newdata$se_range <- newdata$se_upr - newdata$se_lwr
# Fitted probabilities on binary scale using TSS results 
newdata$fit_TSS <-  ifelse(newdata$fit <= TSS, "0", "1") 
newdata$fit_TSS <- as.numeric(newdata$fit_TSS)
# Fitted probabilities on binary scale using Kappa results 
newdata$fit_Kappa <-  ifelse(newdata$fit <= Kappa, "0", "1") 
newdata$fit_Kappa <- as.numeric(newdata$fit_Kappa)

sum(newdata$fit_TSS ) # 133
sum(newdata$fit_Kappa) # 133
head(newdata)

df.sp.pred <- cbind(df[, c(1, 5:6, species)], newdata)

head(df.sp.pred)
write.table(df.sp.pred,"Data analyses/ISOE_LAC/mgcv_gam/ISOE_LAC_GAM_TP_100.txt",sep="\t",row.names=FALSE)


# TP *0.5 (restoration scenario)

newdata = data.frame(Alk = df$Alk,
                     Area_km2 = df$Area_km2,
                     Colour = df$Colour,
                     TP = df$TP * 0.5,
                     deg.sum = df$deg.sum)
head(newdata)

colSums(is.na(newdata)) # no NAs

pgam <- predict(g.selec, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pgam$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pgam$fit + qnorm(0.025)* pgam$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pgam$fit + qnorm(0.975)* pgam$se.fit) # qnorm(0.975) = +1.96
newdata$se_range <- newdata$se_upr - newdata$se_lwr
# Fitted probabilities on binary scale using TSS results
newdata$fit_TSS <-  ifelse(newdata$fit <= TSS, "0", "1")
newdata$fit_TSS <- as.numeric(newdata$fit_TSS)
# Fitted probabilities on binary scale using Kappa results
newdata$fit_Kappa <-  ifelse(newdata$fit <= Kappa, "0", "1")
newdata$fit_Kappa <- as.numeric(newdata$fit_Kappa)

sum(newdata$fit_TSS ) # 629
sum(newdata$fit_Kappa) # 629
head(newdata)

df.sp.pred <- cbind(df[, c(1, 5:6, species)], newdata)

head(df.sp.pred)
write.table(df.sp.pred,"Data analyses/ISOE_LAC/mgcv_gam/ISOE_LAC_GAM_TP_50pc.txt",sep="\t",row.names=FALSE)



# TP *2 (pollution scenario)

newdata = data.frame(Alk = df$Alk,
                     Area_km2 = df$Area_km2,
                     Colour = df$Colour,
                     TP = df$TP * 2,
                     deg.sum = df$deg.sum)
head(newdata)

colSums(is.na(newdata)) # no NAs

pgam <- predict(g.selec, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pgam$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pgam$fit + qnorm(0.025)* pgam$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pgam$fit + qnorm(0.975)* pgam$se.fit) # qnorm(0.975) = +1.96
newdata$se_range <- newdata$se_upr - newdata$se_lwr
# Fitted probabilities on binary scale using TSS results 
newdata$fit_TSS <-  ifelse(newdata$fit <= TSS, "0", "1")
newdata$fit_TSS <- as.numeric(newdata$fit_TSS)
# Fitted probabilities on binary scale using Kappa results
newdata$fit_Kappa <-  ifelse(newdata$fit <= Kappa, "0", "1")
newdata$fit_Kappa <- as.numeric(newdata$fit_Kappa)

sum(newdata$fit_TSS ) # 470
sum(newdata$fit_Kappa) # 470
head(newdata)

df.sp.pred <- cbind(df[, c(1, 5:6, species)], newdata)

head(df.sp.pred)
write.table(df.sp.pred,"Data analyses/ISOE_LAC/mgcv_gam/ISOE_LAC_GAM_TP_200pc.txt",sep="\t",row.names=FALSE)





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

pgam <- predict(g.selec, newdata , se.fit = T, type = "link")

# Add the fit and se to the new data frame
newdata$fit <- plogis(pgam$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$se_lwr <- plogis(pgam$fit + qnorm(0.025)* pgam$se.fit) # qnorm(0.025) = -1.96
newdata$se_upr <- plogis(pgam$fit + qnorm(0.975)* pgam$se.fit) # qnorm(0.975) = +1.96
newdata$se_range <- newdata$se_upr - newdata$se_lwr
# Fitted probabilities on binary scale using TSS results 
newdata$fit_TSS <-  ifelse(newdata$fit <= TSS, "0", "1")
newdata$fit_TSS <- as.numeric(newdata$fit_TSS)
# Fitted probabilities on binary scale using Kappa results
newdata$fit_Kappa <-  ifelse(newdata$fit <= Kappa, "0", "1")
newdata$fit_Kappa <- as.numeric(newdata$fit_Kappa)

sum(newdata$fit_TSS ) # 638
sum(newdata$fit_Kappa) # 638
head(newdata)

species = which( colnames(df)=="ISOE_LAC")
df.sp.pred <- cbind(df.ref[, c(1, 5:6, species)], newdata)

head(df.sp.pred)
write.table(df.sp.pred,"Data analyses/ISOE_LAC/mgcv_gam/ISOE_LAC_GAM_TP_ref.txt",sep="\t",row.names=FALSE)


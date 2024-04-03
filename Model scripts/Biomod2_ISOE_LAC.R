# BIOMOD_Fennoscandia

# -------------------------------------------------------------------------
# Housekeeping

rm(list=ls()) # remove everything currently held in the R memory

graphics.off() # close all open graphic windows

# -------------------------------------------------------------------------

# code developed in January 2024 following Guisan et al 2017, p.367
# with updated names of functions for biomod2 4.2-4

library(pacman)
library(ggplot2)
library(dplyr)
library(tidyr)
library(biomod2) 
library(ggtext)


# ********************************************************************* #
# if working in Windows, set a short worKing directory path for biomod2
# otherwise it won't work
# ********************************************************************* #


setwd("C:/BIOMOD2")
df.biomod <- read.table('Fennoscandia_1081lakes_ref.txt',header=TRUE)
head(df.biomod)

# TP scenarios
df.biomod$log.TP.ref <- log(df.biomod$TP.fit) # reference conditions
df.biomod$log.TP.half <- log(df.biomod$TP*0.5)
df.biomod$log.TP.double <- log(df.biomod$TP*2)
df.biomod$log.TP.5 =  log(5)
df.biomod$log.TP.15 = log(15)
df.biomod$log.TP.50 = log(50)
df.biomod$log.TP.100 = log(100)
head(df.biomod)
# Biomod formating
ISOE_LAC_data <- BIOMOD_FormatingData (
        resp.name = "Isoetes.lacustris",
        resp.var = df.biomod$ISOE_LAC,
        expl.var = cbind(df.biomod[, c(122, 124:125, 127, 130)]),
        resp.xy = cbind(df.biomod[, c(5:6)]),
        dir.name = "C:/BIOMOD2/ISOE_LAC"
        ) # no data has been set aside for model evaluation
ISOE_LAC_data
plot(ISOE_LAC_data)

# Biomod model options
# GLM formulae not using orthogonal polynomial to recover max., optimum & tolerance
glmForm <- 'Isoetes.lacustris ~ poly(log.Alk, degree=2, raw=TRUE) +
                poly(log.TP, degree=2, raw=TRUE) +
                poly(log.Colour, degree=2, raw=TRUE) +
                log.Area_km2 +
                poly(log.deg.sum,degree=2, raw=TRUE)'
glmForm

ISOE_LAC_opt <- BIOMOD_ModelingOptions(
        GLM = list (type = NULL,
                    interaction.level = NULL,
                    myFormula = formula(glmForm), 
                    step = 'AIC',
                    test=NULL, # default for stepwise selection 'AIC'
                    family = 'binomial'(link = 'logit'),
                    mustart = 0.5, # default
                    control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE), # default
        GAM = list (algo = "GAM_mgcv")
                    )
) 

## check changes was done
ISOE_LAC_opt

# Individual models
ISOE_LAC_models <- BIOMOD_Modeling(
        bm.format = ISOE_LAC_data,
        modeling.id = 'AllModels',
        models = c('GLM', 'GAM'), # 'GBM', 'RF', ...
        bm.options = ISOE_LAC_opt,
        CV.strategy = 'random', # 'strat', 'env', ...
    #    CV.strat = 'both', # stratified sampling along x , y geographical coordinates
# species prevalence can be kept constant (p. 287), 
# argument balance set to default "prevalence", in 'strat' or 'env' for k-fold cross validation 
        CV.nb.rep = 10,
        CV.perc = 0.7,  # proportion for calibration 0-1
        CV.do.full.models = FALSE,
        metric.eval = c('TSS','ROC', 'KAPPA'), 
        var.import = 5,
        weights = NULL,
        prevalence = NULL,
        seed.val = 42
)

ISOE_LAC_models
summary(ISOE_LAC_models)

# get list of models
get_built_models(ISOE_LAC_models, full.name = NULL, PA = NULL, run = NULL, algo = NULL)

# MODEL PERFORMANCE
# get model evaluation scores (ROC, TSS, KAPPA, sensitivity for each run)
ISOE_LAC_models_scores <- get_evaluations(ISOE_LAC_models)
ISOE_LAC_models_scores
head(ISOE_LAC_models_scores)

write.table(ISOE_LAC_models_scores,
            "C:/BIOMOD2/ISOE_LAC/Evalx10_glm-gam_ISOE_LAC_Fennoscandia_1081lakes.txt",
            sep="\t",row.names=FALSE)

# Represent evaluation scores 
bm_PlotEvalMean(bm.out = ISOE_LAC_models) # axes titles correct !

# graphics of evaluation scores
bm_PlotEvalMean(bm.out = ISOE_LAC_models,
                dataset = 'calibration',
                group.by = "algo", # group by 'algo', 'run' or 'expl.var*
                metric.eval = c("TSS","ROC"),
                do.plot = TRUE
                ) # axes titles inversed
bm_PlotEvalMean(bm.out = ISOE_LAC_models,
                dataset = 'validation',
                group.by = "algo",
                metric.eval = c("TSS","ROC"),
                ) 
bm_PlotEvalBoxplot(bm.out = ISOE_LAC_models,
                   group.by = c('algo', 'run'),
                   dataset = 'calibration'
                   ) # all datapoints rather than boxplots
bm_PlotEvalBoxplot(bm.out = ISOE_LAC_models,
                   group.by = c('algo','run'),
                   dataset = 'validation'
                   )


# get variable importance
# need to state number of permutations with var.import in model setting above
ISOE_LAC_models_var.import <- get_variables_importance (ISOE_LAC_models)
ISOE_LAC_models_var.import

write.table(ISOE_LAC_models_var.import,
            "C:/BIOMOD2/ISOE_LAC/Var_importance_glm_ISOE_LAC_Fennoscandia_1081lakes.txt",
            sep="\t",row.names=FALSE)

# calculate the mean of variable importance by algorithm
apply(ISOE_LAC_models_var.import, c(1,2), mean) # create an array with NA values

var.import.tibble <- ISOE_LAC_models_var.import %>% 
  group_by(algo, expl.var) %>% 
  summarise(mean_value = mean(var.imp)) # calculate averages
var.import.table <- var.import.tibble %>% 
  pivot_wider(names_from = algo, values_from = mean_value)
var.import.table # produce table (tibble)



# load the produced models

ISOE_LAC_GLM <- get_built_models(ISOE_LAC_models, 
                 full.name = NULL, 
                 PA = NULL, 
                 run = NULL, 
                 algo = "GLM")

partial.plots.mean <- bm_PlotResponseCurves(
  bm.out = ISOE_LAC_models,
  models.chosen = ISOE_LAC_GLM,
  new.env = get_formal_data(ISOE_LAC_models, "expl.var"),
  show.variables = get_formal_data(ISOE_LAC_models, "expl.var.names"),
  fixed.var = "mean", # equivalent to geometric mean since co-variables on log scale
  # alkalinity on geometric mean 0.39+1 which is not good for siliceous species 
  do.bivariate = FALSE,
  do.plot = TRUE,
  do.progress = TRUE
)
partial.plots.mean

# save partial plot lines
df.partial.plots.mean <- data.frame(partial.plots.mean$tab)
head(df.partial.plots.mean)
write.table(df.partial.plots.mean,"C:/BIOMOD2/ISOE_LAC/ISOE_LAC_GLM_partial_plots_mean.txt",
            sep="\t",row.names=FALSE)


# ---------------------------------------------------------------------- #

ISOE_LAC_GBM <- get_built_models(ISOE_LAC_models, 
                                 full.name = NULL, 
                                 PA = NULL, 
                                 run = NULL, 
                                 algo = "GBM")

bm_PlotResponseCurves(
  bm.out = ISOE_LAC_models,
  models.chosen = ISOE_LAC_GBM,
  new.env = get_formal_data(ISOE_LAC_models, "expl.var"),
  show.variables = get_formal_data(ISOE_LAC_models, "expl.var.names"),
  fixed.var = "mean",
  do.bivariate = FALSE,
  do.plot = TRUE,
  do.progress = TRUE
)

# ----------------------------------------------------------------- #

ISOE_LAC_RF <- get_built_models(ISOE_LAC_models, 
                                 full.name = NULL, 
                                 PA = NULL, 
                                 run = NULL, 
                                 algo = "RF")

bm_PlotResponseCurves(
  bm.out = ISOE_LAC_models,
  models.chosen = ISOE_LAC_RF,
  new.env = get_formal_data(ISOE_LAC_models, "expl.var"),
  show.variables = get_formal_data(ISOE_LAC_models, "expl.var.names"),
  fixed.var = "mean",
  do.bivariate = FALSE,
  do.plot = TRUE,
  do.progress = TRUE
)


# ------------------------------------------------------------- #

ISOE_LAC_GAM <- get_built_models(ISOE_LAC_models, 
                                 full.name = NULL, 
                                 PA = NULL, 
                                 run = NULL, 
                                 algo = "GAM")

bm_PlotResponseCurves(
  bm.out = ISOE_LAC_models,
  models.chosen = ISOE_LAC_GAM,
  new.env = get_formal_data(ISOE_LAC_models, "expl.var"),
  show.variables = get_formal_data(ISOE_LAC_models, "expl.var.names"),
  fixed.var = "mean",
  do.bivariate = FALSE,
  do.plot = TRUE,
  do.progress = TRUE
)

# save partial plot lines
df.partial.plots.mean <- data.frame(partial.plots.mean$tab)
head(df.partial.plots.mean)
write.table(df.partial.plots.mean,"C:/BIOMOD2/ISOE_LAC/ISOE_LAC_GAM_partial_plots_mean.txt",
            sep="\t",row.names=FALSE)


#####################################################
################ ENSEMBLE MODEL #####################
#####################################################

# NOTE: the directory pathway must be short otherwise it does not work (BIOMOD2 4.2-4)
# see reported bug: https://github.com/biomodhub/biomod2/issues/170
ISOE_LAC_EM <- BIOMOD_EnsembleModeling(bm.mod = ISOE_LAC_models,
                          models.chosen = 'all',
                          em.by = 'all',
                          em.algo = c('EMwmean', 'EMca', 'EMcv'), # 'EMcv', 'EMci', 'EMmedian', 'EMca', 'EMmean'
                          metric.select = c('TSS'),
                          metric.select.thresh = c(0.6),
                          metric.eval = c('TSS', 'ROC', 'KAPPA'),
                          var.import = 5,
                          EMci.alpha = 0.05,
                       #   EMwmean.decay = 'proportional'
                          )

ISOE_LAC_EM
# provide a list of models computed: mean, cv, ci_inf, ci_sup
# the full model names can be used to plot and save individual model results

# get uncertainties
# help with ci_cv https://github.com/biomodhub/biomod2/issues/345

# Returns list of models long names for BIOMOD.ensemble.models.out
get_built_models( 
  ISOE_LAC_EM,
  full.name = NULL,
  merged.by.algo = NULL,
  merged.by.run = NULL,
  merged.by.PA = NULL,
  filtered.by = NULL,
  algo = NULL
)

# EMwmean (weighted mean by the TSS of the evaluation runs)

get_predictions(
  ISOE_LAC_EM,
  evaluation = FALSE,
  full.name = 'Isoetes.lacustris_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo',
  merged.by.algo = NULL,
  merged.by.run = NULL,
  merged.by.PA = NULL,
  filtered.by = NULL,
  algo = NULL,
  model.as.col = FALSE
)

# EM.wmean with shorter writing and df format
EM.wmean.df <- as.data.frame(get_predictions(ISOE_LAC_EM,
  full.name = 'Isoetes.lacustris_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo'))
EM.wmean.df <- EM.wmean.df %>% rename_at('pred', ~'EM.wmean')
EM.wmean.df <- subset(EM.wmean.df, select = -c(1:7)) 
EM.wmean.df <- EM.wmean.df/1000
head(EM.wmean.df)

# EMca
EM.ca.df <- as.data.frame(get_predictions(ISOE_LAC_EM,
  full.name = 'Isoetes.lacustris_EMcaByTSS_mergedData_mergedRun_mergedAlgo'))
EM.ca.df <- EM.ca.df %>% rename_at('pred', ~'EM.ca')
EM.ca.df  <- subset(EM.ca.df , select = -c(1:7)) 
head(EM.ca.df)

# EMcv
EM.cv.df <- as.data.frame(get_predictions(ISOE_LAC_EM,
  full.name = 'Isoetes.lacustris_EMcvByTSS_mergedData_mergedRun_mergedAlgo'))
EM.cv.df <- EM.cv.df %>% rename_at('pred', ~'EM.cv')
EM.cv.df  <- subset(EM.cv.df , select = -c(1:7)) 
head(EM.cv.df)

# EM.summary
species = which( colnames(df)=="ISOE_LAC") # not working ???
EM.df <- cbind(df.biomod[, c(1, 5:6, 22)]) # ID, X, Y, species
head(EM.df)
EM.proj <- bind_cols(EM.wmean.df, EM.ca.df, EM.cv.df)

EM.proj <- cbind(EM.df, EM.proj)
EM.proj$EM.sd <- EM.proj$EM.cv / 100 * EM.proj$EM.wmean
head(EM.proj)
nrow(EM.proj)

# remove NA values (cv when wmean=0 or tending to zero)
EM.proj <- EM.proj %>% mutate(
                              EM.cv = ifelse(is.na(EM.cv), 0, EM.cv),
                              EM.sd = ifelse(is.na(EM.sd), 0, EM.sd)
                              )
library(ecospat)
kappa100 <- ecospat.max.kappa(EM.proj$EM.wmean, EM.proj$ISOE_LAC)
kappa100
kappa100$max.Kappa # 0.64
Kappa = kappa100$max.threshold # 0.45

TSS100 <- ecospat.max.tss(EM.proj$EM.wmean, EM.proj$ISOE_LAC)
TSS100
TSS100$max.TSS # 0.65
TSS = TSS100$max.threshold # 0.45

# ***************************************************************** #
#   why are the thresholds much higher than for the GLM model ????
# ***************************************************************** #

# this is becasuse the y scale is different, 0-1 but rescaled to plot different algo
# adapted from Elith et al (2005)

# Fitted probabilities on binary scale using TSS results (response threshold)
EM.proj$fit_TSS <-  ifelse(EM.proj$EM.wmean <= TSS, "0", "1") 
EM.proj$fit_TSS <- as.numeric(EM.proj$fit_TSS)
# Fitted probabilities on binary scale using Kappa results (response threshold)
EM.proj$fit_Kappa <-  ifelse(EM.proj$EM.wmean <= Kappa, "0", "1")
EM.proj$fit_Kappa <- as.numeric(EM.proj$fit_Kappa)

head(EM.proj)
write.table(EM.proj,
            "C:/BIOMOD2/ISOE_LAC/EM_ISOE_LAC_Fennoscandia_1081lakes.txt",
            sep="\t",row.names=FALSE)


# just a note for all possible set up in BIOMOD.ensemble.models above
em.algo.long <- c('EMmean' = 'Mean of probabilities', 
                  'EMcv' = 'Coef of variation of probabilities', 
                  'EMciInf' = 'Confidence Interval (Inf)',
                  'EMciSup' = 'Confidence Interval (Sup)', 
                  'EMmedian' = 'Median of probabilities',
                  'EMca' = 'Committee averaging', 
                  'EMwmean' = 'Probabilities weighting mean')


# Get evaluation scores 
ISOE_LAC_EM_scores <- get_evaluations(ISOE_LAC_EM)

write.table(ISOE_LAC_EM_scores,
            "C:/BIOMOD2/ISOE_LAC/EM_eval_glm_ISOE_LAC_Fennoscandia_1081lakes.txt",
            sep="\t",row.names=FALSE)

# Get prediction curve 
get_predictions(ISOE_LAC_EM) # print the data
EM.df <- as.data.frame(get_predictions(ISOE_LAC_EM)) # save the data as df
head(EM.df)

write.table(EM.df,
            "C:/BIOMOD2/ISOE_LAC/EM_fit_glm_gam_ISOE_LAC_Fennoscandia_1081lakes.txt",
            sep="\t",row.names=FALSE)



# Get variables importance
ISOE_LAC_EM_var.import <- get_variables_importance(ISOE_LAC_EM)

write.table(ISOE_LAC_EM_var.import,
      "C:/BIOMOD2/ISOE_LAC/Var_importance_EM_glm_gam_ISOE_LAC_Fennoscandia_1081lakes.txt",
      sep="\t",row.names=FALSE)


# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = ISOE_LAC_EM, group.by = 'full.name')
bm_PlotEvalBoxplot(bm.out = ISOE_LAC_EM, group.by = c('full.name', 'full.name'))
bm_PlotVarImpBoxplot(bm.out = ISOE_LAC_EM, group.by = c('expl.var', 'full.name', 'full.name'))
bm_PlotVarImpBoxplot(bm.out = ISOE_LAC_EM, group.by = c('expl.var', 'algo', 'merged.by.run'))
bm_PlotVarImpBoxplot(bm.out = ISOE_LAC_EM, group.by = c('algo', 'expl.var', 'merged.by.run'))

# Represent response curves and save partial plot lines
partial.plots.mean.EM <- bm_PlotResponseCurves(bm.out = ISOE_LAC_EM, 
                      models.chosen = get_built_models(ISOE_LAC_EM)[c(1, 6, 7)],
                      fixed.var = 'mean') # geometric mean of co-variates since all log transformed
partial.plots.mean.EM

df.partial.plots.mean.EM <- data.frame(partial.plots.mean.EM$tab)
head(df.partial.plots.mean.EM)

write.table(df.partial.plots.mean.EM,
            "C:/BIOMOD2/ISOE_LAC/ISOE_LAC_partial_plots_mean_EM.txt",
            sep="\t",
            row.names=FALSE)

# Predicted model runs with uncertainties for current conditions
ISOE_LAC_runs <- BIOMOD_Projection(
                  bm.mod = ISOE_LAC_models,
                  proj.name = "Potamogeton.lucens.EM",
                  new.env = cbind(df.biomod[, c(122, 124:125, 127, 130)]),
                  new.env.xy = cbind(df.biomod[, c(5:6)]),
                  models.chosen = "all",
                  metric.binary = "TSS",
                  metric.filter = "TSS",
                  compress = TRUE,
                  build.clamping.mask = FALSE, # not an issue here
                  nb.cpu = 1,
                  seed.val = NULL,
                  output.format = ".RData"
                  
)

ISOE_LAC_runs

plot(ISOE_LAC_runs) # plot nice maps of the 10 model runs with continuous probabilities (0-1000)
get_projected_models(ISOE_LAC_runs) # just list the 10 model runs

get_predictions(ISOE_LAC_runs) # print the data
projections.df <- as.data.frame(get_predictions(ISOE_LAC_runs)) # save the data as df
head(projections.df)

write.table(projections.df,
            "C:/BIOMOD2/ISOE_LAC/ISOE_LAC_projections_10runs_glm_gam.txt",
            sep="\t",
            row.names=FALSE)




# Predicted ensemble model with uncertainties for current conditions

# reported bug: https://github.com/biomodhub/biomod2/issues/371 
# had the same issue. Cause had nothing to do with em.by setting above
# either set bm.proj with 'BIOMOD_Projection.out' object, or
# new.env and new.env.xy, but not both
# In any case turn off metric.binary and metric.filter
ISOE_LAC_pEM <- BIOMOD_EnsembleForecasting(
  bm.em = ISOE_LAC_EM, # must be a 'BIOMOD.ensemble.models.out' object
  bm.proj = ISOE_LAC_runs, # must be a 'BIOMOD_Projection.out' object
  proj.name = "Potamogeton.lucens.EM",
#  new.env = cbind(df.biomod[, c(122, 124:125, 127, 130)]),
#  new.env.xy = cbind(df.biomod[, c(5:6)]),
  models.chosen = "all",
  metric.binary = NULL, # suppressing metric.binary and metric.filter allows the model to run
  metric.filter = NULL,
  compress = TRUE,
  nb.cpu = 1,
  na.rm = TRUE,
  output.format = ".RData",
  do.stack = FALSE
)


get_predictions(ISOE_LAC_pEM) # same as ISOE_LAC_EM, all good !



#####################################################
################    RUN SCENARIOS    ################
#####################################################

# Add columns to dataframe to run scenarios (e.g. TP=5, 15, 50, 100 or *0.5, *2)

head(df.biomod[, c(122, 124:125, 127, 130, 136)]) # check co-variates and TP scenario

df.biomod$log.TP <- df.biomod$log.TP.ref # reset log.TP values to scenario log.TP

ISOE_LAC_EM.scenario <- BIOMOD_EnsembleForecasting(
  bm.em = ISOE_LAC_EM, # must be a 'BIOMOD.ensemble.models.out' object
  bm.proj = NULL,
  
  proj.name = "Isoetes.lacustris.EM",
  new.env = cbind(df.biomod[, c(122, 124:125, 127, 130)]), #  TP scenario column removed
  new.env.xy = cbind(df.biomod[, c(5:6)]),
  models.chosen = "all",
  metric.binary = NULL, # suppressing metric.binary and metric.filter allows the model to run
  metric.filter = NULL,
  compress = TRUE,
  nb.cpu = 1,
  na.rm = TRUE,
  output.format = ".RData",
  do.stack = FALSE
)

ISOE_LAC_EM.scenario # provide list of model runs with full names
plot(ISOE_LAC_EM.scenario) # plots three maps_ EM.wmean, EM.cv, EM.ca

get_predictions(ISOE_LAC_EM.scenario ) # print the data

# EM.wmean with shorter writing and df format
EM.wmean.df <- as.data.frame(get_predictions(ISOE_LAC_EM.scenario,
                                             full.name = 'Isoetes.lacustris_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo'))
EM.wmean.df <- EM.wmean.df %>% rename_at('pred', ~'EM.wmean')
EM.wmean.df <- subset(EM.wmean.df, select = -c(1:7)) 
EM.wmean.df <- EM.wmean.df/1000
head(EM.wmean.df)

# EMca
EM.ca.df <- as.data.frame(get_predictions(ISOE_LAC_EM.scenario,
                                          full.name = 'Isoetes.lacustris_EMcaByTSS_mergedData_mergedRun_mergedAlgo'))
EM.ca.df <- EM.ca.df %>% rename_at('pred', ~'EM.ca')
EM.ca.df  <- subset(EM.ca.df , select = -c(1:7)) 
head(EM.ca.df)

# EMcv
EM.cv.df <- as.data.frame(get_predictions(ISOE_LAC_EM.scenario,
                                          full.name = 'Isoetes.lacustris_EMcvByTSS_mergedData_mergedRun_mergedAlgo'))
EM.cv.df <- EM.cv.df %>% rename_at('pred', ~'EM.cv')
EM.cv.df  <- subset(EM.cv.df , select = -c(1:7)) 
head(EM.cv.df)

# EM.summary
EM.df <- cbind(df.biomod[, c(1, 5:6, 22)]) # ID, X, Y, ISOE_LAC
head(EM.df)
EM.scenario <- bind_cols(EM.wmean.df, EM.ca.df, EM.cv.df)

EM.scenario <- cbind(EM.df, EM.scenario)
EM.scenario$EM.sd <- EM.scenario$EM.cv / 100 * EM.scenario$EM.wmean

# remove NA values (cv when wmean=0 or tending to zero)
EM.scenario <- EM.scenario %>% mutate(
                    EM.cv = ifelse(is.na(EM.cv), 0, EM.cv),
                    EM.sd = ifelse(is.na(EM.sd), 0, EM.sd)
                    )
# use same TSS and Kappa threshold as for the fitted EM 
# Fitted probabilities on binary scale using TSS results (response threshold)
EM.scenario$fit_TSS <-  ifelse(EM.scenario$EM.wmean <= TSS, "0", "1") 
EM.scenario$fit_TSS <- as.numeric(EM.scenario$fit_TSS)
# Fitted probabilities on binary scale using Kappa results (response threshold)
EM.scenario$fit_Kappa <-  ifelse(EM.scenario$EM.wmean <= Kappa, "0", "1") 
EM.scenario$fit_Kappa <- as.numeric(EM.scenario$fit_Kappa)
head(EM.scenario)

write.table(EM.scenario,
            "C:/BIOMOD2/ISOE_LAC/ISOE_LAC_EM.scenario_TP_ref.txt",
            sep="\t",
            row.names=FALSE)


#####################################################
################ RANGE ENVELOP MODEL ################
#####################################################


# Compute SRE for several quantile values
sre.100 <- bm_SRE(resp.var = df.biomod$ISOE_LAC,
                  expl.var = cbind(df.biomod[, c(122, 124:125, 127, 130)]),
                  new.env = cbind(df.biomod[, c(122, 124:125, 127, 130)]),
                  quant = 0)
sre.095 <- bm_SRE(resp.var = df.biomod$ISOE_LAC,
                  expl.var = cbind(df.biomod[, c(122, 124:125, 127, 130)]),
                  new.env = cbind(df.biomod[, c(122, 124:125, 127, 130)]),
                  quant = 0.025)
sre.090 <- bm_SRE(resp.var = df.biomod$ISOE_LAC,
                  expl.var = cbind(df.biomod[, c(122, 124:125, 127, 130)]),
                  new.env = cbind(df.biomod[, c(122, 124:125, 127, 130)]),
                  quant = 0.05)

# Visualize results
res <- c(df.biomod$ISOE_LAC, sre.100, sre.095, sre.090)
names(res) <- c("Original distribution", "Full data calibration", "Over 95 percent", "Over 90 percent")
plot(res, zlim = c(0, 1)) # warning: "zlim" is not a graphical parameter


# SDM maps gradient of pollution pressure

##################################################


# ------------------------------------------------------------------------ #
# Housekeeping

rm(list=ls()) # remove everything currently held in the R memory

graphics.off() # close all open graphic windows
# ------------------------------------------------------------------------ #



library(ggplot2)
library(dplyr)
library(sf)
library(raster)

setwd("~/Fennoscandia")

#####################
##### ISOE_LAC  #####
#####################


x.width = 2400 * 1.5

x.height = 2400 * 1.5

x.data = read.delim("Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_predictions.txt")
head(x.data)
x.data.sf = x.data %>% st_as_sf(coords = c("X_3035", "Y_3035")) %>% st_set_crs(3035)


### Base maps


st_layers("Data/gadm41_NOR.gpkg")

x.gadm.no = read_sf("Data/gadm41_NOR.gpkg", layer = "ADM_ADM_0") 
x.gadm.se = read_sf("Data/gadm41_SWE.gpkg", layer = "ADM_ADM_0")
x.gadm.fi = read_sf("Data/gadm41_FIN.gpkg", layer = "ADM_ADM_0")

x.gadm = rbind(x.gadm.no, x.gadm.se, x.gadm.fi)

x.gadm.1 = x.gadm %>%
  st_transform(3035)


### ISOE_LAC observation

x.data.sf.obs  = x.data.sf %>%
  arrange(ISOE_LAC) %>% # to bring points with highest values to the front
  mutate(ISOE_LAC_fac = as.factor(ISOE_LAC))


x.plot = ggplot() +
  
  geom_sf(data = x.gadm.1 %>% st_as_sfc(), fill = "white") +
  
  geom_sf(data = x.data.sf.obs,
          
          aes( colour = ISOE_LAC_fac , size = 1)) + # observed, control size of point
  
  guides(size = "none") + # size and colour arguments could not just be separated by a comma
  
  guides(colour = guide_legend(override.aes = list(size=5))) +
  
  scale_colour_viridis_d(direction=-1)   + # with a d for discrete
  
  labs(colour = "occurrence") + # edit legend title
  
  theme(axis.text=element_text(size=40), # edit size of text and elements
        legend.text=element_text(size=40),
        legend.title=element_text(size=40),
        legend.key.size = unit(1, "cm")) 


x.plot




ggsave(filename = "Data analyses/ISOE_LAC/ISOE_LAC_observation.png", 
       
       plot = x.plot,
       
       units = "px",
       
       width = x.width,
       
       height = x.height)

####################################################
### ISOE_LAC fitted 

x.data.sf.fit  = x.data.sf %>% # fitted
  arrange(fit) # fitted

x.plot = ggplot() +
  
  geom_sf(data = x.gadm.1 %>% st_as_sfc(), fill = "white") +
  
  geom_sf(data = x.data.sf.fit, # fitted
          
          aes( colour = fit  , size = 1)) + # fitted
  
  guides(size = "none") + 
  
  scale_colour_viridis_c(direction=-1, n.breaks = 3) + # with a c for continuous
  
  labs(colour = "probability") + # fitted
  
  theme(axis.text=element_text(size=40), 
        legend.text=element_text(size=40),
        legend.title=element_text(size=40)
  ) 

x.plot



ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_fit.png", 
       
       plot = x.plot,
       
       units = "px",
       
       width = x.width,
       
       height = x.height)

#################################################
### ISOE_LAC se_range (95% confidence interval)

x.data.sf.unc  = x.data.sf %>%  # uncertainty
  arrange(se_range) # uncertainty

x.plot = ggplot() +
  
  geom_sf(data = x.gadm.1 %>% st_as_sfc(), fill = "white") +
  
  geom_sf(data = x.data.sf.unc, # uncertainty
          
          aes( colour = se_range  , size = 1)) + # uncertainty
  
  guides(size = "none") + 
  
  scale_colour_viridis_c(direction=-1, n.breaks = 3) + # with a c for continuous
  
  labs(colour = "uncertainty") + # uncertainty
  
  theme(axis.text=element_text(size=40), 
        legend.text=element_text(size=40),
        legend.title=element_text(size=40)
  ) 



x.plot



ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_uncertainty.png", 
       
       plot = x.plot,
       
       units = "px",
       
       width = x.width,
       
       height = x.height)


#################################
##### ISOE_LAC  TSS         #####
#################################


x.width = 2400 * 1.5

x.height = 2400 * 1.5

x.data = read.delim("Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_predictions.txt")
head(x.data)
x.data.sf = x.data %>% st_as_sf(coords = c("X_3035", "Y_3035")) %>% st_set_crs(3035)


### Base maps


st_layers("Data/gadm41_NOR.gpkg")

x.gadm.no = read_sf("Data/gadm41_NOR.gpkg", layer = "ADM_ADM_0") 
x.gadm.se = read_sf("Data/gadm41_SWE.gpkg", layer = "ADM_ADM_0")
x.gadm.fi = read_sf("Data/gadm41_FIN.gpkg", layer = "ADM_ADM_0")

x.gadm = rbind(x.gadm.no, x.gadm.se, x.gadm.fi)

x.gadm.1 = x.gadm %>%
  st_transform(3035)

### ISOE_LAC fitted probability of occurrence filtered by TSS (binary)

x.data.sf.TSS  = x.data.sf %>%
  arrange(fit_TSS) %>% 
  mutate(fit_TSS_fac = as.factor(fit_TSS))

x.plot = ggplot() +
  
  geom_sf(data = x.gadm.1 %>% st_as_sfc(), fill = "white") +
  
  geom_sf(data = x.data.sf.TSS,
          
          aes( colour = fit_TSS_fac  , size = 1)) + 
  
  guides(size = "none") + 
  
  guides(colour = guide_legend(override.aes = list(size=5))) +
  
  scale_colour_viridis_d(direction=-1)   + 
  
  labs(colour = "fit TSS") + 
  
  theme(axis.text=element_text(size=40), 
        legend.text=element_text(size=40),
        legend.title=element_text(size=40),
        legend.key.size = unit(1, "cm")) 

x.plot



ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_fit_TSS.png", 
       
       plot = x.plot,
       
       units = "px",
       
       width = x.width,
       
       height = x.height)


#################################
##### ISOE_LAC  KAPPA       #####
#################################


x.width = 2400 * 1.5

x.height = 2400 * 1.5

x.data = read.delim("Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_predictions.txt")
head(x.data)
x.data.sf = x.data %>% st_as_sf(coords = c("X_3035", "Y_3035")) %>% st_set_crs(3035)


### Base maps


st_layers("Data/gadm41_NOR.gpkg")

x.gadm.no = read_sf("Data/gadm41_NOR.gpkg", layer = "ADM_ADM_0") 
x.gadm.se = read_sf("Data/gadm41_SWE.gpkg", layer = "ADM_ADM_0")
x.gadm.fi = read_sf("Data/gadm41_FIN.gpkg", layer = "ADM_ADM_0")

x.gadm = rbind(x.gadm.no, x.gadm.se, x.gadm.fi)

x.gadm.1 = x.gadm %>%
  st_transform(3035)

### ISOE_LAC fitted probability of occurrence filtered by Kappa (binary)

x.data.sf.Kappa  = x.data.sf %>%
  arrange(fit_Kappa) %>% 
  mutate(fit_Kappa_fac = as.factor(fit_Kappa))

x.plot = ggplot() +
  
  geom_sf(data = x.gadm.1 %>% st_as_sfc(), fill = "white") +
  
  geom_sf(data = x.data.sf.Kappa,
          
          aes( colour = fit_Kappa_fac  , size = 1)) + 
  
  guides(size = "none") + 
  
  guides(colour = guide_legend(override.aes = list(size=5))) +
  
  scale_colour_viridis_d(direction=-1)   + 
  
  labs(colour = "fit Kappa") + 
  
  theme(axis.text=element_text(size=40), 
        legend.text=element_text(size=40),
        legend.title=element_text(size=40),
        legend.key.size = unit(1, "cm")) 

x.plot



ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_fit_Kappa.png", 
       
       plot = x.plot,
       
       units = "px",
       
       width = x.width,
       
       height = x.height)


###############################
##### ISOE_LAC  TP= 5 ppb #####
###############################


x.width = 2400 * 1.5

x.height = 2400 * 1.5

x.data = read.delim("Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_5.txt")
head(x.data)
x.data.sf = x.data %>% st_as_sf(coords = c("X_3035", "Y_3035")) %>% st_set_crs(3035)


### Base maps


st_layers("Data/gadm41_NOR.gpkg")

x.gadm.no = read_sf("Data/gadm41_NOR.gpkg", layer = "ADM_ADM_0") 
x.gadm.se = read_sf("Data/gadm41_SWE.gpkg", layer = "ADM_ADM_0")
x.gadm.fi = read_sf("Data/gadm41_FIN.gpkg", layer = "ADM_ADM_0")

x.gadm = rbind(x.gadm.no, x.gadm.se, x.gadm.fi)

x.gadm.1 = x.gadm %>%
  st_transform(3035)

### ISOE_LAC fitted probability of occurrence filtered by TSS (binary)

x.data.sf.TSS  = x.data.sf %>%
  arrange(fit_TSS) %>% 
  mutate(fit_TSS_fac = as.factor(fit_TSS))

x.plot = ggplot() +
  
  geom_sf(data = x.gadm.1 %>% st_as_sfc(), fill = "white") +
  
  geom_sf(data = x.data.sf.TSS,
          
          aes( colour = fit_TSS_fac  , size = 1)) + 
  
  guides(size = "none") + 
  
  guides(colour = guide_legend(override.aes = list(size=5))) +
  
  scale_colour_viridis_d(direction=-1)   + 
  
  labs(colour = "fit TSS") + 
  
  theme(axis.text=element_text(size=40), 
        legend.text=element_text(size=40),
        legend.title=element_text(size=40),
        legend.key.size = unit(1, "cm")) 

x.plot



ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_fit_TSS_TP_5ppb.png", 
       
       plot = x.plot,
       
       units = "px",
       
       width = x.width,
       
       height = x.height)


###############################
##### ISOE_LAC  TP= 15 ppb #####
###############################


x.width = 2400 * 1.5

x.height = 2400 * 1.5

x.data = read.delim("Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_15.txt")
head(x.data)
x.data.sf = x.data %>% st_as_sf(coords = c("X_3035", "Y_3035")) %>% st_set_crs(3035)


### Base maps


st_layers("Data/gadm41_NOR.gpkg")

x.gadm.no = read_sf("Data/gadm41_NOR.gpkg", layer = "ADM_ADM_0") 
x.gadm.se = read_sf("Data/gadm41_SWE.gpkg", layer = "ADM_ADM_0")
x.gadm.fi = read_sf("Data/gadm41_FIN.gpkg", layer = "ADM_ADM_0")

x.gadm = rbind(x.gadm.no, x.gadm.se, x.gadm.fi)

x.gadm.1 = x.gadm %>%
  st_transform(3035)

### ISOE_LAC fitted probability of occurrence filtered by TSS (binary)

x.data.sf.TSS  = x.data.sf %>%
  arrange(fit_TSS) %>% 
  mutate(fit_TSS_fac = as.factor(fit_TSS))

x.plot = ggplot() +
  
  geom_sf(data = x.gadm.1 %>% st_as_sfc(), fill = "white") +
  
  geom_sf(data = x.data.sf.TSS,
          
          aes( colour = fit_TSS_fac  , size = 1)) + 
  
  guides(size = "none") + 
  
  guides(colour = guide_legend(override.aes = list(size=5))) +
  
  scale_colour_viridis_d(direction=-1)   + 
  
  labs(colour = "fit TSS") + 
  
  theme(axis.text=element_text(size=40), 
        legend.text=element_text(size=40),
        legend.title=element_text(size=40),
        legend.key.size = unit(1, "cm")) 

x.plot



ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_fit_TSS_TP_15ppb.png", 
       
       plot = x.plot,
       
       units = "px",
       
       width = x.width,
       
       height = x.height)

################################
##### ISOE_LAC  TP= 50 ppb #####
################################


x.width = 2400 * 1.5

x.height = 2400 * 1.5

x.data = read.delim("Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_50.txt")
head(x.data)
x.data.sf = x.data %>% st_as_sf(coords = c("X_3035", "Y_3035")) %>% st_set_crs(3035)


### Base maps


st_layers("Data/gadm41_NOR.gpkg")

x.gadm.no = read_sf("Data/gadm41_NOR.gpkg", layer = "ADM_ADM_0") 
x.gadm.se = read_sf("Data/gadm41_SWE.gpkg", layer = "ADM_ADM_0")
x.gadm.fi = read_sf("Data/gadm41_FIN.gpkg", layer = "ADM_ADM_0")

x.gadm = rbind(x.gadm.no, x.gadm.se, x.gadm.fi)

x.gadm.1 = x.gadm %>%
  st_transform(3035)

### ISOE_LAC fitted probability of occurrence filtered by TSS (binary)

x.data.sf.TSS  = x.data.sf %>%
  arrange(fit_TSS) %>% 
  mutate(fit_TSS_fac = as.factor(fit_TSS))

x.plot = ggplot() +
  
  geom_sf(data = x.gadm.1 %>% st_as_sfc(), fill = "white") +
  
  geom_sf(data = x.data.sf.TSS,
          
          aes( colour = fit_TSS_fac  , size = 1)) + 
  
  guides(size = "none") + 
  
  guides(colour = guide_legend(override.aes = list(size=5))) +
  
  scale_colour_viridis_d(direction=-1)   + 
  
  labs(colour = "fit TSS") + 
  
  theme(axis.text=element_text(size=40), 
        legend.text=element_text(size=40),
        legend.title=element_text(size=40),
        legend.key.size = unit(1, "cm")) 

x.plot



ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_fit_TSS_TP_50ppb.png", 
       
       plot = x.plot,
       
       units = "px",
       
       width = x.width,
       
       height = x.height)


#################################
##### ISOE_LAC  TP= 100 ppb #####
#################################


x.width = 2400 * 1.5

x.height = 2400 * 1.5

x.data = read.delim("Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_100.txt")
head(x.data)
x.data.sf = x.data %>% st_as_sf(coords = c("X_3035", "Y_3035")) %>% st_set_crs(3035)


### Base maps


st_layers("Data/gadm41_NOR.gpkg")

x.gadm.no = read_sf("Data/gadm41_NOR.gpkg", layer = "ADM_ADM_0") 
x.gadm.se = read_sf("Data/gadm41_SWE.gpkg", layer = "ADM_ADM_0")
x.gadm.fi = read_sf("Data/gadm41_FIN.gpkg", layer = "ADM_ADM_0")

x.gadm = rbind(x.gadm.no, x.gadm.se, x.gadm.fi)

x.gadm.1 = x.gadm %>%
  st_transform(3035)

### ISOE_LAC fitted probability of occurrence filtered by TSS (binary)

x.data.sf.TSS  = x.data.sf %>%
  arrange(fit_TSS) %>% 
  mutate(fit_TSS_fac = as.factor(fit_TSS))

x.plot = ggplot() +
  
  geom_sf(data = x.gadm.1 %>% st_as_sfc(), fill = "white") +
  
  geom_sf(data = x.data.sf.TSS,
          
          aes( colour = fit_TSS_fac  , size = 1)) + 
  
  guides(size = "none") + 
  
  guides(colour = guide_legend(override.aes = list(size=5))) +
  
  scale_colour_viridis_d(direction=-1)   + 
  
  labs(colour = "fit TSS") + 
  
  theme(axis.text=element_text(size=40), 
        legend.text=element_text(size=40),
        legend.title=element_text(size=40),
        legend.key.size = unit(1, "cm")) 

x.plot



ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_fit_TSS_TP_100ppb.png", 
       
       plot = x.plot,
       
       units = "px",
       
       width = x.width,
       
       height = x.height)



###############################
##### ISOE_LAC  TP 50 %   #####
###############################


x.width = 2400 * 1.5

x.height = 2400 * 1.5

x.data = read.delim("Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_50pc.txt")
head(x.data)
x.data.sf = x.data %>% st_as_sf(coords = c("X_3035", "Y_3035")) %>% st_set_crs(3035)


### Base maps


st_layers("Data/gadm41_NOR.gpkg")

x.gadm.no = read_sf("Data/gadm41_NOR.gpkg", layer = "ADM_ADM_0") 
x.gadm.se = read_sf("Data/gadm41_SWE.gpkg", layer = "ADM_ADM_0")
x.gadm.fi = read_sf("Data/gadm41_FIN.gpkg", layer = "ADM_ADM_0")

x.gadm = rbind(x.gadm.no, x.gadm.se, x.gadm.fi)

x.gadm.1 = x.gadm %>%
  st_transform(3035)

### ISOE_LAC fitted probability of occurrence filtered by TSS (binary)

x.data.sf.TSS  = x.data.sf %>%
  arrange(fit_TSS) %>% 
  mutate(fit_TSS_fac = as.factor(fit_TSS))

x.plot = ggplot() +
  
  geom_sf(data = x.gadm.1 %>% st_as_sfc(), fill = "white") +
  
  geom_sf(data = x.data.sf.TSS,
          
          aes( colour = fit_TSS_fac  , size = 1)) + 
  
  guides(size = "none") + 
  
  guides(colour = guide_legend(override.aes = list(size=5))) +
  
  scale_colour_viridis_d(direction=-1)   + 
  
  labs(colour = "fit TSS") + 
  
  theme(axis.text=element_text(size=40), 
        legend.text=element_text(size=40),
        legend.title=element_text(size=40),
        legend.key.size = unit(1, "cm")) 

x.plot



ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_fit_TSS_TP_50pc.png", 
       
       plot = x.plot,
       
       units = "px",
       
       width = x.width,
       
       height = x.height)



#################################
##### ISOE_LAC  TP= 200 %   #####
#################################


x.width = 2400 * 1.5

x.height = 2400 * 1.5

x.data = read.delim("Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_200pc.txt")
head(x.data)
x.data.sf = x.data %>% st_as_sf(coords = c("X_3035", "Y_3035")) %>% st_set_crs(3035)


### Base maps


st_layers("Data/gadm41_NOR.gpkg")

x.gadm.no = read_sf("Data/gadm41_NOR.gpkg", layer = "ADM_ADM_0") 
x.gadm.se = read_sf("Data/gadm41_SWE.gpkg", layer = "ADM_ADM_0")
x.gadm.fi = read_sf("Data/gadm41_FIN.gpkg", layer = "ADM_ADM_0")

x.gadm = rbind(x.gadm.no, x.gadm.se, x.gadm.fi)

x.gadm.1 = x.gadm %>%
  st_transform(3035)

### ISOE_LAC fitted probability of occurrence filtered by TSS (binary)

x.data.sf.TSS  = x.data.sf %>%
  arrange(fit_TSS) %>% 
  mutate(fit_TSS_fac = as.factor(fit_TSS))

x.plot = ggplot() +
  
  geom_sf(data = x.gadm.1 %>% st_as_sfc(), fill = "white") +
  
  geom_sf(data = x.data.sf.TSS,
          
          aes( colour = fit_TSS_fac  , size = 1)) + 
  
  guides(size = "none") + 
  
  guides(colour = guide_legend(override.aes = list(size=5))) +
  
  scale_colour_viridis_d(direction=-1)   + 
  
  labs(colour = "fit TSS") + 
  
  theme(axis.text=element_text(size=40), 
        legend.text=element_text(size=40),
        legend.title=element_text(size=40),
        legend.key.size = unit(1, "cm")) 

x.plot



ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_fit_TSS_TP_200pc.png", 
       
       plot = x.plot,
       
       units = "px",
       
       width = x.width,
       
       height = x.height)



##############################################
##### ISOE_LAC  TP reference conditions  #####
##############################################


x.width = 2400 * 1.5

x.height = 2400 * 1.5

x.data = read.delim("Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_ref.txt")
head(x.data)
x.data.sf = x.data %>% st_as_sf(coords = c("X_3035", "Y_3035")) %>% st_set_crs(3035)


### Base maps


st_layers("Data/gadm41_NOR.gpkg")

x.gadm.no = read_sf("Data/gadm41_NOR.gpkg", layer = "ADM_ADM_0") 
x.gadm.se = read_sf("Data/gadm41_SWE.gpkg", layer = "ADM_ADM_0")
x.gadm.fi = read_sf("Data/gadm41_FIN.gpkg", layer = "ADM_ADM_0")

x.gadm = rbind(x.gadm.no, x.gadm.se, x.gadm.fi)

x.gadm.1 = x.gadm %>%
  st_transform(3035)

### ISOE_LAC fitted probability of occurrence filtered by TSS (binary)

x.data.sf.TSS  = x.data.sf %>%
  arrange(fit_TSS) %>% 
  mutate(fit_TSS_fac = as.factor(fit_TSS))

x.plot = ggplot() +
  
  geom_sf(data = x.gadm.1 %>% st_as_sfc(), fill = "white") +
  
  geom_sf(data = x.data.sf.TSS,
          
          aes( colour = fit_TSS_fac  , size = 1)) + 
  
  guides(size = "none") + 
  
  guides(colour = guide_legend(override.aes = list(size=5))) +
  
  scale_colour_viridis_d(direction=-1)   + 
  
  labs(colour = "fit TSS") + 
  
  theme(axis.text=element_text(size=40), 
        legend.text=element_text(size=40),
        legend.title=element_text(size=40),
        legend.key.size = unit(1, "cm")) 

x.plot



ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_ref.png", 
       
       plot = x.plot,
       
       units = "px",
       
       width = x.width,
       
       height = x.height)


##################################
##### ISOE_LAC  status gain  #####
##################################


x.width = 2400 * 1.5

x.height = 2400 * 1.5

x.data = read.delim("Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_predictions.txt")
head(x.data)
x.data.sf = x.data %>% st_as_sf(coords = c("X_3035", "Y_3035")) %>% st_set_crs(3035)


### Base maps


st_layers("Data/gadm41_NOR.gpkg")

x.gadm.no = read_sf("Data/gadm41_NOR.gpkg", layer = "ADM_ADM_0") 
x.gadm.se = read_sf("Data/gadm41_SWE.gpkg", layer = "ADM_ADM_0")
x.gadm.fi = read_sf("Data/gadm41_FIN.gpkg", layer = "ADM_ADM_0")

x.gadm = rbind(x.gadm.no, x.gadm.se, x.gadm.fi)

x.gadm.1 = x.gadm %>%
  st_transform(3035)

### ISOE_LAC fitted probability of occurrence filtered by TSS (binary)

x.data.sf.TSS  = x.data.sf %>%
  arrange(gain) %>% 
  mutate(gain_fac = as.factor(gain))

x.plot = ggplot() +
  
  geom_sf(data = x.gadm.1 %>% st_as_sfc(), fill = "white") +
  
  geom_sf(data = x.data.sf.TSS,
          
          aes( colour = gain_fac  , size = 1)) + 
  
  guides(size = "none") + 
  
  guides(colour = guide_legend(override.aes = list(size=5))) +
  
  scale_colour_viridis_d(direction=-1)   + 
  
  labs(colour = "gain") + 
  
  theme(axis.text=element_text(size=40), 
        legend.text=element_text(size=40),
        legend.title=element_text(size=40),
        legend.key.size = unit(1, "cm")) 

x.plot



ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_gain_ref.png", 
       
       plot = x.plot,
       
       units = "px",
       
       width = x.width,
       
       height = x.height)



##################################
##### ISOE_LAC  status loss  #####
##################################


x.width = 2400 * 1.5

x.height = 2400 * 1.5

x.data = read.delim("Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_predictions.txt")
head(x.data)
x.data.sf = x.data %>% st_as_sf(coords = c("X_3035", "Y_3035")) %>% st_set_crs(3035)


### Base maps


st_layers("Data/gadm41_NOR.gpkg")

x.gadm.no = read_sf("Data/gadm41_NOR.gpkg", layer = "ADM_ADM_0") 
x.gadm.se = read_sf("Data/gadm41_SWE.gpkg", layer = "ADM_ADM_0")
x.gadm.fi = read_sf("Data/gadm41_FIN.gpkg", layer = "ADM_ADM_0")

x.gadm = rbind(x.gadm.no, x.gadm.se, x.gadm.fi)

x.gadm.1 = x.gadm %>%
  st_transform(3035)

### ISOE_LAC fitted probability of occurrence filtered by TSS (binary)

x.data.sf.TSS  = x.data.sf %>%
  arrange(loss) %>% 
  mutate(loss_fac = as.factor(loss))

x.plot = ggplot() +
  
  geom_sf(data = x.gadm.1 %>% st_as_sfc(), fill = "white") +
  
  geom_sf(data = x.data.sf.TSS,
          
          aes( colour = loss_fac  , size = 1)) + 
  
  guides(size = "none") + 
  
  guides(colour = guide_legend(override.aes = list(size=5))) +
  
  scale_colour_viridis_d(direction=-1)   + 
  
  labs(colour = "loss") + 
  
  theme(axis.text=element_text(size=40), 
        legend.text=element_text(size=40),
        legend.title=element_text(size=40),
        legend.key.size = unit(1, "cm")) 

x.plot



ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_loss_ref.png", 
       
       plot = x.plot,
       
       units = "px",
       
       width = x.width,
       
       height = x.height)




####################################
##### ISOE_LAC  status common  #####
####################################


x.width = 2400 * 1.5

x.height = 2400 * 1.5

x.data = read.delim("Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_predictions.txt")
head(x.data)
x.data.sf = x.data %>% st_as_sf(coords = c("X_3035", "Y_3035")) %>% st_set_crs(3035)


### Base maps


st_layers("Data/gadm41_NOR.gpkg")

x.gadm.no = read_sf("Data/gadm41_NOR.gpkg", layer = "ADM_ADM_0") 
x.gadm.se = read_sf("Data/gadm41_SWE.gpkg", layer = "ADM_ADM_0")
x.gadm.fi = read_sf("Data/gadm41_FIN.gpkg", layer = "ADM_ADM_0")

x.gadm = rbind(x.gadm.no, x.gadm.se, x.gadm.fi)

x.gadm.1 = x.gadm %>%
  st_transform(3035)

### ISOE_LAC fitted probability of occurrence filtered by TSS (binary)

x.data.sf.TSS  = x.data.sf %>%
  arrange(common) %>% 
  mutate(common_fac = as.factor(common))

x.plot = ggplot() +
  
  geom_sf(data = x.gadm.1 %>% st_as_sfc(), fill = "white") +
  
  geom_sf(data = x.data.sf.TSS,
          
          aes( colour = common_fac  , size = 1)) + 
  
  guides(size = "none") + 
  
  guides(colour = guide_legend(override.aes = list(size=5))) +
  
  scale_colour_viridis_d(direction=-1)   + 
  
  labs(colour = "common") + 
  
  theme(axis.text=element_text(size=40), 
        legend.text=element_text(size=40),
        legend.title=element_text(size=40),
        legend.key.size = unit(1, "cm")) 

x.plot



ggsave(filename = "Data analyses/ISOE_LAC/glm/ISOE_LAC_GLM_TP_common_ref.png", 
       
       plot = x.plot,
       
       units = "px",
       
       width = x.width,
       
       height = x.height)




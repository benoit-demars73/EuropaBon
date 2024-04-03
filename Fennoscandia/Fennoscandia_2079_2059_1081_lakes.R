
# -------------------------------------------------------------------------
# Housekeeping

rm(list=ls()) # remove everything currently held in the R memory

graphics.off() # close all open graphic windows

# -------------------------------------------------------------------------

# code written in June 2023, updated 5th October 2023

library(sf)

setwd("~/Fennoscandia")

df <- read.table('species_climate_nosefi_v3.txt',header=TRUE)
nrow(df) # 2079 lakes
head(df) 
str(df)
df$Country <- as.factor(df$Country)

# some remaining duplicates which need deleting
# as well as some neighbour lakes with the same coordinates
# AND_426 and AND_427 should be merged, whole database should be reviewed
# recall Moran.I for autocorrelation analyses cannot have distance
# between sites =0

df %>% filter(!grepl('MD_MVM_25663', EuropaBon_ID) & 
                !grepl('MD_MVM_45505', EuropaBon_ID) & 
                !grepl('MD_MVM_29198', EuropaBon_ID) & 
                !grepl('MD_MVM_29199', EuropaBon_ID) & 
                !grepl('NOR_1058_1389', EuropaBon_ID) & 
                !grepl('NOR_1131_609', EuropaBon_ID) & 
                !grepl('NOR_1258_1375', EuropaBon_ID) & 
                !grepl('NOR_153630_443', EuropaBon_ID) & 
                !grepl('NOR_NA_823', EuropaBon_ID) & 
                !grepl('NOR_5013_826', EuropaBon_ID) & 
                !grepl('NOR_557_731', EuropaBon_ID) & 
                !grepl('NOR_5771-1227', EuropaBon_ID) & 
                !grepl('NOR_158709_1210', EuropaBon_ID) & 
                !grepl('NOR_4910_854', EuropaBon_ID) & 
                !grepl('NOR_196447_857', EuropaBon_ID) & 
                !grepl('NOR_1994_380', EuropaBon_ID) & 
                !grepl('NOR_2181_235', EuropaBon_ID) &  
                !grepl('NOR_31186_1007', EuropaBon_ID) & 
                !grepl('NOR_4517_695', EuropaBon_ID) & 
                !grepl('NOR_4517_696', EuropaBon_ID) & 
                !grepl('NOR_4527_853', EuropaBon_ID) 
) -> df1
nrow(df1) # 2059
write.table(df1,"Fennoscandia_2059_lakes.txt",sep="\t",row.names=FALSE)




##################################
#######  subset lake data    #####
##################################

# remove rows with missing data for key predictive variables
df1 %>% 
  filter(!(is.na(Area_km2) | is.na(Elevation) | 
             is.na(Alk) | is.na(Colour) | 
             is.na(TN)| is.na(TP)| is.na(deg.sum) | is.na(rr.av))) -> df2
nrow(df2) # 1104 lakes
head(df2)
# remove rows with unknown year of survey 
df2 %>% 
  filter(!(is.na(Year))) -> df3
nrow(df3) # 1104 lakes
# remove years prior 1980
df3 %>%
  filter(Year > 1980) -> df4
nrow(df4) # 1084

# remove duplicate rows in data frame
df4 %>%
  distinct()-> df5
nrow(df5) # 1084, no duplicates found 

# option, remove x y duplicates
df5 %>% distinct(X_3035, Y_3035, .keep_all = TRUE) -> df6
nrow(df6) # 1081

write.table(df6,"Fennoscandia_1081lakes.txt",sep="\t",row.names=FALSE)

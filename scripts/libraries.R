##########################################################
##   Automatizing phenomic research: ALPHA3D pipeline   ##
##                                                      ##
##                     LIBRARIES                        ## 
##                                                      ##
##                  Irene Zanandrea                     ##
##########################################################

# set working directory
PATH = "[PATH_TO_MAIN_FOLDER]"
setwd(PATH)

# install packages
#install.packages("tidyverse")
#install.packages("jsonlite")
#install.packages("evolqg")
#install.packages("lmerTest")

# load packages
library(dplyr)
library(jsonlite)
library(geomorph)  
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(writexl)
library(tidyr)
library(tidyverse)
library(scales) 
library(evolqg)

#for ANOVAS
library(lmerTest)
library(abind)
conflicted::conflicts_prefer(lmerTest::lmer)
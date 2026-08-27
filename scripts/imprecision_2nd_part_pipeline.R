##########################################################
##   Automatizing phenomic research: ALPHA3D pipeline   ##
##                                                      ##
##  TESTING IMPRECISION IN THE OTHER PART OF PIPELINE   ##      
##                                                      ##
##                  Irene Zanandrea                     ##
##########################################################


## LOAD packages, functions, landmarks, lists and dataframe
source("./scripts/libraries.R")
source("./scripts/functions_accuracy.R")
source("./scripts/landmarks_traits.R")
source("./scripts/list_nested_anova.R")



#-----------------------------------------------
# IMPRECISION IN THE 2ND HALF OF ALPHA3D
#-----------------------------------------------

######  IMPRECISION FOTO2  TS1 & TS2   --------

log_imprecision_foto2 <- calculate_imprecision(
  data1 = log_foto2_ts1_distances_df, 
  data2 = log_foto2_ts2_distances_df, 
  id_col = "ID", 
  new_suffix = "-foto2") 

log_imprecision_foto2_long <- log_imprecision_foto2 %>%  
  pivot_longer(cols = starts_with("imprecision_"),  
               names_to = "Trait",  
               values_to = "Variance_foto2")  

log_imprecision_foto2_long$Trait <- gsub("_mean", "", 
                                         log_imprecision_foto2_long$Trait)
log_imprecision_foto2_long$Trait <- gsub("^imprecision_", "", 
                                         log_imprecision_foto2_long$Trait)

log_imprecision_foto2_long <- log_imprecision_foto2_long %>%
  mutate(Trait = factor(Trait, levels = order_traits))

ggplot(log_imprecision_foto2_long, aes(x = factor(Trait), y = Variance_foto2)) +  
  geom_boxplot() +  
  stat_summary(fun = median, geom = "point", shape = 20, size = 3, 
               color = "orange", fill = "orange") +  
  labs(title = " Foto2 different TS Imprecision", 
       x = "Linear Measurament", 
       y = "Variance (%)" )+  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  


######  IMPRECISION FOTO2 BIS TS1 & TS2   --------

log_imprecision_foto2bis <- calculate_imprecision(
  data1 = log_foto2bis_ts1_distances_df, 
  data2 = log_foto2bis_ts2_distances_df, 
  id_col = "ID", 
  new_suffix = "-foto2bis") 

log_imprecision_foto2bis_long <- log_imprecision_foto2bis %>%  
  pivot_longer(cols = starts_with("imprecision_"),  
               names_to = "Trait",  
               values_to = "Variance_foto2bis")  

log_imprecision_foto2bis_long$Trait <- gsub("_mean", "", 
                                            log_imprecision_foto2bis_long$Trait)
log_imprecision_foto2bis_long$Trait <- gsub("^imprecision_", "", 
                                            log_imprecision_foto2bis_long$Trait)

log_imprecision_foto2bis_long <- log_imprecision_foto2bis_long %>%
  mutate(Trait = factor(Trait, levels = order_traits))

ggplot(log_imprecision_foto2bis_long, aes(x = factor(Trait), y = Variance_foto2bis)) +  
  geom_boxplot() +  
  stat_summary(fun = median, geom = "point", shape = 20, size = 3, 
               color = "violet", fill = "violet") +  
  labs(title = " Foto2 model 2 different TS Imprecision", 
       x = "Linear Measurament", 
       y = "Variance (%)")+   
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  




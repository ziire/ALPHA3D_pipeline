##########################################################
##   Automatizing phenomic research: ALPHA3D pipeline   ##
##                                                      ##
##               IMPRECISION AS VARIANCE                ## 
##                                                      ##
##                  Irene Zanandrea                     ##
##########################################################


## load packages, functions and lists and dataframe
source("./scripts/libraries.R")
source("./scripts/functions_accuracy.R")
source("./scripts/landmarks_traits.R")
source("./scripts/create_lists_and_dataframe.R")


# Create a new column for ordinal numbers from 1 to n
log_foto2_distances_df_long$TraitOrdinal <- as.numeric(log_foto2_distances_df_long$Trait)
log_foto2_distances_df_long$TraitOrdinal <- factor(log_foto2_distances_df_long$TraitOrdinal)


#---------------------------
#  TOTAL IMPRECISION
#---------------------------

#Calculate imprecision (function defined in functions_accuracy.R)
#Gives back df -> I call it log_imprecision_total
log_imprecision_total <- calculate_imprecision(
            data1 = log_foto1_distances_df, 
            data2 = log_foto2_distances_df, 
            id_col = "ID", 
            new_suffix = "-total") 
        
#Make the dataset in the long format to plot it
log_imprecision_total_long <- log_imprecision_total %>%  
    pivot_longer(cols = starts_with("imprecision_"),  
                 names_to = "Trait",  
                 values_to = "Variance_total")

#Remove parts not needed
log_imprecision_total_long$Trait <- gsub("_mean", "", log_imprecision_total_long$Trait)
log_imprecision_total_long$Trait <- gsub("^imprecision_", "", log_imprecision_total_long$Trait)
  
#Reorder the traits (order_traits defined in landmarks_and_traits.R)
  log_imprecision_total_long <- log_imprecision_total_long %>%
    mutate(Trait = factor(Trait,levels = order_traits))
  
#Plot
plot_1 <- ggplot(
  log_imprecision_total_long, 
  aes(x = factor(Trait), y = Variance_total)) +  
  geom_boxplot() +  
  stat_summary(fun = median, geom = "point", shape = 20, size = 3, 
               color = "mediumvioletred", fill = "mediumvioletred") +  
  labs(title = "TOTAL Imprecision", x = "Linear Measurament", 
       y = "Logged Variance (%)") +  
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
plot_1
  
#Error   
  log_tot_error_index <- which(log_imprecision_total_long$Variance_total >= 0.04)
  log_tot_error_index_values <- log_imprecision_total_long[log_tot_error_index, ]
  print(log_tot_error_index_values)
  

#---------------------------------------
# IMPRECISION CAUSED BY 3D AND AL
#---------------------------------------
  
log_imprecision_3D_AL <- calculate_imprecision(
                           data1 = log_foto1_distances_df, 
                           data2 = log_fotobis_distances_df, 
                           id_col = "ID", 
                           new_suffix = "-3D_AL") 

log_imprecision_3D_AL_long <- log_imprecision_3D_AL %>%  
              pivot_longer(cols = starts_with("imprecision_"),  
                           names_to = "Trait",  
                           values_to = "Variance_3D_AL")  

log_imprecision_3D_AL_long$Trait <- gsub("_mean", "",log_imprecision_3D_AL_long$Trait)
log_imprecision_3D_AL_long$Trait <- gsub("^imprecision_", "",log_imprecision_3D_AL_long$Trait)

log_imprecision_3D_AL_long <- log_imprecision_3D_AL_long %>%
                           mutate(Trait = factor(Trait, levels = order_traits))

plot_2 <- ggplot(
  log_imprecision_3D_AL_long, 
  aes(x = factor(Trait), y = Variance_3D_AL)) + 
  geom_boxplot() +  
  stat_summary(fun = median, geom = "point", shape = 20, size = 3, 
               color = "darkgoldenrod", fill = "darkgoldenrod") +  
  labs(title = "3D Modeling and AL Imprecision", 
       x = "Linear Measurament", 
       y = "Logged Variance (%)") +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
plot_2

error_index_3D_AL <- which(log_imprecision_3D_AL_long$Variance_3D_AL >= 0.06)
error_index_values_3D_AL <- log_imprecision_3D_AL_long[error_index_3D_AL, ]
print(error_index_values_3D_AL)


#-----------------------------------
# IMPRECISION CAUSED ONLY BY AL
#-----------------------------------

log_imprecision_AL <- calculate_imprecision(
                data1 = log_foto1_distances_df, 
                data2 = log_fotoextra_distances_df, 
                id_col = "ID", 
                new_suffix = "-AL") 

log_imprecision_AL_long <- log_imprecision_AL %>%  
             pivot_longer(cols = starts_with("imprecision_"),  
                          names_to = "Trait",  
                          values_to = "Variance_AL")  

log_imprecision_AL_long$Trait <- gsub("_mean", "", 
                                      log_imprecision_AL_long$Trait)
log_imprecision_AL_long$Trait <- gsub("^imprecision_", "", 
                                  log_imprecision_AL_long$Trait)

log_imprecision_AL_long <- log_imprecision_AL_long %>%
                        mutate(Trait = factor(Trait, levels = order_traits))

plot_3 <- ggplot(
  log_imprecision_AL_long, 
  aes(x = factor(Trait), y = Variance_AL)) +  
  geom_boxplot() +  
  stat_summary(fun = median, geom = "point", shape = 20, size = 3, 
               color = "forestgreen", fill = "forestgreen") +  
  labs(title = " AL (Automated Landmarking) Imprecision", 
       x = "Linear Measurament", 
       y = "Logged Variance (%)") +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
plot_3

error_index_AL <- which(log_imprecision_AL_long$Variance_AL >= 0.03)
error_index_values_AL <- log_imprecision_AL_long[error_index_AL, ]
print(error_index_values_AL)


#-----------------------------------
# MEAN AND MEDIAN of variance
#-----------------------------------

log_mean_imprecision_total <- mean(log_imprecision_total_long$Variance_total, na.rm = TRUE)
log_mean_imprecision_3D_AL <- mean(log_imprecision_3D_AL_long$Variance_3D_AL, na.rm = TRUE)
log_mean_imprecision_AL <- mean(log_imprecision_AL_long$Variance_AL, na.rm = TRUE)

log_median_imprecision_total <- median(log_imprecision_total_long$Variance_total, na.rm = TRUE)
log_median_imprecision_3D_AL <- median(log_imprecision_3D_AL_long$Variance_3D_AL, na.rm = TRUE)
log_median_imprecision_AL <- median(log_imprecision_AL_long$Variance_AL, na.rm = TRUE)

### Table 1. Mean and median values of imprecision
log_mean_median_df  <- data.frame(
  Imprecision_total = c(log_mean_imprecision_total, log_median_imprecision_total),
  Imprecision_3D_AL = c(log_mean_imprecision_3D_AL, log_median_imprecision_3D_AL),
  Imprecision_AL = c(log_mean_imprecision_AL,log_median_imprecision_AL))
rownames(log_mean_median_df) <- c("Mean", "Median")  # Impostazione dei nomi delle righe



#---------------------------------------------
# PLOT TOGETHER ALL THE THREE IMPRECISIONS
#---------------------------------------------

##### 1 - in the 3 datasets
#Removing parts
log_imprecision_total_long$ID <- gsub("-total", "", log_imprecision_total_long$ID)
log_imprecision_3D_AL_long$ID <- gsub("-3D_AL", "", log_imprecision_3D_AL_long$ID)
log_imprecision_AL_long$ID <- gsub("-AL", "", log_imprecision_AL_long$ID)

##### 2 - Combine 
#First the first two data frames
log_combined_df_1 <- merge(log_imprecision_total_long, log_imprecision_3D_AL_long, 
                           by = c ("Trait", "ID"), all = TRUE) 
log_combined_df <- merge(log_combined_df_1, log_imprecision_AL_long, 
                         by = c ("Trait", "ID"), all = TRUE)   
log_combined_df <- log_combined_df %>%
              mutate(Trait = factor(Trait, levels = order_traits))

##### 3 - Create a summary
summary_log_combined_df <- log_combined_df %>%
  group_by(Trait) %>%
  summarise(
    Mean_Imprecision_total  = round(mean(Variance_total), 4),#arrotonda a 4 decimali
    Median_Imprecision_total= round(median(Variance_total), 4),
    Mean_Imprecision_3D_AL  = round(mean(Variance_3D_AL), 4),
    Median_Imprecision_3D_AL= round(median(Variance_3D_AL), 4),
    Mean_Imprecision_AL     = round(mean(Variance_AL), 4),
    Median_Imprecision_AL   = round(median(Variance_AL), 4))
#Write the csv
#write.csv(summary_log_combined_df, "log_mean_imprecision_all_traits.csv", row.names = FALSE)

##### 4- Plot the combined dataset
# Convert Trait to a factor
log_combined_df$Trait <- factor(log_combined_df$Trait, levels = order_traits)  
# Create a new column for ordinal numbers from 1 to n
log_combined_df$TraitOrdinal <- as.numeric(log_combined_df$Trait)
# Convert the ordinal column to a factor
log_combined_df$TraitOrdinal <- factor(log_combined_df$TraitOrdinal)
#Make the long format
log_combined_df_long <- log_combined_df %>%
           pivot_longer(cols = tidyselect::matches("Variance_"),
                        names_to = "Imprecision",
                        values_to = "Variance")

log_combined_df_long <- log_combined_df_long %>%
        dplyr::mutate(Imprecision = factor(Imprecision,
                      levels = c("Variance_total","Variance_3D_AL","Variance_AL")))
#Create ggplot
ggplot(log_combined_df_long, 
       aes(x = TraitOrdinal, y = Variance, color = Imprecision)) +
      geom_boxplot(aes(fill = Imprecision), alpha = 0.5) +
      labs(x = "Linear Measurements", y = "Logged Variance (%)",color = "Imprecision") +
      # defining manual colors and legend labels
      #1- colors for the LEGEND
      scale_fill_manual(name   = "Imprecision",       #Sets the legend title to "Imprecision"
                        breaks = c("Variance_total",  #Uses levels "Var_total", "Var_3D_AL", "Var_AL"
                                   "Variance_3D_AL", "Variance_AL"),
                        labels = c("Total", "3D and AL", "AL"), #how I want the (Total, AL)
                        values = c("Variance_total" = "mediumvioletred", #these are the colours
                                   "Variance_3D_AL" = "darkgoldenrod",
                                   "Variance_AL"    = "forestgreen")) +
      #2- colors for the POINTS in the plot. Assigns the same colors to each level for line/point colors 
      scale_colour_manual(values = c("Variance_total" = "mediumvioletred",
                                     "Variance_3D_AL" = "darkgoldenrod",
                                     "Variance_AL"    = "forestgreen")) +
      theme_minimal() +
      theme(legend.position = "bottom",
            legend.title = element_text(size = 18, face = "bold"),
            legend.text  = element_text(size = 21),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.grid.major.x = element_blank(),
            axis.text.x  = element_text(angle = 45, hjust = 1, size = 21),
            axis.text.y  = element_text(hjust = 1, size = 21),
            axis.title.x = element_text(size = 22),
            axis.title.y = element_text(size = 22)) +
      guides(color = "none")


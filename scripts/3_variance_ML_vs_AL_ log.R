##########################################################
##   Automatizing phenomic research: ALPHA3D pipeline   ##
##                                                      ##
##         AUTOMATED vs MANUAL LANDMARKING              ## 
##                                                      ##
##                  Irene Zanandrea                     ##
##########################################################


## load packages, functions and lists and dataframe
source("./scripts/libraries.R")
source("./scripts/functions_accuracy.R")
source("./scripts/landmarks_traits.R")



#-------------------------------------
# AUTOMATED LANDMARKING DATAFRAME
#-------------------------------------

######  FOTO 1 TS1   ---------------

foto1.list <- list.files(path = "./landmark_ml_vs_al/0_AL1",
                         full.names = TRUE, recursive = TRUE) 
foto1.append = "-foto1"       
foto1_landmarks_array <- lapply(foto1.list, process_json, landmarks_n = 34, coord_n = 3)
# extract BILATERAL traits
foto1_distances_bil <- c()
foto1_distances_bil_df <- data.frame()
for(j in 1:length(foto1_landmarks_array)){ #for each skull
 for(i in seq_along(landmark_pairs)){     #find each pair
   foto1_distances_bil[i] <- euclidean_distance(x = foto1_landmarks_array[[j]], 
                            point1 = landmark_pairs[[i]][1],point2 = landmark_pairs[[i]][2])}
   foto1_distances_bil_df <- rbind(foto1_distances_bil_df, foto1_distances_bil)}
ID.foto1 <- foto1.list %>%   
        basename() %>% gsub("_median\\.mrk\\.json$", "", .) %>%  paste0(., foto1.append)                  
colnames(foto1_distances_bil_df) <- traits
foto1_distances_bil_df$ID <- ID.foto1
  #avarage for bilateral
   foto1_distances_avg_df <- avg_columns(foto1_distances_bil_df)
   colnames(foto1_distances_avg_df)[1] <- "ID"          
   colnames(foto1_distances_avg_df)[-1] <- name_traits  
# extract OTHER traits
foto1_distances_other <- c()
foto1_distances_other_df <- data.frame()
for(j in 1:length(foto1_landmarks_array)){#for each skull
  for(i in seq_along(landmark_pairs1)){    # LANDMARKS_PAIRS_1(for other traits!) 
    foto1_distances_other[i] <- euclidean_distance(x = foto1_landmarks_array[[j]], 
                    point1 = landmark_pairs1[[i]][1], point2 = landmark_pairs1[[i]][2])} 
    foto1_distances_other_df <- rbind(foto1_distances_other_df, foto1_distances_other)}
ID.foto1 <- foto1.list %>%                  # Extract the ID
        basename() %>% gsub("_median\\.mrk\\.json$", "", .) %>%  paste0(., foto1.append)                   
colnames(foto1_distances_other_df) <- traits1  #TRAITS1 (for other traits!!) !!!!!!
foto1_distances_other_df$ID <- ID.foto1
# Unite the two
foto1_distances_df <- merge(foto1_distances_avg_df,foto1_distances_other_df, by = c ("ID"))
# Logtransform
log_foto1_distances_df <- cbind(foto1_distances_df[, 1, drop = FALSE],
                                log(foto1_distances_df[, 2:36]))
log_foto1_distances_df$ID <- gsub("-foto1", "", log_foto1_distances_df$ID)
 

######  FOTOEXTRA (FOTO1 TS2)    ---------------

fotoextra.list <- list.files(path = "./landmark_ml_vs_al/0_AL2",
                             pattern = "\\.json",full.names = TRUE, recursive = TRUE) 
fotoextra.append = "-fotoextra"
fotoextra_landmarks_array <- lapply(fotoextra.list, process_json,landmarks_n = 34,  coord_n = 3)
# extract BILATERAL traits 
fotoextra_distances_bil <- c()
fotoextra_distances_bil_df <- data.frame()
for(j in 1:length(fotoextra_landmarks_array)){ #for each skull
  for(i in seq_along(landmark_pairs)){ #find each pair
    fotoextra_distances_bil[i] <- euclidean_distance(x = fotoextra_landmarks_array[[j]], 
                                point1 = landmark_pairs[[i]][1], point2 = landmark_pairs[[i]][2])}
    fotoextra_distances_bil_df <- rbind(fotoextra_distances_bil_df,fotoextra_distances_bil)}
ID.fotoextra <- fotoextra.list %>%
        basename() %>% gsub("_median\\.mrk\\.json$", "", .) %>% paste0(., fotoextra.append)  
colnames(fotoextra_distances_bil_df) <- traits
fotoextra_distances_bil_df$ID <- ID.fotoextra
  #average for bilateral
   fotoextra_distances_avg_df <- avg_columns(fotoextra_distances_bil_df)
   colnames(fotoextra_distances_avg_df)[1] <- "ID"          
   colnames(fotoextra_distances_avg_df)[-1] <- name_traits  
# extract OTHER traits
fotoextra_distances_other <- c()
fotoextra_distances_other_df <- data.frame()
for(j in 1:length(fotoextra_landmarks_array)){ #for each skull
  for(i in seq_along(landmark_pairs1)){ # LANDMARKS_PAIRS_1 (for other traits!!) 
    fotoextra_distances_other[i] <- euclidean_distance(x = fotoextra_landmarks_array[[j]], 
                                   point1 = landmark_pairs1[[i]][1], point2 = landmark_pairs1[[i]][2])}
    fotoextra_distances_other_df <- rbind(fotoextra_distances_other_df,fotoextra_distances_other)}
ID.fotoextra <- fotoextra.list %>%
      basename() %>%  gsub("_median\\.mrk\\.json$", "", .) %>%  paste0(., fotoextra.append)  
colnames(fotoextra_distances_other_df) <- traits1 # TRAITS1 (for other traits!!)
fotoextra_distances_other_df$ID <- ID.fotoextra
# Unite the two
fotoextra_distances_df <- merge(fotoextra_distances_avg_df,fotoextra_distances_other_df, by = c ("ID"))
# Logtransform
log_fotoextra_distances_df <- cbind(fotoextra_distances_df[, 1, drop = FALSE],
                                log(fotoextra_distances_df[, 2:36]))
log_fotoextra_distances_df$ID <- gsub("-fotoextra", "", log_fotoextra_distances_df$ID)
 


#----------------------------------------
# MANUAL LANDMARKING DATAFRAME
#----------------------------------------

######  ML 1st ROUND ML    ------------

manual1.list <- list.files(path = "./landmark_ml_vs_al/1_ML1",
                           full.names = TRUE,recursive = TRUE) 
manual1.append = "-manual1" 
manual1_landmarks_array <- lapply(manual1.list, process_json, 
                                  landmarks_n = 34,coord_n = 3)
### extract BILATERAL traits
manual1_distances_bil <- c()
manual1_distances_bil_df <- data.frame()
for(j in 1:length(manual1_landmarks_array)){ #for each skull
  for(i in seq_along(landmark_pairs)){ #find each pair
    manual1_distances_bil[i] <- euclidean_distance(x = manual1_landmarks_array[[j]], 
                                point1 = landmark_pairs[[i]][1],point2 = landmark_pairs[[i]][2])}
    manual1_distances_bil_df <- rbind(manual1_distances_bil_df, manual1_distances_bil)}
ID.manual1 <- manual1.list %>%      # Extract the ID
        basename() %>% gsub(".mrk\\.json$", "", .) %>% paste0(., manual1.append)  
colnames(manual1_distances_bil_df) <- traits
manual1_distances_bil_df$ID <- ID.manual1
 # average for bilateral
 manual1_distances_avg_df <- avg_columns(manual1_distances_bil_df)
 colnames(manual1_distances_avg_df)[1] <- "ID"         
 colnames(manual1_distances_avg_df)[-1] <- name_traits  
### extract OTHER traits 
manual1_distances_other <- c()
manual1_distances_other_df <- data.frame()
for(j in 1:length(manual1_landmarks_array)){ #for each skull
  for(i in seq_along(landmark_pairs1)){ #find each pair
    manual1_distances_other[i] <- euclidean_distance(x = manual1_landmarks_array[[j]], 
                               point1 = landmark_pairs1[[i]][1], point2 = landmark_pairs1[[i]][2])}
    manual1_distances_other_df <- rbind(manual1_distances_other_df, manual1_distances_other)}
ID.manual1 <- manual1.list %>%         # Extract the ID
        basename() %>%  gsub(".mrk\\.json$", "", .) %>% paste0(., manual1.append)        
colnames(manual1_distances_other_df) <- traits1 #TRAITS 1!!!
manual1_distances_other_df$ID <- ID.manual1
### Unite the two 
manual1_distances_df <- merge(manual1_distances_avg_df,manual1_distances_other_df, by = c ("ID"))
# Logtransform
log_manual1_distances_df <- cbind(manual1_distances_df[, 1, drop = FALSE],
                                    log(manual1_distances_df[, 2:36]))
log_manual1_distances_df$ID <- gsub("-manual1", "", log_manual1_distances_df$ID)


######  ML 2nd ROUND ML    ------------

manual2.list <- list.files(path = "./landmark_ml_vs_al/1_ML2",
                           full.names = TRUE, recursive = TRUE) 
manual2.append = "-manual2"
manual2_landmarks_array <- lapply(manual2.list, process_json,landmarks_n = 34, coord_n = 3)

### extract BILATERAL traits
manual2_distances_bil <- c()
manual2_distances_bil_df <- data.frame()
for(j in 1:length(manual2_landmarks_array)){ #for each skull
  for(i in seq_along(landmark_pairs)){ #find each pair
    manual2_distances_bil[i] <- euclidean_distance(x = manual2_landmarks_array[[j]], 
                          point1 = landmark_pairs[[i]][1],point2 = landmark_pairs[[i]][2])}
    manual2_distances_bil_df <- rbind(manual2_distances_bil_df, manual2_distances_bil)}
ID.manual2 <- manual2.list %>%                  # Extract the ID
          basename() %>%  gsub(".mrk\\.json$", "", .) %>%  paste0(., manual2.append)                  
colnames(manual2_distances_bil_df) <- traits
manual2_distances_bil_df$ID <- ID.manual2
 # average for bilateral
   manual2_distances_avg_df <- avg_columns(manual2_distances_bil_df)
   colnames(manual2_distances_avg_df)[1] <- "ID"          
   colnames(manual2_distances_avg_df)[-1] <- name_traits  
### extract OTHER traits
manual2_distances_other <- c()
manual2_distances_other_df <- data.frame()
for(j in 1:length(manual2_landmarks_array)){ #for each skull
  for(i in seq_along(landmark_pairs1)){      #find each pair
   manual2_distances_other[i] <- euclidean_distance(x = manual2_landmarks_array[[j]], 
                                  point1 = landmark_pairs1[[i]][1],point2 = landmark_pairs1[[i]][2])}
   manual2_distances_other_df <- rbind(manual2_distances_other_df, manual2_distances_other)}
ID.manual2 <- manual2.list %>%                 
       basename() %>% gsub(".mrk\\.json$", "", .) %>%  paste0(., manual2.append)                   
colnames(manual2_distances_other_df) <- traits1
manual2_distances_other_df$ID <- ID.manual2
### Unite the two
manual2_distances_df <- merge(manual2_distances_avg_df, manual2_distances_other_df, by = c ("ID"))
# Logtransform
log_manual2_distances_df <- cbind(manual2_distances_df[, 1, drop = FALSE],
                                  log(manual2_distances_df[, 2:36]))
log_manual2_distances_df$ID <- gsub("-manual2", "", log_manual2_distances_df$ID)



#----------------------------------------
# CALCULATE VARIANCE
#----------------------------------------

####  AUTOMATED LANDMARKING  -----
log_imprecision_AL <- calculate_imprecision(data1 = log_foto1_distances_df, 
                          data2 = log_fotoextra_distances_df, 
                          id_col = "ID",
                          new_suffix = "-AL") 
log_imprecision_AL_long <- log_imprecision_AL %>%  
  pivot_longer(cols = starts_with("imprecision_"),  
               names_to = "Trait",  values_to = "Variance_AL") 

####  MANUAL LANDMARKING  -----
log_imprecision_ML <- calculate_imprecision(data1 = log_manual1_distances_df, 
                            data2 = log_manual2_distances_df, 
                            id_col = "ID", 
                            new_suffix = "-ML") 
log_imprecision_ML_long <- log_imprecision_ML %>%  
  pivot_longer(cols = starts_with("imprecision_"),  
               names_to = "Trait",  values_to = "Variance_ML")  


#----------------------------------------
# PREPARE DATASETS 
#----------------------------------------

## Removing parts
log_imprecision_AL_long$ID <- gsub("-AL", "", log_imprecision_AL_long$ID)
log_imprecision_AL_long$Trait <- gsub("_mean", "", log_imprecision_AL_long$Trait)
log_imprecision_AL_long$Trait <- gsub("^imprecision_", "", log_imprecision_AL_long$Trait)

log_imprecision_ML_long$ID <- gsub("-ML", "", log_imprecision_ML_long$ID)
log_imprecision_ML_long$Trait <- gsub("_mean", "", log_imprecision_ML_long$Trait)
log_imprecision_ML_long$Trait <- gsub("^imprecision_", "", log_imprecision_ML_long$Trait)

log_imprecision_ML_long <- log_imprecision_ML_long %>%
                        mutate(Trait = factor(Trait, levels = order_traits))
# Plot MANUAL
ggplot(log_imprecision_ML_long, aes(x = factor(Trait),y = Variance_ML)) +  
      geom_boxplot() +  
      stat_summary(fun = median, geom = "point", shape = 20, size = 3, 
                   color = "deepskyblue4", fill = "deepskyblue4") +  
      labs(title = " ML (Manual Landmarking) Imprecision", 
           x = "Linear Measurament", y = "Variance (%)") +  
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# Plot AUTOMATIC
ggplot(log_imprecision_AL_long, aes(x = factor(Trait), y = Variance_AL)) +  
      geom_boxplot() +  
      stat_summary(fun = median, geom = "point", shape = 20, size = 3, 
                   color = "forestgreen", fill = "forestgreen") +  
      labs(title = " AL (Automatic Landmarking) Imprecision", 
           x = "Linear Measurament", 
           y = "Variance (%)")
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
    
    
    
#------------------------------------------------------
# Calculate the MEAN and MEDIAN for both datasets 
#------------------------------------------------------    

mean_imprecision_AL_long<- mean(log_imprecision_AL_long$Variance_AL, na.rm = TRUE)
mean_imprecision_ML_long <- mean(log_imprecision_ML_long$Variance_ML, na.rm = TRUE)
median_imprecision_AL_long <- median(log_imprecision_AL_long$Variance_AL, na.rm = TRUE)
median_imprecision_ML_long<- median(log_imprecision_ML_long$Variance_ML, na.rm = TRUE)

## Combine and long format
imprecision_AL_ML <- merge(log_imprecision_AL_long, log_imprecision_ML_long, 
                     by = c ("Trait", "ID"), all = TRUE)  
#New names and numbers 
imprecision_AL_ML$Trait <- factor(imprecision_AL_ML$Trait,levels = order_traits) 


##### MEAN PER TRAIT ------

## Table S5: Comparison of MEAN for each trait. 
mean_imprecision_AL <- imprecision_AL_ML %>%  # dataframe mean imprecision AL
              group_by(Trait) %>%
              summarise(Mean_Variance_AL = round(mean(Variance_AL, na.rm = TRUE),4))

mean_imprecision_ML <- imprecision_AL_ML %>%  # dataframe mean imprecision AL
  group_by(Trait) %>%
  summarise(Mean_Variance_ML = round(mean(Variance_ML, na.rm = TRUE),4))
mean_imprecision_ML_AL <- merge(mean_imprecision_AL, mean_imprecision_ML, by = c ("Trait"), all = TRUE)
mean_imprecision_ML_AL$Trait <- factor(mean_imprecision_ML_AL$Trait, 
                                       levels = order_traits)  # new order for the traits
mean_imprecision_ML_AL <- mean_imprecision_ML_AL %>%  arrange(Trait)
mean_imprecision_ML_AL$Trait <- factor(mean_imprecision_ML_AL$Trait,   levels = order_traits)
# Save  dataset mean_ML_AL in a CSV file 
#write.csv(mean_imprecision_ML_AL, "mean_variance_log_AL_ML.csv", row.names = FALSE)

## Table S6: Comparison of MEDIAN  for each trait. 
median_imprecision_AL <- imprecision_AL_ML %>%  # dataframe mean imprecision AL
  group_by(Trait) %>%
  summarise(Median_Variance_AL = round(median(Variance_AL, na.rm = TRUE),4))

median_imprecision_ML <- imprecision_AL_ML %>%  # dataframe mean imprecision AL
  group_by(Trait) %>%
  summarise(Median_Variance_ML = round(median(Variance_ML, na.rm = TRUE),4))

median_imprecision_ML_AL <- merge(median_imprecision_AL, 
                                  median_imprecision_ML, by = c ("Trait"), all = TRUE)
median_imprecision_ML_AL$Trait <- factor(median_imprecision_ML_AL$Trait, 
                                         levels = order_traits)  # Imposta di nuovo l'ordine dei traits
median_imprecision_ML_AL <- median_imprecision_ML_AL %>% # Ordina il dataframe f
  arrange(Trait)
median_imprecision_ML_AL$Trait <- factor(median_imprecision_ML_AL$Trait,levels = order_traits)
# Save  dataset median_imprecision_ML_AL 
#write.csv(median_imprecision_ML_AL, "median_variance_log_AL_ML.csv", row.names = FALSE)



#----------------------------------------
# PLOT TOGETHER 
#----------------------------------------

imprecision_long <- imprecision_AL_ML %>%
  pivot_longer(cols = tidyselect::matches("^Variance_"),
               names_to = "Imprecision",
               values_to = "Variance")
# Create a new column for ordinal numbers from 1 to n
imprecision_long$TraitOrdinal <- as.numeric(imprecision_long$Trait)
# Convert the ordinal column to a factor
imprecision_long$TraitOrdinal <- factor(imprecision_long$TraitOrdinal)

#ggplot
ggplot(imprecision_long, aes(x = TraitOrdinal, y = Variance, color = Imprecision)) +
      geom_boxplot(aes(fill = Imprecision), alpha = 0.5) +
      labs(title = "Comparison of imprecision between manual and automated landmarking", 
           x = "Linear Measurements", 
           y = "Logged Variance (%)", 
           color = "Imprecision") +
      scale_fill_manual(name   = "Imprecision",
                        breaks = c("Variance_AL","Variance_ML"),
                        labels = c("AL", "ML"),
                        values = c("Variance_AL"= "forestgreen",
                                   "Variance_ML"= "deepskyblue4")) +
      scale_colour_manual(values = c("Variance_AL"= "forestgreen",
                                     "Variance_ML"= "deepskyblue4")) +
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

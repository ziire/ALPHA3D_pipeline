##########################################################
##   Automatizing phenomic research: ALPHA3D pipeline   ##
##                                                      ##
##                       BIAS                           ## 
##                                                      ##
##                  Irene Zanandrea                     ##
##########################################################


## load packages, functions and lists and dataframe
source("./scripts/libraries.R")
source("./scripts/functions_accuracy.R")
source("./scripts/landmarks_traits.R")

#---------------------------
#  MANUAL LANDMARKS
#---------------------------

manual.list <- list.files(path = "./landmark_bias/1_ML",
                          pattern = "\\.json",full.names = TRUE, recursive = TRUE) 
manual.list
manual.append = "-ml"
manual_landmarks_array <- lapply(manual.list, process_json, landmarks_n = 34, coord_n = 3)

# extract BILATERAL traits
  manual_distances_bil <- c()                     
  manual_distances_bil_df <- data.frame()
  for(j in 1:length(manual_landmarks_array)){ #for each skull
    for(i in seq_along(landmark_pairs)){      #find each pair
      manual_distances_bil[i] <- euclidean_distance(x = manual_landmarks_array[[j]], 
                            point1 = landmark_pairs[[i]][1],point2 = landmark_pairs[[i]][2])}
  manual_distances_bil_df <- rbind(manual_distances_bil_df,manual_distances_bil)}
  ID.manual <- manual.list %>%
     basename() %>% gsub("_median\\.mrk\\.json$", "", .) %>% paste0(., manual.append) 
  colnames(manual_distances_bil_df) <- traits
  manual_distances_bil_df$ID <- ID.manual
  #average
  manual_distances_avg_df <- avg_columns(manual_distances_bil_df)
  colnames(manual_distances_avg_df)[1] <- "ID"         
  colnames(manual_distances_avg_df)[-1] <- name_traits  
# extract OTHER traits
  manual_distances_other <- c()                   
  manual_distances_other_df <- data.frame()
  for(j in 1:length(manual_landmarks_array)){ #for each skull
    for(i in seq_along(landmark_pairs1)){      #find each pair
      manual_distances_other[i] <- euclidean_distance(x = manual_landmarks_array[[j]], 
                         point1 = landmark_pairs1[[i]][1],point2 = landmark_pairs1[[i]][2])}
      manual_distances_other_df <- rbind(manual_distances_other_df,manual_distances_other)}
  ID.manual <- manual.list %>%
        basename() %>% gsub("_median\\.mrk\\.json$", "", .) %>%  paste0(., manual.append)           
  colnames(manual_distances_other_df) <- traits1
  manual_distances_other_df$ID <- ID.manual
# Unite the two
  manual_distances_df <- merge(manual_distances_avg_df, manual_distances_other_df, by = c ("ID"))
# Logtransform
  log_manual_distances_df <- cbind(manual_distances_df[, 1, drop = FALSE],
                                   log(manual_distances_df[, 2:36]))
  log_manual_distances_df$ID <- gsub(".mrk.json-ml", "", log_manual_distances_df$ID)
  
      
#---------------------------
#  AUTOMATED LANDMARKS
#---------------------------

auto.list <- list.files(path = "./landmark_bias/2_AL",
                        pattern = "\\.json", full.names = TRUE,recursive = TRUE) 
auto.list
auto.append = "-al"
auto_landmarks_array <- lapply(auto.list, process_json, landmarks_n = 34, coord_n = 3)

# extract BILATERAL traits
  auto_distances_bil <- c()
  auto_distances_bil_df <- data.frame()
  for(j in 1:length(auto_landmarks_array)){  #for each skull
    for(i in seq_along(landmark_pairs)){     #find each pair
      auto_distances_bil[i] <- euclidean_distance(x = auto_landmarks_array[[j]], 
                                point1 = landmark_pairs[[i]][1],point2 = landmark_pairs[[i]][2])}
      auto_distances_bil_df <- rbind(auto_distances_bil_df, auto_distances_bil)}
  ID.auto <- auto.list %>%
     basename() %>% gsub("_median\\.mrk\\.json$", "", .) %>% paste0(., auto.append)     
  colnames(auto_distances_bil_df) <- traits
  auto_distances_bil_df$ID <- ID.auto
  #average
  auto_distances_avg_df <- avg_columns(auto_distances_bil_df)
  colnames(auto_distances_avg_df)[1] <- "ID"        
  colnames(auto_distances_avg_df)[-1] <- name_traits  
# extract OTHER traits
  auto_distances_other <- c()
  auto_distances_other_df <- data.frame()
  for(j in 1:length(auto_landmarks_array)){  
    for(i in seq_along(landmark_pairs1)){      
     auto_distances_other[i] <- euclidean_distance(x = auto_landmarks_array[[j]], 
                            point1 = landmark_pairs1[[i]][1],point2 = landmark_pairs1[[i]][2])}
     auto_distances_other_df <- rbind(auto_distances_other_df, auto_distances_other)}
  ID.auto <- auto.list %>%
        basename() %>%  gsub("_median\\.mrk\\.json$", "", .) %>% paste0(., auto.append)   
  colnames(auto_distances_other_df) <- traits1
  auto_distances_other_df$ID <- ID.auto
# Unite the two
  auto_distances_df <- merge(auto_distances_avg_df, auto_distances_other_df, by = c ("ID"))
# Logtransform    
  log_auto_distances_df <- cbind(auto_distances_df[, 1, drop = FALSE],
                                     log(auto_distances_df[, 2:36]))
  log_auto_distances_df$ID <- gsub("-al", "", log_auto_distances_df$ID)



  
  
#---------------------------
#   BIAS COMPUTATION
#---------------------------

### 1- calculate differences for traits for each cranium between auto and manual
diff_auto_man <- bias_function(log_auto_distances_df, log_manual_distances_df)
names(diff_auto_man) <- gsub("_mean$", "", names(diff_auto_man))#remove "_mean" from columns
### 2- long format
diff_long <- diff_auto_man %>%  
              pivot_longer(cols = -ID,  names_to = "Trait",  
                           values_to = "Difference")  
diff_long$Trait <- gsub("_mean", "", diff_long$Trait)

### Table S.5 Mean bias analysis for each trait on a subset of 35 skulls. 
column_means_1 <- round(
  colMeans(diff_auto_man[ , -which(names(diff_auto_man) == "ID")],
           na.rm = TRUE), 3)
means_bias <- data.frame(Trait = names(column_means_1), 
                         Mean = column_means_1)
means_bias$Trait <- factor(means_bias$Trait,levels = order_traits)
#write_csv(means_bias, "bias_log.csv")

### Table S.7 Mean bias SQUARED
column_means_squared <- round(colMeans(
                          diff_auto_man[ , -which(names(diff_auto_man) == "ID")],
                          na.rm = TRUE)^2, 3) #SQUARED
means_bias_squared <- data.frame(Trait = names(column_means_squared), 
                         Mean = column_means_squared)
means_bias_squared$Trait <- factor(means_bias_squared$Trait, 
                           levels = order_traits)
write_csv(means_bias_squared, "bias_squared_log.csv")

### 3-  plot
diff_long$Trait <- factor(diff_long$Trait, 
                          levels = order_traits)# Convert Trait to a factor with levels
#Convert traits in number
diff_long$TraitOrdinal <- as.numeric(diff_long$Trait)    # new column for ordinal numbers 
diff_long$TraitOrdinal <- factor(diff_long$TraitOrdinal) # Convert ordinal column to a factor
#ggplot 
ggplot(diff_long, aes(x = TraitOrdinal, y = Difference)) +
      geom_boxplot() +
      labs(x = "Linear Measurements",y = "Logged Bias (%)") +
      theme_minimal() + 
      theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.grid.major.x = element_blank(),  # Remove major vertical grid lines
            axis.text.x  = element_text(angle = 45, hjust = 1, size = 21),
            axis.text.y  = element_text(hjust = 1, size = 21),
            axis.title.x = element_text(size = 22),
            axis.title.y = element_text(size = 22)) +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1.9)
  
  

#-------------------------------------------
#   BIAS SQUARED for INACCURACY ANALYSIS
#-------------------------------------------  
  
bias_long <- diff_long
bias_long$ID <- gsub("-al", "", bias_long$ID)
bias_long <- subset(bias_long,select = -c(TraitOrdinal))
bias_long$Bias_2 <- bias_long$Difference^2

ggplot(bias_long, aes(x = factor(Trait), y = Bias_2)) +
  geom_boxplot() +
  labs(title = "Bias squared of automated landmarking relative to manual landmarking", 
       x = "Linear Measurements", 
       y = expression(paste("Bias squared (cm"^2, ")"))) + 
  theme_minimal() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid.major.x = element_blank(),  # Remove major vertical grid lines
        axis.text.x = element_text(angle = 45, hjust = 1, size = 27),
        axis.text.y = element_text(hjust = 1, size = 25),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(size = 23),
        plot.title = element_text(size = 25, face = "bold")) 

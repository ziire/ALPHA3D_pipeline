##########################################################
##   Automatizing phenomic research: ALPHA3D pipeline   ##
##                                                      ##
##                 DATA FOR ANOVAS                      ## 
##          CREATE LISTS and DATAFRAMES                 ##
##                                                      ##
##                  Irene Zanandrea                     ##
##########################################################


## load packages, functions and landmark names
source("./scripts/libraries.R")
source("./scripts/functions_accuracy.R")
source("./scripts/landmarks_traits.R")



#-------------------------
# PROCESS FOTO1
#-------------------------


####### MODEL1  ----------

### MODEL1 TrainingSet1 
    foto1.ts1.list <- list.files(path = "./landmarks_anova/1_foto1_ts1",
      full.names = TRUE, recursive = TRUE) 
      foto1.ts1.append = "_foto1_ts1" 
      foto1_ts1_landmarks_list <- lapply(foto1.ts1.list, process_json, 
                                         landmarks_n = 34, coord_n = 3)
    # 1-extract BILATERAL traits
    foto1_ts1_distances_bil <- c()
    foto1_ts1_distances_bil_df <- data.frame()
    for(j in 1:length(foto1_ts1_landmarks_list)){  #for each skull
      for(i in seq_along(landmark_pairs)){ 
       foto1_ts1_distances_bil[i]<- euclidean_distance(x = foto1_ts1_landmarks_list[[j]], #find each pair
                        point1 = landmark_pairs[[i]][1], point2 = landmark_pairs[[i]][2])}
     foto1_ts1_distances_bil_df  <- rbind(foto1_ts1_distances_bil_df, 
                                     foto1_ts1_distances_bil)}
    ID.foto1.ts1 <- foto1.ts1.list %>%         # Extract the ID
      basename() %>%                           # Extract only file name 
      gsub("_median\\.mrk\\.json$", "", .) %>% # Remove "_median"
      paste0(., foto1.ts1.append)              # Add "_foto1_ts1"
      colnames(foto1_ts1_distances_bil_df) <- traits
      foto1_ts1_distances_bil_df$ID <- ID.foto1.ts1
      #do the average for the bilateral traits
      foto1_ts1_distances_avg_df <- avg_columns(foto1_ts1_distances_bil_df)
      colnames(foto1_ts1_distances_avg_df)[1] <- "ID"          # Add first column back 
      colnames(foto1_ts1_distances_avg_df)[-1] <- name_traits  # Rename columns BUT no ID
    # 2- extract OTHER traits
    foto1_ts1_distances_other <- c()
    foto1_ts1_distances_other_df <- data.frame()
    for(j in 1:length(foto1_ts1_landmarks_list)){ 
       for(i in seq_along(landmark_pairs1))
         { foto1_ts1_distances_other[i] <- euclidean_distance(x = foto1_ts1_landmarks_list[[j]], 
                   point1 = landmark_pairs1[[i]][1], point2 = landmark_pairs1[[i]][2])}
        foto1_ts1_distances_other_df <- rbind(foto1_ts1_distances_other_df,
                                              foto1_ts1_distances_other)}
    ID.foto1.ts1 <- foto1.ts1.list %>%           # Extract the ID
       basename() %>%                            # Extract only the name of the file and not the path 
       gsub("_median\\.mrk\\.json$", "", .) %>%  # Removes "_median"
       paste0(., foto1.ts1.append)              
    colnames(foto1_ts1_distances_other_df) <- traits1
    foto1_ts1_distances_other_df$ID <- ID.foto1.ts1
    # 3- Unite the two
    foto1_ts1_distances_df <- merge(foto1_ts1_distances_avg_df,foto1_ts1_distances_other_df, by = c ("ID"))
    # 4- log the dataset 
    log_foto1_ts1_distances_df <- cbind(foto1_ts1_distances_df[, 1, drop = FALSE],
                                        log(foto1_ts1_distances_df[, 2:36]))
    log_foto1_ts1_distances_df$ID <- gsub("-foto1", "", log_foto1_ts1_distances_df$ID)
    # 5- plot
    log_foto1_ts1_distances_df_long <- log_foto1_ts1_distances_df %>%  
       pivot_longer(cols = -ID,  names_to = "Trait",  values_to = "Trait_Length")  
    log_foto1_ts1_distances_df_long$Trait <- gsub("_mean", "", log_foto1_ts1_distances_df_long$Trait)
    log_foto1_ts1_distances_df_long <- log_foto1_ts1_distances_df_long %>%
    mutate(Trait = factor(Trait, levels = order_traits))
              
      
### MODEL1 TrainingSet2 
    foto1.ts2.list <- list.files(path = "./landmarks_anova/2_foto1_ts2",
          pattern = "\\.json", full.names = TRUE, recursive = TRUE) 
          foto1.ts2.list
          foto1.ts2.append = "_foto1_ts2"
          foto1_ts2_landmarks_list <- lapply(foto1.ts2.list, process_json,
                                             landmarks_n = 34, coord_n = 3)
    # 1-extract BILATERAL traits
    foto1_ts2_distances_bil <- c()
    foto1_ts2_distances_bil_df <- data.frame()
    for(j in 1:length(foto1_ts2_landmarks_list))
      { for(i in seq_along(landmark_pairs)){ #find each pair
            foto1_ts2_distances_bil[i] <- euclidean_distance(x = foto1_ts2_landmarks_list[[j]], 
                               point1 = landmark_pairs[[i]][1],point2 = landmark_pairs[[i]][2])}
        foto1_ts2_distances_bil_df <- rbind(foto1_ts2_distances_bil_df, foto1_ts2_distances_bil)}
    ID.foto1.ts2 <- foto1.ts2.list %>%
       basename() %>%  
       gsub("_median\\.mrk\\.json$", "", .) %>%  
       paste0(., foto1.ts2.append)  
    colnames(foto1_ts2_distances_bil_df) <- traits
    foto1_ts2_distances_bil_df$ID <- ID.foto1.ts2
    #average
    foto1_ts2_distances_avg_df <- avg_columns(foto1_ts2_distances_bil_df)
    colnames(foto1_ts2_distances_avg_df)[1] <- "ID"         
    colnames(foto1_ts2_distances_avg_df)[-1] <- name_traits   
    # 2- extract OTHER traits
    foto1_ts2_distances_other <- c()
    foto1_ts2_distances_other_df <- data.frame()
    for(j in 1:length(foto1_ts2_landmarks_list))#for each skull
      {for(i in seq_along(landmark_pairs1)){ #find each pair
         foto1_ts2_distances_other[i] <- euclidean_distance(x = foto1_ts2_landmarks_list[[j]], 
                              point1 = landmark_pairs1[[i]][1],point2 = landmark_pairs1[[i]][2])}
       foto1_ts2_distances_other_df <- rbind(foto1_ts2_distances_other_df, foto1_ts2_distances_other)}
    ID.foto1.ts2 <- foto1.ts2.list %>%
       basename() %>%  gsub("_median\\.mrk\\.json$", "", .) %>% paste0(., foto1.ts2.append)  
    colnames(foto1_ts2_distances_other_df) <- traits1
    foto1_ts2_distances_other_df$ID <- ID.foto1.ts2
    ## 3- Unite the two
    foto1_ts2_distances_df <- merge(foto1_ts2_distances_avg_df, foto1_ts2_distances_other_df, by = c ("ID"))
    # 4- log the dataset 
    log_foto1_ts2_distances_df <- cbind(foto1_ts2_distances_df[, 1, drop = FALSE],
                                        log(foto1_ts2_distances_df[, 2:36]))
    log_foto1_ts2_distances_df$ID <- gsub("-foto1", "", log_foto1_ts2_distances_df$ID)
    ## 5- plot
    log_foto1_ts2_distances_df_long <- log_foto1_ts2_distances_df %>%  
        pivot_longer(cols = -ID,  names_to = "Trait",  values_to = "Trait_Length")  
    log_foto1_ts2_distances_df_long$Trait <- gsub("_mean", "", log_foto1_ts2_distances_df_long$Trait)
    log_foto1_ts2_distances_df_long <- log_foto1_ts2_distances_df_long %>%
    mutate(Trait = factor(Trait, levels = order_traits))
        
          
    
  
###### MODEL1 BIS   ----------

### MODEL1 BIS TrainingSet1 
    foto1bis.ts1.list <- list.files(path = "./landmarks_anova/3_foto1bis_ts1",
          pattern = "\\.json", full.names = TRUE, recursive = TRUE) 
          foto1bis.ts1.list
          foto1bis.ts1.append = "_foto1bis_ts2"
          foto1bis_ts1_landmarks_list <- lapply(foto1bis.ts1.list, process_json, landmarks_n = 34, coord_n = 3)
    # 1-extract BILATERAL traits
    foto1bis_ts1_distances_bil <- c()
    foto1bis_ts1_distances_bil_df <- data.frame()
    for(j in 1:length(foto1bis_ts1_landmarks_list))
      {for(i in seq_along(landmark_pairs)){ 
         foto1bis_ts1_distances_bil[i] <- euclidean_distance(x = foto1bis_ts1_landmarks_list[[j]], 
                               point1 = landmark_pairs[[i]][1],point2 = landmark_pairs[[i]][2])}
       foto1bis_ts1_distances_bil_df <- rbind(foto1bis_ts1_distances_bil_df,foto1bis_ts1_distances_bil)}
    ID.foto1bis.ts1 <- foto1bis.ts1.list %>%
       basename() %>%  gsub("_median\\.mrk\\.json$", "", .) %>%  paste0(., foto1bis.ts1.append)  
    colnames(foto1bis_ts1_distances_bil_df) <- traits
    foto1bis_ts1_distances_bil_df$ID <- ID.foto1bis.ts1
    #average
    foto1bis_ts1_distances_avg_df <- avg_columns(foto1bis_ts1_distances_bil_df)
    colnames(foto1bis_ts1_distances_avg_df)[1] <- "ID"         
    colnames(foto1bis_ts1_distances_avg_df)[-1] <- name_traits  
    # 2- extract OTHER traits
    foto1bis_ts1_distances_other <- c()
    foto1bis_ts1_distances_other_df <- data.frame()
    for(j in 1:length(foto1bis_ts1_landmarks_list))
      {for(i in seq_along(landmark_pairs1)){ 
          foto1bis_ts1_distances_other[i] <- euclidean_distance(x = foto1bis_ts1_landmarks_list[[j]], 
                                  point1 = landmark_pairs1[[i]][1],point2 = landmark_pairs1[[i]][2])}
       foto1bis_ts1_distances_other_df <- rbind(foto1bis_ts1_distances_other_df,foto1bis_ts1_distances_other)}
    ID.foto1bis.ts1 <- foto1bis.ts1.list %>%
      basename() %>% gsub("_median\\.mrk\\.json$", "", .) %>%  paste0(., foto1bis.ts1.append)  
    colnames(foto1bis_ts1_distances_other_df) <- traits1
    ## 3- Unite the two
    foto1bis_ts1_distances_df <- merge(foto1bis_ts1_distances_avg_df,foto1bis_ts1_distances_other_df, by = c ("ID"))
    # 4- log the dataset 
    log_foto1bis_ts1_distances_df <- cbind(foto1bis_ts1_distances_df[, 1, drop = FALSE],
                                     log(foto1bis_ts1_distances_df[, 2:36]))
    log_foto1bis_ts1_distances_df$ID <- gsub("-foto1", "", log_foto1bis_ts1_distances_df$ID)
    ## 5- plot
    log_foto1bis_ts1_distances_df_long <- log_foto1bis_ts1_distances_df %>%  
            pivot_longer(cols = -ID,  names_to = "Trait",  values_to = "Trait_Length")  
    log_foto1bis_ts1_distances_df_long$Trait <- gsub("_mean", "", log_foto1bis_ts1_distances_df_long$Trait)
    log_foto1bis_ts1_distances_df_long <- log_foto1bis_ts1_distances_df_long %>%
        mutate(Trait = factor(Trait, levels = order_traits))
      
      
### MODEL1 BIS TrainingSet2 
    foto1bis.ts2.list <- list.files(path = "./landmarks_anova/4_foto1bis_ts2",
          pattern = "\\.json", full.names = TRUE, recursive = TRUE) 
          foto1bis.ts2.list
          foto1bis.ts2.append = "_foto1bis_ts2"
          foto1bis_ts2_landmarks_list <- lapply(foto1bis.ts2.list, process_json, 
                                           landmarks_n = 34, coord_n = 3)
    # 1-extract BILATERAL traits
    foto1bis_ts2_distances_bil <- c()
    foto1bis_ts2_distances_bil_df <- data.frame()
    for(j in 1:length(foto1bis_ts2_landmarks_list)){ 
      for(i in seq_along(landmark_pairs)){ 
         foto1bis_ts2_distances_bil[i] <- euclidean_distance(x = foto1bis_ts2_landmarks_list[[j]], 
                                    point1 = landmark_pairs[[i]][1], point2 = landmark_pairs[[i]][2])}
      foto1bis_ts2_distances_bil_df <- rbind(foto1bis_ts2_distances_bil_df,foto1bis_ts2_distances_bil)}
    ID.foto1bis.ts2 <- foto1bis.ts2.list %>%
      basename() %>% gsub("_median\\.mrk\\.json$", "", .) %>%  paste0(., foto1bis.ts2.append)  
    colnames(foto1bis_ts2_distances_bil_df) <- traits
    foto1bis_ts2_distances_bil_df$ID <- ID.foto1bis.ts2
    #average
    foto1bis_ts2_distances_avg_df <- avg_columns(foto1bis_ts2_distances_bil_df)
    colnames(foto1bis_ts2_distances_avg_df)[1] <- "ID"         
    colnames(foto1bis_ts2_distances_avg_df)[-1] <- name_traits  
    ## 2- extract OTHER traits
    foto1bis_ts2_distances_other <- c()
    foto1bis_ts2_distances_other_df <- data.frame()
    for(j in 1:length(foto1bis_ts2_landmarks_list)){ 
     for(i in seq_along(landmark_pairs1)){ 
       foto1bis_ts2_distances_other[i] <- euclidean_distance(x = foto1bis_ts2_landmarks_list[[j]], 
                           point1 = landmark_pairs1[[i]][1], point2 = landmark_pairs1[[i]][2])}
     foto1bis_ts2_distances_other_df <- rbind(foto1bis_ts2_distances_other_df,foto1bis_ts2_distances_other)}
    ID.foto1bis.ts2 <- foto1bis.ts2.list %>%
       basename() %>% gsub("_median\\.mrk\\.json$", "", .) %>%  paste0(., foto1bis.ts2.append)  
    colnames(foto1bis_ts2_distances_other_df) <- traits1
    foto1bis_ts2_distances_other_df$ID <- ID.foto1bis.ts2
    ## 3- Unite the two
    foto1bis_ts2_distances_df <- merge(foto1bis_ts2_distances_avg_df,foto1bis_ts2_distances_other_df, by = c ("ID"))
    # 4- log the dataset 
    log_foto1bis_ts2_distances_df <- cbind(foto1bis_ts2_distances_df[, 1, drop = FALSE],
                                     log(foto1bis_ts2_distances_df[, 2:36]))
    log_foto1bis_ts2_distances_df$ID <- gsub("-foto1", "", log_foto1bis_ts2_distances_df$ID)
    ## 5- plot
    log_foto1bis_ts2_distances_df_long <- log_foto1bis_ts2_distances_df %>%  
       pivot_longer(cols = -ID,  names_to = "Trait",  values_to = "Trait_Length")  
    log_foto1bis_ts2_distances_df_long$Trait <- gsub("_mean", "", log_foto1bis_ts2_distances_df_long$Trait)
    log_foto1bis_ts2_distances_df_long <- log_foto1bis_ts2_distances_df_long %>%
    mutate(Trait = factor(Trait, levels = order_traits))
     
      
  
    
#-------------------------
# PROCESS FOTO2
#-------------------------
      
###### MODEL2 -----------
      
### MODEL2 TrainingSet1 
    foto2.ts1.list <- list.files(path = "./landmarks_anova/5_foto2_ts1",
        pattern = "\\.json", full.names = TRUE, recursive = TRUE) 
        foto2.ts1.append = "_foto2_ts2"
        foto2_ts1_landmarks_list <- lapply(foto2.ts1.list, process_json, landmarks_n = 34, coord_n = 3)
    ## 1- extract BILATERAL traits
    foto2_ts1_distances_bil <- c()
    foto2_ts1_distances_bil_df <- data.frame()
    for(j in 1:length(foto2_ts1_landmarks_list)){  
      for(i in seq_along(landmark_pairs)){      
        foto2_ts1_distances_bil[i] <- euclidean_distance(x = foto2_ts1_landmarks_list[[j]], 
                          point1 = landmark_pairs[[i]][1],point2 = landmark_pairs[[i]][2])}
      foto2_ts1_distances_bil_df <- rbind(foto2_ts1_distances_bil_df,foto2_ts1_distances_bil)}
    ID.foto2.ts1 <- foto2.ts1.list %>%
      basename() %>%  gsub("_median\\.mrk\\.json$", "", .) %>% paste0(., foto2.ts1.append)                 
    colnames(foto2_ts1_distances_bil_df) <- traits
    foto2_ts1_distances_bil_df$ID <- ID.foto2.ts1
    # average
    foto2_ts1_distances_avg_df <- avg_columns(foto2_ts1_distances_bil_df)
    colnames(foto2_ts1_distances_avg_df)[1] <- "ID"        
    colnames(foto2_ts1_distances_avg_df)[-1] <- name_traits 
    ## 2- extract OTHER traits
    foto2_ts1_distances_other <- c()
    foto2_ts1_distances_other_df <- data.frame()
    for(j in 1:length(foto2_ts1_landmarks_list)){  #for each skull
      for(i in seq_along(landmark_pairs1)){      #find each pair
       foto2_ts1_distances_other[i] <- euclidean_distance(x = foto2_ts1_landmarks_list[[j]], 
                            point1 = landmark_pairs1[[i]][1],point2 = landmark_pairs1[[i]][2])}
      foto2_ts1_distances_other_df <- rbind(foto2_ts1_distances_other_df, foto2_ts1_distances_other)}
    ID.foto2.ts1 <- foto2.ts1.list %>%
      basename() %>%  gsub("_median\\.mrk\\.json$", "", .) %>%paste0(., foto2.ts1.append)                 
    colnames(foto2_ts1_distances_other_df) <- traits1
    foto2_ts1_distances_other_df$ID <- ID.foto2.ts1
    ## 3 - Unite the two
    foto2_ts1_distances_df <- merge(foto2_ts1_distances_avg_df, foto2_ts1_distances_other_df, by = c ("ID"))
    # 4- log the dataset  
    log_foto2_ts1_distances_df <- cbind(foto2_ts1_distances_df[, 1, drop = FALSE],
                                        log(foto2_ts1_distances_df[, 2:36]))
    log_foto2_ts1_distances_df$ID <- gsub("-foto1", "", log_foto2_ts1_distances_df$ID)
    ## 5 - plot
    log_foto2_ts1_distances_df_long <- log_foto2_ts1_distances_df %>%  
            pivot_longer(cols = -ID,  names_to = "Trait",values_to = "Trait_Length")  
    log_foto2_ts1_distances_df_long$Trait <- gsub("_mean", "", log_foto2_ts1_distances_df_long$Trait)
    log_foto2_ts1_distances_df_long <- log_foto2_ts1_distances_df_long %>%
       mutate(Trait = factor(Trait, levels = order_traits))
      
      
### MODEL2 TrainingSet1 
    foto2.ts2.list <- list.files(path = "./landmarks_anova/6_foto2_ts2",
      pattern = "\\.json", full.names = TRUE, recursive = TRUE) 
      foto2.ts2.append = "_foto2_ts2"
      foto2_ts2_landmarks_list <- lapply(foto2.ts2.list, process_json,landmarks_n = 34, coord_n = 3)
    ## 1- extract BILATERAL traits
    foto2_ts2_distances_bil <- c()
    foto2_ts2_distances_bil_df <- data.frame()
    for(j in 1:length(foto2_ts2_landmarks_list)){  
      for(i in seq_along(landmark_pairs)){      
         foto2_ts2_distances_bil[i] <- euclidean_distance(x = foto2_ts2_landmarks_list[[j]], 
                              point1 = landmark_pairs[[i]][1],point2 = landmark_pairs[[i]][2])}
      foto2_ts2_distances_bil_df <- rbind(foto2_ts2_distances_bil_df,foto2_ts2_distances_bil)}
    ID.foto2.ts2 <- foto2.ts2.list %>%
      basename() %>%  gsub("_median\\.mrk\\.json$", "", .) %>%  paste0(., foto2.ts2.append)                 
    colnames(foto2_ts2_distances_bil_df) <- traits
    foto2_ts2_distances_bil_df$ID <- ID.foto2.ts2
    # avarage
    foto2_ts2_distances_avg_df <- avg_columns(foto2_ts2_distances_bil_df)
    colnames(foto2_ts2_distances_avg_df)[1] <- "ID"        
    colnames(foto2_ts2_distances_avg_df)[-1] <- name_traits 
    ## 2- extract OTHER traits
    foto2_ts2_distances_other <- c()
    foto2_ts2_distances_other_df <- data.frame()
    for(j in 1:length(foto2_ts2_landmarks_list)){  #for each skull
      for(i in seq_along(landmark_pairs1)){      #find each pair
        foto2_ts2_distances_other[i] <- euclidean_distance(x = foto2_ts2_landmarks_list[[j]], 
                              point1 = landmark_pairs1[[i]][1], point2 = landmark_pairs1[[i]][2])}
      foto2_ts2_distances_other_df <- rbind(foto2_ts2_distances_other_df,foto2_ts2_distances_other)}
    ID.foto2.ts2 <- foto2.ts2.list %>%
      basename() %>% gsub("_median\\.mrk\\.json$", "", .) %>%  paste0(., foto2.ts2.append)                 
    colnames(foto2_ts2_distances_other_df) <- traits1
    foto2_ts2_distances_other_df$ID <- ID.foto2.ts2
    ## 3 - Unite the two
    foto2_ts2_distances_df <- merge(foto2_ts2_distances_avg_df,foto2_ts2_distances_other_df, by = c ("ID"))
    ## 4- log the dataset  
    log_foto2_ts2_distances_df <- cbind(foto2_ts2_distances_df[, 1, drop = FALSE],
                                              log(foto2_ts2_distances_df[, 2:36]))
    log_foto2_ts2_distances_df$ID <- gsub("-foto1", "", log_foto2_ts2_distances_df$ID)
    ## 5- plot
    log_foto2_ts2_distances_df_long <- log_foto2_ts2_distances_df %>%  
            pivot_longer(cols = -ID,  names_to = "Trait", values_to = "Trait_Length")  
    log_foto2_ts2_distances_df_long$Trait <- gsub("_mean", "",log_foto2_ts2_distances_df_long$Trait)
    log_foto2_ts2_distances_df_long <- log_foto2_ts2_distances_df_long %>%
        mutate(Trait = factor(Trait, levels = order_traits))
    
         
      
####### MODEL2 BIS -----------
      
### MODEL2 BIS TrainingSet1 
    foto2bis.ts1.list <- list.files( path = "./landmarks_anova/7_foto2bis_ts1",
          pattern = "\\.json", full.names = TRUE, recursive = TRUE) 
          foto2bis.ts1.list
          foto2bis.ts1.append = "-foto2bis_ts1"
          foto2bis_ts1_landmarks_list <- lapply(foto2bis.ts1.list, process_json,
                                            landmarks_n = 34, coord_n = 3)
    ## 1- BILATERAL traits
    foto2bis_ts1_distances_bil <- c()
    foto2bis_ts1_distances_bil_df <- data.frame()
    for(j in 1:length(foto2bis_ts1_landmarks_list)){ #for each skull
      for(i in seq_along(landmark_pairs)){ #find each pair
       foto2bis_ts1_distances_bil[i] <- euclidean_distance(x = foto2bis_ts1_landmarks_list[[j]], 
                      point1 = landmark_pairs[[i]][1], point2 = landmark_pairs[[i]][2])}
      foto2bis_ts1_distances_bil_df <- rbind(foto2bis_ts1_distances_bil_df,foto2bis_ts1_distances_bil)}
    ID.foto2bis.ts1 <- foto2bis.ts1.list %>%
      basename() %>%  gsub("_median\\.mrk\\.json$", "", .) %>% paste0(., foto2bis.ts1.append)  
    colnames(foto2bis_ts1_distances_bil_df) <- traits
    foto2bis_ts1_distances_bil_df$ID <- ID.foto2bis.ts1
    #average
    foto2bis_ts1_distances_avg_df <- avg_columns(foto2bis_ts1_distances_bil_df)
    colnames(foto2bis_ts1_distances_avg_df)[1] <- "ID"         
    colnames(foto2bis_ts1_distances_avg_df)[-1] <- name_traits  
    ## 2- extract OTHER traits
    foto2bis_ts1_distances_other <- c()
    foto2bis_ts1_distances_other_df <- data.frame()
    for(j in 1:length(foto2bis_ts1_landmarks_list)){ #for each skull
     for(i in seq_along(landmark_pairs1)){ #find each pair
       foto2bis_ts1_distances_other[i] <- euclidean_distance(x = foto2bis_ts1_landmarks_list[[j]], 
                               point1 = landmark_pairs1[[i]][1],point2 = landmark_pairs1[[i]][2])}
     foto2bis_ts1_distances_other_df <- rbind(foto2bis_ts1_distances_other_df, foto2bis_ts1_distances_other)}
    ID.foto2bis.ts1 <- foto2bis.ts1.list %>%
       basename() %>% gsub("_median\\.mrk\\.json$", "", .) %>%  paste0(., foto2bis.ts1.append)  
    colnames(foto2bis_ts1_distances_other_df) <- traits1
    foto2bis_ts1_distances_other_df$ID <- ID.foto2bis.ts1
    ## 3- Unite the two
    foto2bis_ts1_distances_df <- merge(foto2bis_ts1_distances_avg_df,foto2bis_ts1_distances_other_df, by = c ("ID"))
    ## 4- log the dataset 
    log_foto2bis_ts1_distances_df <- cbind(foto2bis_ts1_distances_df[, 1, drop = FALSE],
                                           log(foto2bis_ts1_distances_df[, 2:36]))
    log_foto2bis_ts1_distances_df$ID <- gsub("-foto1", "", log_foto2bis_ts1_distances_df$ID)
    ## 5- plot
    log_foto2bis_ts1_distances_df_long <- log_foto2bis_ts1_distances_df %>%  
          pivot_longer(cols = -ID,  names_to = "Trait",  values_to = "Trait_Length")  
    log_foto2bis_ts1_distances_df_long$Trait <- gsub("_mean", "", log_foto2bis_ts1_distances_df_long$Trait)
    log_foto2bis_ts1_distances_df_long <- log_foto2bis_ts1_distances_df_long %>%
         mutate(Trait = factor(Trait, levels = order_traits))
          
      
### MODEL2 BIS TrainingSet2       
    foto2bis.ts2.list <- list.files(path = "./landmarks_anova/8_foto2bis_ts2",
          pattern = "\\.json", full.names = TRUE, recursive = TRUE) 
          foto2bis.ts2.list
          foto2bis.ts2.append = "-foto2bis_ts2"
          foto2bis_ts2_landmarks_list <- lapply(foto2bis.ts2.list, process_json,
                                            landmarks_n = 34, coord_n = 3)
    ## 1- BILATERAL traits
    foto2bis_ts2_distances_bil <- c()
    foto2bis_ts2_distances_bil_df <- data.frame()
    for(j in 1:length(foto2bis_ts2_landmarks_list)){ #for each skull
     for(i in seq_along(landmark_pairs)){ #find each pair
       foto2bis_ts2_distances_bil[i] <- euclidean_distance(x = foto2bis_ts2_landmarks_list[[j]], 
                              point1 = landmark_pairs[[i]][1],point2 = landmark_pairs[[i]][2])}
     foto2bis_ts2_distances_bil_df <- rbind(foto2bis_ts2_distances_bil_df,foto2bis_ts2_distances_bil)}
    ID.foto2bis.ts2 <- foto2bis.ts2.list %>%
       basename() %>% gsub("_median\\.mrk\\.json$", "", .) %>% paste0(., foto2bis.ts2.append)  
    colnames(foto2bis_ts2_distances_bil_df) <- traits
    foto2bis_ts2_distances_bil_df$ID <- ID.foto2bis.ts2
    #average
    foto2bis_ts2_distances_avg_df <- avg_columns(foto2bis_ts2_distances_bil_df)
    colnames(foto2bis_ts2_distances_avg_df)[1] <- "ID"         
    colnames(foto2bis_ts2_distances_avg_df)[-1] <- name_traits  
    ## 2- extract OTHER traits
    foto2bis_ts2_distances_other <- c()
    foto2bis_ts2_distances_other_df <- data.frame()
    for(j in 1:length(foto2bis_ts2_landmarks_list)){ #for each skull
     for(i in seq_along(landmark_pairs1)){ #find each pair
       foto2bis_ts2_distances_other[i] <- euclidean_distance(x = foto2bis_ts2_landmarks_list[[j]], 
                                  point1 = landmark_pairs1[[i]][1], point2 = landmark_pairs1[[i]][2])}
     foto2bis_ts2_distances_other_df <- rbind(foto2bis_ts2_distances_other_df, foto2bis_ts2_distances_other)}
    ID.foto2bis.ts2 <- foto2bis.ts2.list %>%
       basename() %>% gsub("_median\\.mrk\\.json$", "", .) %>%paste0(., foto2bis.ts2.append)  
    colnames(foto2bis_ts2_distances_other_df) <- traits1
    foto2bis_ts2_distances_other_df$ID <- ID.foto2bis.ts2
    ## 3- Unite the two
    foto2bis_ts2_distances_df <- merge(foto2bis_ts2_distances_avg_df,foto2bis_ts2_distances_other_df, by = c ("ID"))
    ## 4- log the dataset 
    log_foto2bis_ts2_distances_df <- cbind(foto2bis_ts2_distances_df[, 1, drop = FALSE],
                                           log(foto2bis_ts2_distances_df[, 2:36]))
    log_foto2bis_ts2_distances_df$ID <- gsub("-foto1", "", log_foto2bis_ts2_distances_df$ID)
    ## 5- plot
    log_foto2bis_ts2_distances_df_long <- log_foto2bis_ts2_distances_df %>%  
            pivot_longer(cols = -ID,  names_to = "Trait",  values_to = "Trait_Length")  
    log_foto2bis_ts2_distances_df_long$Trait <- gsub("_mean", "",log_foto2bis_ts2_distances_df_long$Trait)
    log_foto2bis_ts2_distances_df_long <- log_foto2bis_ts2_distances_df_long %>%
        mutate(Trait = factor(Trait, levels = order_traits))
          
          
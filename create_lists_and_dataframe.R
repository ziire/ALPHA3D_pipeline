##########################################################
##   Automatizing phenomic research: ALPHA3D pipeline   ##
##                                                      ##
##          CREATE LISTS AND DATAFRAMES                 ## 
##                                                      ##
##                  Irene Zanandrea                     ##
##########################################################

## load packages, functions and landmark names
source("./scripts/libraries.R")
source("./scripts/functions_accuracy.R")
source("./scripts/landmarks_traits.R")


#-----------------------------------------------------
# DATAFRAMES OF FOTO 1 AND FOTO 2 (TOTAL IMPRECISION)
#-----------------------------------------------------

###### A- FOTO 1   ---------------

foto1.list <- list.files(path = "./landmark/1_Etot_v2/FOTO1", full.names = TRUE, recursive = TRUE) 
foto1.append = "-foto1" 
foto1_landmarks_list <- lapply(foto1.list, 
                               process_json, 
                               landmarks_n = 34, 
                               coord_n = 3)
    #list of n skulls that contains n skulls that are 3D array 
    #containing 34 landmarks with 3 coordinates  -> 3D array with dimensions[34,3,1] (1json=1 specimen) 
    #34 rows = landmarks, 3 columns = coordinates (X, Y, Z), 1 slice=1 specimen 
length(foto1_landmarks_list)     # number of JSON files processed(number of skulls)
dim(foto1_landmarks_list[[1]])   # should be (34, 3, 1)
## 1- extract BILATERAL traits
foto1_distances_bil <- c()
foto1_distances_bil_df <- data.frame()
for(j in 1:length(foto1_landmarks_list)) #for each skull
  { for(i in seq_along(landmark_pairs)){ #find each pair
    foto1_distances_bil[i]<- euclidean_distance(x = foto1_landmarks_list[[j]], 
                            point1 = landmark_pairs[[i]][1], point2 = landmark_pairs[[i]][2])
    }
  foto1_distances_bil_df  <- rbind(foto1_distances_bil_df,foto1_distances_bil)}
ID.foto1 <- foto1.list %>%                 # Extract the ID
  basename() %>%                           # Extract only file name 
  gsub("_median\\.mrk\\.json$", "", .) %>% # Remove "_median"
  paste0(., foto1.append)                  # Add "-foto1"
colnames(foto1_distances_bil_df) <- traits
foto1_distances_bil_df$ID <- ID.foto1
length(foto1_distances_bil_df)             # number of columns; last column is ID 
# compute the average
foto1_distances_avg_df <- avg_columns(foto1_distances_bil_df)
colnames(foto1_distances_avg_df)[1] <- "ID"          # Add first column back 
colnames(foto1_distances_avg_df)[-1] <- name_traits  # Rename columns BUT no ID
length(foto1_distances_avg_df)
## 2- extract OTHER traits
foto1_distances_other <- c()
foto1_distances_other_df <- data.frame()
for(j in 1:length(foto1_landmarks_list))
  { for(i in seq_along(landmark_pairs1)){ #find each pair
    foto1_distances_other[i] <- euclidean_distance(x = foto1_landmarks_list[[j]], 
                          point1 = landmark_pairs1[[i]][1], point2 = landmark_pairs1[[i]][2])}
  foto1_distances_other_df <- rbind(foto1_distances_other_df, foto1_distances_other)}
ID.foto1 <- foto1.list %>%                  # Extract the ID
  basename() %>%                            # Estrae solo il nome del file e non il percorso
  gsub("_median\\.mrk\\.json$", "", .) %>%  # Rimuove "_median"
  paste0(., foto1.append)                   # Aggiunge "-foto1" al numero
colnames(foto1_distances_other_df) <- traits1
foto1_distances_other_df$ID <- ID.foto1
length(foto1_distances_other_df)
### 3 - Unite the two
foto1_distances_df <- merge(foto1_distances_avg_df,foto1_distances_other_df,by = c ("ID"))
length(foto1_distances_df)
### 4- log the dataset
log_foto1_distances_df <- cbind(foto1_distances_df[, 1, drop = FALSE], log(foto1_distances_df[, 2:36]))
log_foto1_distances_df$ID <- gsub("-foto1", "", log_foto1_distances_df$ID)



###### B - FOTO 2   ---------------
foto2.list <- list.files(path = "./landmark/1_Etot_v2/FOTO2",
                         pattern = "\\.json", full.names = TRUE, recursive = TRUE) 
foto2.append = "-foto2"
foto2_landmarks_list <- lapply(foto2.list, process_json, landmarks_n = 34, coord_n = 3)
## 1- extract BILATERAL traits
foto2_distances_bil <- c()
foto2_distances_bil_df <- data.frame()
for(j in 1:length(foto2_landmarks_list)){  
  for(i in seq_along(landmark_pairs)){      
    foto2_distances_bil[i] <- euclidean_distance(x = foto2_landmarks_list[[j]], 
                      point1 = landmark_pairs[[i]][1],  point2 = landmark_pairs[[i]][2])}
  foto2_distances_bil_df <- rbind(foto2_distances_bil_df, 
                                  foto2_distances_bil)}
ID.foto2 <- foto2.list %>%
  basename() %>%  gsub("_median\\.mrk\\.json$", "", .) %>%  paste0(., foto2.append)                 
colnames(foto2_distances_bil_df) <- traits
foto2_distances_bil_df$ID <- ID.foto2
# compute the average
foto2_distances_avg_df <- avg_columns(foto2_distances_bil_df)
colnames(foto2_distances_avg_df)[1] <- "ID"        
colnames(foto2_distances_avg_df)[-1] <- name_traits 
## 2- extract OTHER traits
foto2_distances_other <- c()
foto2_distances_other_df <- data.frame()
for(j in 1:length(foto2_landmarks_list)){  #for each skull
  for(i in seq_along(landmark_pairs1)){      #find each pair
    foto2_distances_other[i] <- euclidean_distance(x = foto2_landmarks_list[[j]], 
                        point1 = landmark_pairs1[[i]][1], point2 = landmark_pairs1[[i]][2])}
  foto2_distances_other_df <- rbind(foto2_distances_other_df,foto2_distances_other)}
ID.foto2 <- foto2.list %>%
  basename() %>%  gsub("_median\\.mrk\\.json$", "", .) %>% paste0(., foto2.append)                 
colnames(foto2_distances_other_df) <- traits1
foto2_distances_other_df$ID <- ID.foto2
## 3 - Unite the two
foto2_distances_df <- merge(foto2_distances_avg_df,foto2_distances_other_df, by = c ("ID"))
### 4- log the dataset
log_foto2_distances_df <- cbind(foto2_distances_df[, 1, drop = FALSE],
                                log(foto2_distances_df[, 2:36]))
log_foto2_distances_df$ID <- gsub("-foto2", "", log_foto2_distances_df$ID)




#-------------------------------------------------------------
# DATAFRAME OF FOTO 1BIS (IMPRECISION CAUSED BY 3D AND AL)
#-------------------------------------------------------------

fotobis.list <- list.files(path = "./landmark/2_E3D_AL_v2",
                           pattern = "\\.json", full.names = TRUE, recursive = TRUE) 
fotobis.append = "-fotobis"
fotobis_landmarks_list <- lapply(fotobis.list, process_json,landmarks_n = 34, coord_n = 3)
####### 1- BILATERAL traits
fotobis_distances_bil <- c()
fotobis_distances_bil_df <- data.frame()
for(j in 1:length(fotobis_landmarks_list)){ 
  for(i in seq_along(landmark_pairs)){ 
    fotobis_distances_bil[i] <- euclidean_distance(x = fotobis_landmarks_list[[j]], 
                     point1 = landmark_pairs[[i]][1],point2 = landmark_pairs[[i]][2])}
  fotobis_distances_bil_df <- rbind(fotobis_distances_bil_df, fotobis_distances_bil)}
ID.fotobis <- fotobis.list %>%
  basename() %>% gsub("_median\\.mrk\\.json$", "", .) %>%  paste0(., fotobis.append)  
colnames(fotobis_distances_bil_df) <- traits
fotobis_distances_bil_df$ID <- ID.fotobis
#average
fotobis_distances_avg_df <- avg_columns(fotobis_distances_bil_df)
colnames(fotobis_distances_avg_df)[1] <- "ID"         
colnames(fotobis_distances_avg_df)[-1] <- name_traits  
## 2- extract OTHER traits
fotobis_distances_other <- c()
fotobis_distances_other_df <- data.frame()
for(j in 1:length(fotobis_landmarks_list)){ 
  for(i in seq_along(landmark_pairs1)){ 
    fotobis_distances_other[i] <- euclidean_distance(x = fotobis_landmarks_list[[j]], 
                       point1 = landmark_pairs1[[i]][1],point2 = landmark_pairs1[[i]][2])}
  fotobis_distances_other_df <- rbind(fotobis_distances_other_df, fotobis_distances_other)}
ID.fotobis <- fotobis.list %>%
  basename() %>% gsub("_median\\.mrk\\.json$", "", .) %>%paste0(., fotobis.append)  
colnames(fotobis_distances_other_df) <- traits1
fotobis_distances_other_df$ID <- ID.fotobis
## 3- Unite the two
fotobis_distances_df <- merge(fotobis_distances_avg_df,fotobis_distances_other_df, by = c ("ID"))
## 4- log the dataset
log_fotobis_distances_df <- cbind(fotobis_distances_df[, 1, drop = FALSE],
                                  log(fotobis_distances_df[, 2:36]))
log_fotobis_distances_df$ID <- gsub("-fotobis", "", log_fotobis_distances_df$ID)



#--------------------------------------------------------------
# DATAFRAME OF FOTO 1EXTRA (IMPRECISION CAUSED ONLY BY AL)
#--------------------------------------------------------------

######## A- FOTOEXTRA
fotoextra.list <- list.files(path = "./landmark/3_EAL_v2",
                             pattern = "\\.json", full.names = TRUE, recursive = TRUE) 
fotoextra.list
fotoextra.append = "-fotoextra"
fotoextra_landmarks_list <- lapply(fotoextra.list, process_json,
                                    landmarks_n = 34, coord_n = 3)
## 1- BILATERAL traits
fotoextra_distances_bil <- c()
fotoextra_distances_bil_df <- data.frame()
for(j in 1:length(fotoextra_landmarks_list)){ #for each skull
  for(i in seq_along(landmark_pairs)){ #find each pair
    fotoextra_distances_bil[i] <- euclidean_distance(x = fotoextra_landmarks_list[[j]], 
                      point1 = landmark_pairs[[i]][1], point2 = landmark_pairs[[i]][2])}
  fotoextra_distances_bil_df <- rbind(fotoextra_distances_bil_df, fotoextra_distances_bil)}
ID.fotoextra <- fotoextra.list %>%
  basename() %>%gsub("_median\\.mrk\\.json$", "", .) %>%paste0(., fotoextra.append)  
colnames(fotoextra_distances_bil_df) <- traits
fotoextra_distances_bil_df$ID <- ID.fotoextra
#average
fotoextra_distances_avg_df <- avg_columns(fotoextra_distances_bil_df)
colnames(fotoextra_distances_avg_df)[1] <- "ID"         
colnames(fotoextra_distances_avg_df)[-1] <- name_traits  
## 2- extract OTHER traits
fotoextra_distances_other <- c()
fotoextra_distances_other_df <- data.frame()
for(j in 1:length(fotoextra_landmarks_list)){ #for each skull
  for(i in seq_along(landmark_pairs1)){ #find each pair
    fotoextra_distances_other[i] <- euclidean_distance(x = fotoextra_landmarks_list[[j]], 
                    point1 = landmark_pairs1[[i]][1],point2 = landmark_pairs1[[i]][2])}
  fotoextra_distances_other_df <- rbind(fotoextra_distances_other_df,fotoextra_distances_other)}
ID.fotoextra <- fotoextra.list %>%
  basename() %>% gsub("_median\\.mrk\\.json$", "", .) %>% paste0(., fotoextra.append)  
colnames(fotoextra_distances_other_df) <- traits1
fotoextra_distances_other_df$ID <- ID.fotoextra
## 3- Unite the two
fotoextra_distances_df <- merge(fotoextra_distances_avg_df,fotoextra_distances_other_df, by = c ("ID"))
## 4- log the dataset
log_fotoextra_distances_df <- cbind(fotoextra_distances_df[, 1, drop = FALSE],
                                    log(fotoextra_distances_df[, 2:36]))
log_fotoextra_distances_df$ID <- gsub("-fotoextra", "", log_fotoextra_distances_df$ID)




##########################################################
##   Automatizing phenomic research: ALPHA3D pipeline   ##
##                                                      ##
##           CS CALCULATIONS and ANOVAS                 ## 
##                                                      ##
##                  Irene Zanandrea                     ##
##########################################################


## LOAD packages, functions, landmarks, lists and dataframe
source("./scripts/libraries.R")
source("./scripts/functions_accuracy.R")
source("./scripts/landmarks_traits.R")
source("./scripts/list_nested_anova.R")
source("./scripts/imprecision_2nd_part_pipeline")




#Nested ANOVA on CS: 
#analysis of size  -> compute CS as the response.

#Procrustes ANOVA: 
#analysis of shape coordinates -> the response is the Procrustes coordinates.



#-------------------------------------
# CENTROID SIZE CALCULATION
#--------------------------------------

# to compute the centroid size in geomorph -> convert the list in an array

## 1. Extract the ID from my dataframe
ID_foto <- log_foto1_ts1_distances_df$ID
ID_foto <- gsub("_foto1_ts1", "", ID_foto)
ID_foto
## 2. Compute CS

###### FOTO 1   ---------------

## MODEL1 
   ## TS1 (FOTO1)
   n_spec_foto1_ts1 <- length(foto1_ts1_landmarks_list)
   coord_array_foto1_ts1 <- array(data = unlist(foto1_ts1_landmarks_list),
                                 dim  = c(34, 3, n_spec_foto1_ts1),
                                 dimnames = list(landmark = NULL,
                                                 axis     = c("x", "y", "z"),
                                                 specimen = ID_foto))
   dim(coord_array_foto1_ts1) # should be 34x3xn_spec_foto1_TS1 (p×k×n_landmark array)
   #calculate CS for GPA (General Procrustes Analysis)
   Y.gpa.foto1.ts1 <- gpagen(coord_array_foto1_ts1)
   CS_foto1_ts1    <- Y.gpa.foto1.ts1$Csize
   CS_foto1_ts1
   CS_foto1_ts1_log <- log(CS_foto1_ts1)
   
   ## TS2 (FOTOEXTRA)
   n_spec_foto1_ts2 <- length(foto1_ts2_landmarks_list)
   coord_array_foto1_ts2 <- array(data = unlist(foto1_ts2_landmarks_list),
                dim  = c(34, 3, n_spec_foto1_ts2),dimnames = list(
                landmark = NULL, axis = c("x", "y", "z"),specimen = ID_foto))
   Y.gpa.foto1.ts2 <- gpagen(coord_array_foto1_ts2)
   CS_foto1_ts2    <- Y.gpa.foto1.ts2$Csize
   CS_foto1_ts2
   CS_foto1_ts2_log <- log(CS_foto1_ts2)
      
   
## MODEL1BIS
   ## TS1 
   n_spec_foto1bis_ts1 <- length(foto1bis_ts1_landmarks_list)
   coord_array_foto1bis_ts1 <- array(data = unlist(foto1bis_ts1_landmarks_list),
                dim  = c(34, 3, n_spec_foto1bis_ts1),dimnames = list(
                landmark = NULL, axis= c("x", "y", "z"),specimen = ID_foto))
   Y.gpa.foto1bis_ts1 <- gpagen(coord_array_foto1bis_ts1)
   CS_foto1bis_ts1    <- Y.gpa.foto1bis_ts1$Csize    
   CS_foto1bis_ts1
   CS_foto1bis_ts1_log <- log(CS_foto1bis_ts1)
   
   ## TS2 (FOTOBIS)
   n_spec_foto1bis_ts2<- length(foto1bis_ts2_landmarks_list)
   coord_array_foto1bis_ts2 <- array(data = unlist(foto1bis_ts2_landmarks_list),
                dim  = c(34, 3, n_spec_foto1bis_ts2),dimnames = list(
                landmark = NULL, axis= c("x", "y", "z"),specimen = ID_foto))
   Y.gpa.foto1bis.ts2 <- gpagen(coord_array_foto1bis_ts2)
   CS_foto1bis_ts2    <- Y.gpa.foto1bis.ts2$Csize    
   CS_foto1bis_ts2
   CS_foto1bis_ts2_log <- log(CS_foto1bis_ts2)
 
      

###### FOTO 2   ---------------
      
## MODEL2
   ## TS1 
   n_spec_foto2_ts1 <- length(foto2_ts1_landmarks_list)
   coord_array_foto2_ts1 <- array(data = unlist(foto2_ts1_landmarks_list),
                 dim  = c(34, 3, n_spec_foto2_ts1),dimnames = list(
                landmark = NULL, axis= c("x", "y", "z"), specimen = ID_foto))
   dim(coord_array_foto2_ts1) 
   Y.gpa.foto2.ts1 <- gpagen(coord_array_foto2_ts1)
   CS_foto2_ts1    <- Y.gpa.foto2.ts1$Csize
   CS_foto2_ts1
   CS_foto2_ts1_log <- log(CS_foto2_ts1)
      
   ## TS2 (FOTO2)
   n_spec_foto2_ts2 <- length(foto2_ts2_landmarks_list)
   coord_array_foto2_ts2 <- array(data = unlist(foto2_ts2_landmarks_list),
                 dim  = c(34, 3, n_spec_foto2_ts2),dimnames = list(
                 landmark = NULL, axis= c("x", "y", "z"), specimen = ID_foto))
   dim(coord_array_foto2_ts2) 
   Y.gpa.foto2.ts2 <- gpagen(coord_array_foto2_ts2)
   CS_foto2_ts2    <- Y.gpa.foto2.ts2$Csize
   CS_foto2_ts2
   CS_foto2_ts2_log <- log(CS_foto2_ts2)
      
   
## MODEL2 BIS
   ## TS1  
   n_spec_foto2bis_ts1 <- length(foto2bis_ts1_landmarks_list)
   coord_array_foto2bis_ts1 <- array(data = unlist(foto2bis_ts1_landmarks_list),
                 dim  = c(34, 3, n_spec_foto2bis_ts1), dimnames = list(
                 landmark = NULL, axis= c("x", "y", "z"),specimen = ID_foto))
   Y.gpa.foto2bis.ts1 <- gpagen(coord_array_foto2bis_ts1)
   CS_foto2bis_ts1    <- Y.gpa.foto2bis.ts1$Csize
   CS_foto2bis_ts1
   CS_foto2bis_ts1_log <- log(CS_foto2bis_ts1)
   
   ## TS2 
   n_spec_foto2bis_ts2 <- length(foto2bis_ts2_landmarks_list)
   coord_array_foto2bis_ts2 <- array(data = unlist(foto2bis_ts2_landmarks_list),
                 dim  = c(34, 3, n_spec_foto2bis_ts2), dimnames = list(
                 landmark = NULL, axis= c("x", "y", "z"),specimen = ID_foto))
   Y.gpa.foto2bis.ts2 <- gpagen(coord_array_foto2bis_ts2)
   CS_foto2bis_ts2   <- Y.gpa.foto2bis.ts2$Csize
   CS_foto2bis_ts2
   CS_foto2bis_ts2_log <- log(CS_foto2bis_ts2)
  
  
  
   
#-----------------------------------------
# COMBINE THE CS IN ONE DATAFRAME
#------------------------------------------    

# Individual = biological variation
# Photo set = foto1 vs foto2
# Replicate = foto1 vs foto1bis (and foto2 vs foto2bis)
# Training set = ts1 vs ts2

##  1. Combine all the CS into one data frame -> long format
n <- length(CS_foto1_ts1)   # Make sure n is correct -> should be 33

#Combine CS into a matrix: rows = skulls, columns = 8 conditions
CS_matrix<- cbind(CS_foto1_ts1,
                  CS_foto1_ts2,
                  CS_foto1bis_ts1,
                  CS_foto1bis_ts2,
                  CS_foto2_ts1,
                  CS_foto2_ts2,
                  CS_foto2bis_ts1,
                  CS_foto2bis_ts2)   
dim(CS_matrix) #should be 33 x 8

# make it as a vector -> for each skull I put its 8 replicates in sequence
Csize <- as.vector(t(CS_matrix))   # length = 8xn = 264

# Factors in the same order
# id = individual
  id <- factor(rep(ID_foto, each = 8))   # 33 skulls x 8 = 264 rows
# foto = the photosets
  foto <- factor(rep(c("foto1", "foto1", "foto1", "foto1",
                       "foto2", "foto2", "foto2", "foto2"),times = n))
# replicate = original or fotobis
  replicate <- factor(rep(c("orig", "orig", "bis", "bis",
                            "orig", "orig", "bis", "bis"),times = n))
# ts = training set
  ts <- factor(rep(c("ts1", "ts2", "ts1", "ts2", "ts1", "ts2", "ts1", "ts2"),
                   times = n))
#build the dataframe
df <- data.frame(Csize     = Csize,
                 id        = id,
                 foto      = foto,
                 replicate = replicate,
                 ts        = ts)

table(df$id)  # Check counts per id -> they should all be 8
with(df, table(id, foto, replicate, ts))# Check that each id has exactly one observation per combination

#Check some statistics :)
tapply(df$Csize, df$id, mean)       # mean Csize per skull
tapply(df$Csize, df$id, sd)         # sd per skull
summary(tapply(df$Csize, df$id, mean))

### mean, min, max and standard error per cranium
# Function for min, max and standard error
f_stats <- function(x) {
  c(mean = mean(x, na.rm = TRUE), #calculate the MEAN value of CS
    min  = min(x, na.rm = TRUE),  #calculate the MINIMUM value of CS
    max  = max(x, na.rm = TRUE),  #calculate the MAXIMUM value of CS
    se   = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))} #calculate the SE

tab_CS <- tapply(df$Csize, df$id, f_stats)

# to have a nice dataframe:
tab_CS_df <- do.call(rbind, tab_CS)
tab_CS_df
#write.csv(tab_CS_df, file = "Csize_stats_per_skull.csv", row.names = TRUE)





#---------------------------
#    NESTED ANOVA
#---------------------------

###In the nested ANOVA I use CS

#make sure they are factor
  df$id        <- factor(df$id)
  df$foto      <- factor(df$foto)      # "foto1", "foto2"
  df$replicate <- factor(df$replicate) # "orig", "bis"
  df$ts        <- factor(df$ts)        # "ts1", "ts2"

#nested ANOVA
nested_anova <- aov(Csize ~ 1 + Error(id/foto/replicate/ts),     # random nesting structure                 
                    data = df)
aov_tbl <- summary(nested_anova) 
  
  
##### get F-value and p-values
  
  ## FUNCTION to   extract the residual row (df and MS) from a table

  get_residual <- function(tbl) {     # each tbl is a data.frame; the row we need is called "Residuals"
    tbl["Residuals", c("Df", "Sum Sq", "Mean Sq")]  }

  # Individual 
  df_id_tbl <- aov_tbl[[1]][[1]]           # the data.frame inside the list
  res_id <- get_residual(df_id_tbl)
  SS_id <- as.numeric(res_id["Sum Sq"])
  df_id <- as.numeric(res_id["Df"])
  MS_id<- as.numeric(res_id["Mean Sq"])
  # Photoset
  df_foto_tbl <- aov_tbl[[2]][[1]]
  res_foto <- get_residual(df_foto_tbl)
  SS_foto <- as.numeric(res_foto["Sum Sq"])
  df_foto <- as.numeric(res_foto["Df"])
  MS_foto <- as.numeric(res_foto["Mean Sq"])
  # Replicate 
  df_rep_tbl<- aov_tbl[[3]][[1]]
  res_rep <- get_residual(df_rep_tbl)
  SS_rep <- as.numeric(res_rep["Sum Sq"])
  df_rep <- as.numeric(res_rep["Df"])
  MS_rep <- as.numeric(res_rep["Mean Sq"])
  # Training set 
  df_ts_tbl <- aov_tbl[[4]][[1]]
  res_ts<- get_residual(df_ts_tbl)
  SS_ts <- as.numeric(res_ts["Sum Sq"])
  df_ts <- as.numeric(res_ts["Df"])
  MS_ts <- as.numeric(res_ts["Mean Sq"])
  
  ##### %SS
  total_SS <- SS_id + SS_foto + SS_rep + SS_ts 
  pct_id <- 100 * SS_id   / total_SS
  pct_foto <- 100 * SS_foto / total_SS
  pct_rep <- 100 * SS_rep  / total_SS
  pct_ts <- 100 * SS_ts   / total_SS
 
  ##### F and p-values
  # Individual vs. Photoset
  F_id <- MS_id   / MS_foto
  p_id <- pf(F_id, df1 = df_id,   df2 = df_foto, lower.tail = FALSE)
  # Photoset vs. Replicate
  F_foto <- MS_foto / MS_rep
  p_foto <- pf(F_foto, df1 = df_foto, df2 = df_rep, lower.tail = FALSE)
  # Replicate vs. Training set
  F_rep <- MS_rep  / MS_ts
  p_rep <- pf(F_rep, df1 = df_rep, df2 = df_ts, lower.tail = FALSE)
  
##final table
final_ANOVA <- data.frame(
    Effect  = c("Individual", "Photoset", "Replicate","Training set"),
    SS      = c(SS_id,   SS_foto,   SS_rep,   SS_ts),
    `%SS`   = c(pct_id,  pct_foto,  pct_rep,  pct_ts),
    MS      = c(MS_id,   MS_foto,   MS_rep,   MS_ts),
    df      = c(df_id,   df_foto,   df_rep,   df_ts),
    F       = c(F_id,    F_foto,    F_rep,    NA),
    `Pr(>F)`= c(p_id,    p_foto,    p_rep,    NA) )
final_ANOVA
  
  




#--------------------------------
#     PROCRUSTES ANOVA
#--------------------------------

# In the procrustes ANOVA I study the variation in landmark positions 
    
#Combine all the coordinates into one big 3D array (34×3×256)
coords_all <- abind(coord_array_foto1_ts1,
                    coord_array_foto1_ts2,
                    coord_array_foto1bis_ts1,
                    coord_array_foto1bis_ts2,
                    coord_array_foto2_ts1,
                    coord_array_foto2_ts2,
                    coord_array_foto2bis_ts1,
                    coord_array_foto2bis_ts2,
                    along = 3)
dim(coords_all) # should be: 34 3 272
 
#### Generalized Procrustes Analysis (GPA)
#library(geomorph)
#gpagen: performs the generalized Procrustes superimposition,
         #centers, scales, rotates,
         #outputs Procrustes shape coordinates and centroid size.
gpa_all <- gpagen(coords_all)  #gpagen takes landmarks and performs a GPA
str(gpa_all)  # $coords:Procrustes coordinates (shape).
              # $Csize:centroid size for each specimen.
              # $p, $k, $n: number of landmarks, dimensions, specimens.

plot(gpa_all)             # SUPER COOL 3D visual check of the superimposition

Csize   <- gpa_all$Csize
#If I want to check distribution of centroid size
hist(Csize,
     breaks = 16, # più barre distanziate
     col    = "lightblue",
     border = "white",
     main   = "Distribtion of Centroid size",
     xlab   = "Centroid size",
     ylab   = "Frequency")               

# Extract Procrustes coordinates 
  coords_all_gpa <- gpa_all$coords   # 34 x 3 x 256
  
#Procrustes ANOVA: procD.lm
procr_anova <- procD.lm(coords_all_gpa ~
                          id +                     # biological variation
                          id:foto +                # error within each photoset
                          id:foto:replicate,       # replicate error
                        data = df,                 # the residuals are the ts!!!
                        iter = 10000, turbo = TRUE)# permutation tests (I chose 1000 permutations)
summary(procr_anova)





#---------------------------
#    LOG NESTED ANOVA
#---------------------------

##### LOG CS

###  1. I need to combine all the CS into one data frame -> long format
n_log <- length(CS_foto1_ts1_log)   # Make sure n is correct -> should be 33

#Combine CS into a matrix: rows = skulls, columns = 8 conditions
CS_matrix_log <-  cbind(CS_foto1_ts1_log,
                        CS_foto1_ts2_log,
                        CS_foto1bis_ts1_log,
                        CS_foto1bis_ts2_log,
                        CS_foto2_ts1_log,
                        CS_foto2_ts2_log,
                        CS_foto2bis_ts1_log,
                        CS_foto2bis_ts2_log)   #dim(CS_mat) should be 33 x 8

# I amke it as a vector -> for each skull I put its 8 replicates in sequence
Csize_log <- as.vector(t(CS_matrix_log))   # length = 8xn = 264

# Factors in the same order
#id= individual
id <- factor(rep(ID_foto, each = 8))   # 33 skulls x 8 = 264 rows
#foto= the photosets
foto <- factor(rep(c("foto1", "foto1", "foto1", "foto1",
                     "foto2", "foto2", "foto2", "foto2"),
                   times = n))
#replicate = original or fotobis
replicate <- factor(rep(c("orig", "orig", "bis", "bis",
                          "orig", "orig", "bis", "bis"),
                        times = n_log))
#ts= training set
ts <- factor(rep(c("ts1", "ts2", "ts1", "ts2", "ts1", "ts2", "ts1", "ts2"),
                 times = n_log))
#build the dataframe
df_log <- data.frame(Csize_log = Csize_log,
                     id        = id,
                     foto      = foto,
                     replicate = replicate,
                     ts        = ts)


# Check counts per id -> they should all be 8
table(df_log$id)
# Check design per skull
with(df_log, table(id, foto, replicate, ts))
# Check that each id has exactly one observation per combination
with(df_log, table(id, foto, replicate, ts))

#Check some statistics :)
tapply(df_log$Csize_log, df_log$id, mean)       # mean Csize per skull
tapply(df_log$Csize_log, df_log$id, sd)         # sd per skull
summary(tapply(df_log$Csize_log, df_log$id, mean))


#make sure they are factor
df_log$id        <- factor(df_log$id)
df_log$foto      <- factor(df_log$foto)      # "foto1", "foto2"
df_log$replicate <- factor(df_log$replicate) # "orig", "bis"
df_log$ts        <- factor(df_log$ts)        # "ts1", "ts2"

#nested ANOVA
nested_anova_log <- aov(Csize_log ~ 1 +                         # overall mean only (no fixed effects)
                          Error(id/foto/replicate/ts),  # random nesting structure
                        data = df)
aov_tbl_log <- summary(nested_anova_log) 


##### get F-value and p-values

## FUNCTION to   extract the residual row (df and MS) from a table
get_residual <- function(tbl) {     # each tbl is a data.frame; the row we need is called "Residuals"
  tbl["Residuals", c("Df", "Sum Sq", "Mean Sq")]
}
# Individual 
df_id_tbl_log<- aov_tbl_log[[1]][[1]]           # the data.frame inside the list
res_id_log   <- get_residual(df_id_tbl_log)
SS_id_log    <- as.numeric(res_id_log["Sum Sq"])
df_id_log    <- as.numeric(res_id_log["Df"])
MS_id_log    <- as.numeric(res_id_log["Mean Sq"])
# Photoset
df_foto_tbl_log <- aov_tbl_log[[2]][[1]]
res_foto_log    <- get_residual(df_foto_tbl_log)
SS_foto_log     <- as.numeric(res_foto_log["Sum Sq"])
df_foto_log     <- as.numeric(res_foto_log["Df"])
MS_foto_log     <- as.numeric(res_foto_log["Mean Sq"])
# Replicate 
df_rep_tbl_log <- aov_tbl_log[[3]][[1]]
res_rep_log    <- get_residual(df_rep_tbl_log)
SS_rep_log     <- as.numeric(res_rep_log["Sum Sq"])
df_rep_log     <- as.numeric(res_rep_log["Df"])
MS_rep_log     <- as.numeric(res_rep_log["Mean Sq"])
# Training set 
df_ts_tbl_log <- aov_tbl_log[[4]][[1]]
res_ts_log    <- get_residual(df_ts_tbl_log)
SS_ts_log     <- as.numeric(res_ts_log["Sum Sq"])
df_ts_log     <- as.numeric(res_ts_log["Df"])
MS_ts_log     <- as.numeric(res_ts_log["Mean Sq"])

##### %SS
total_SS_log <- SS_id_log + SS_foto_log + SS_rep_log + SS_ts_log 
pct_id_log   <- 100 * SS_id_log   / total_SS_log
pct_foto_log <- 100 * SS_foto_log / total_SS_log
pct_rep_log  <- 100 * SS_rep_log  / total_SS_log
pct_ts_log   <- 100 * SS_ts_log   / total_SS_log

##### F and p-values
# Individual vs. Photoset
F_id_log   <- MS_id_log   / MS_foto_log
p_id_log   <- pf(F_id_log, df1 = df_id_log,   df2 = df_foto_log, lower.tail = FALSE)
# Photoset vs. Replicate
F_foto_log <- MS_foto_log / MS_rep_log
p_foto_log <- pf(F_foto_log, df1 = df_foto_log, df2 = df_rep_log, lower.tail = FALSE)
# Replicate vs. Training set
F_rep_log  <- MS_rep_log  / MS_ts_log
p_rep_log  <- pf(F_rep_log, df1 = df_rep_log, df2 = df_ts_log, lower.tail = FALSE)

##final table
final_ANOVA_log <- data.frame(
  Effect    = c("Individual", "Photoset", "Replicate","Training set"),
  SS        = c(SS_id_log,   SS_foto_log,   SS_rep_log,   SS_ts_log),
  `%SS`     = c(pct_id_log,  pct_foto_log,  pct_rep_log,  pct_ts_log),
  MS        = c(MS_id_log,   MS_foto_log,   MS_rep_log,   MS_ts_log),
  df        = c(df_id_log,   df_foto_log,   df_rep_log,   df_ts_log),
  F         = c(F_id_log,    F_foto_log,    F_rep_log,    NA),
  `Pr(>F)`  = c(p_id_log,    p_foto_log,    p_rep_log,    NA)
)
final_ANOVA_log

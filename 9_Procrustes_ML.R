##########################################################
##   Automatizing phenomic research: ALPHA3D pipeline   ##
##                                                      ##
##         PROCRUSTES ANOVA on MANUAL LANDMARKS         ## 
##                                                      ##
##                  Irene Zanandrea                     ##
##########################################################


## LOAD packages, functions, landmarks, lists and dataframe
source("./scripts/3_variance_ML_vs_AL_ log.R")
#load library
library(abind)


# Extract the ID from my dataframe
ID_foto <- log_manual1_distances_df$ID
ID_foto <- gsub("-manual1", "", ID_foto)
ID_foto

####  Manual 1  -----
n_spec_ml1 <- length(manual1_landmarks_array)
coord_array_ml1 <- array(data = unlist(manual1_landmarks_array),
                         dim  = c(34, 3, n_spec_ml1),
                         dimnames = list(landmark = NULL,
                                         axis= c("x", "y", "z"),
                                         specimen = ID_foto))
dim(coord_array_ml1) # should be 34x3x35 (p×k×n_landmark array)
#calculate CS for GPA (General Procrustes Analysis)
Y.gpa.ml1 <- gpagen(coord_array_ml1)
CS_ml1    <- Y.gpa.ml1$Csize
CS_ml1

####  Manual 2  -----
n_spec_ml2 <- length(manual2_landmarks_array)
coord_array_ml2 <- array(data = unlist(manual2_landmarks_array),
                         dim  = c(34, 3, n_spec_ml2),
                         dimnames = list(landmark = NULL,
                                         axis     = c("x", "y", "z"),
                                         specimen = ID_foto))
dim(coord_array_ml2) # should be 34x3x35 (p×k×n_landmark array)
Y.gpa.ml2 <- gpagen(coord_array_ml2)
CS_ml2    <- Y.gpa.ml2$Csize
CS_ml2

#Combine all the coordinates into one big 3D array (34×3×256)
coords_all <- abind(coord_array_ml1, coord_array_ml2, along = 3)
dim(coords_all) # should be: 34 3 70


###### Generalized Procrustes Analysis (GPA)  ---
#library(geomorph)

gpa_all <- gpagen(coords_all)  #gpagen takes landmarks and performs a GPA
str(gpa_all)  # $coords:Procrustes coordinates (shape).

# SUPER COOL 3D visual check of the superimposition
plot(gpa_all)             

Csize   <- gpa_all$Csize
hist(Csize,# distribution of centroid size
     breaks = 16,                     # più barre
     col    = "lightblue",
     border = "white",
     main   = "Distribtion of Centroid size",
     xlab   = "Centroid size",
     ylab   = "Frequency")           



##  Combine all the CS into one data frame -> long format
n <- length(CS_ml1) 
#Combine CS into a matrix: rows = skulls, columns = 2 conditions
CS_matrix<- cbind(CS_ml1,CS_ml2)   #dim(CS_matrix) should be 33 x 2

# I make it as a vector -> for each skull I put its 2 replicates in sequence
Csize <- as.vector(t(CS_matrix))   # length = 2xn = 264


# Factors in the same order
#id= individual
id <- factor(rep(ID_foto, each = 2))   # 33 skulls x 2 = 264 rows
#maual= the different times mML
manual <- factor(rep(c("ml1", "ml2"),
                   times = n))
#build the dataframe
df <- data.frame(Csize = Csize,
                 id    = id,
                 manual= manual)


# Extract Procrustes coordinates 
coords_all_gpa <- gpa_all$coords   # 34 x 3 x 256

df$id        <- as.factor(df$id)
df$manual    <- as.factor(df$manual)

#Procrustes ANOVA: procD.lm
procr_anova <- procD.lm(coords_all_gpa ~
                          id +             # biological variation
                          id:manual,       # manual error
                        data = df,        
                        iter = 10000, turbo = TRUE)# permutation tests 
summary(procr_anova)


procr_anova_1 <- procD.lm(coords_all_gpa ~
                          id ,       # manual error
                        data = df,        
                        iter = 10000, turbo = TRUE)# permutation tests 
summary(procr_anova_1)

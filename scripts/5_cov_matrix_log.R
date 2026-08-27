##########################################################
##   Automatizing phenomic research: ALPHA3D pipeline   ##
##                                                      ##
##                COVARIANCE MATRIX                     ## 
##                                                      ##
##                  Irene Zanandrea                     ##
##########################################################


## load MANUAL AND AUTOMATED DATASET
source("./scripts/libraries.R")
source("./scripts/landmarks_traits.R")
source("./scripts/functions_accuracy.R")
source("./scripts/create_lists_and_dataframe.R")
source("./scripts/2_bias_log.R")



## All the the dataframes are already longtransformed 
## Reorder the traits 
new_order_traits <- c('ISPM_mean','ISPNS','ISNSL','PMZS_mean','PMZI_mean','PMMT_mean','NSLNA','NAPNS',  
                  'NSLZS_mean','NSLZI_mean','MTPNS_mean','ZSZI_mean','ZIMT_mean','ZIZYGO_mean','ZITSP_mean',
                  'EAMZYGO_mean','ZYGOTSP_mean','PTZYGO_mean','PTTSP_mean','NABR', 'BRPT_mean','BRAPET_mean', 
                  'PTAPET_mean','PTBA_mean','PTEAM_mean','LDAS_mean','BRLD','OPILD','PTAS_mean','PNSAPET_mean', 
                  'APETBA_mean','APETTS_mean','BAEAM_mean', 'JPAS_mean','BAOPI')

# Convert Trait to a factor with the specified levels
new_log_auto_distances_df <- log_auto_distances_df%>% select(all_of(new_order_traits)) 
new_log_manual_distances_df <- log_manual_distances_df%>% select(all_of(new_order_traits))  



#-------------------------------------
# COVARIANCE MATRIX
#-------------------------------------

####### MANUAL
manual.cov = cov(new_log_manual_distances_df)
manual.cov
    
isSymmetric.matrix(manual.cov) #TRUE
manual.variances = diag(manual.cov)
manual.variances
    
manual.eigen <- eigen(manual.cov)$values
man.eig.per = manual.eigen/sum(manual.eigen)
pc.m <- c(2:36)
man.eig.per.mat <- cbind(pc.m, man.eig.per)
    
man.pc.dist <- ggplot(man.eig.per.mat,aes(x = pc.m, y = man.eig.per)) +
                      geom_line() +
                      geom_point() +
                      scale_y_continuous("%Variation in the PC") +
                      scale_x_discrete("Principal component rank") +
                      theme(text = element_text(size = 16),
                            legend.position = "none",
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(), 
                            axis.line = element_line(colour = "black"),
                            plot.background = element_rect(fill = 'transparent', color = NA)) +
                      theme_linedraw()+ 
                      labs(title = "Manual - all trait")+
                      ylim(0,0.85)
man.pc.dist
    
####### AUTOMATIC
auto.cov = cov(new_log_auto_distances_df)
auto.cov
    
isSymmetric.matrix(auto.cov) #TRUE
auto.variances = diag(auto.cov)
auto.variances 
    
auto.eigen <- eigen(auto.cov)$values
aut.eig.per = auto.eigen/sum(auto.eigen)
pc <- c(2:36)
aut.eig.per.mat <- cbind(pc, aut.eig.per)
    
aut.pc.dist <- ggplot(aut.eig.per.mat,aes(x = pc, y = aut.eig.per)) +
                      geom_line() +
                      geom_point() +
                      scale_y_continuous("%Variation in the PC") +
                      scale_x_discrete("Principal component rank") +
                      theme(text = element_text(size = 16),
                            legend.position = "none",
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(), 
                            axis.line = element_line(colour = "black"),
                            plot.background = element_rect(fill = 'transparent', color = NA)) +
                      theme_linedraw() + 
                      labs(title = "Automated - all trait")+
                      ylim(0,0.85)
aut.pc.dist



#------------------------------------------------------
# CHOOSE PRECISE TRAITS TO BUILD COVARIANCE MATRICES
#------------------------------------------------------  

## 4 traits VERY PRECISE
   ### MANUAL 
   manual.cov.4b = cov(new_log_manual_distances_df[,c(2,8,17,33)]) #4b=4 best
   isSymmetric.matrix(manual.cov.4b) #TRUE          
   ### AUTOMATIC                    
   auto.cov.4b = cov(new_log_auto_distances_df[,c(2,8,17,33)]) 
   isSymmetric.matrix(auto.cov.4b) #TRUE    
        
## 4 traits NO PRECISE
   ### MANUAL  
   manual.cov.4w = cov(new_log_manual_distances_df[,c(7,19,26,34)]) #4w=4 worst
   isSymmetric.matrix(manual.cov.4w)
   ### AUTOMATIC
   auto.cov.4w = cov(new_log_auto_distances_df[,c(7,19,26,34)]) 
   isSymmetric.matrix(auto.cov.4w) #TRUE
        
## 8 traits VERY PRECISE
   ### MANUAL
   manual.cov.8b = cov(new_log_manual_distances_df[,c(2,4,5,8,9,10,17,33)]) 
   isSymmetric.matrix(manual.cov.8b) #TRUE      
   ### AUTOMATIC
   auto.cov.8b = cov(new_log_auto_distances_df[,c(2,4,5,8,9,10,17,33)])
   isSymmetric.matrix(auto.cov.8b) #TRUE
        
## 8 traits NO accurate
   ### MANUAL
   manual.cov.8w = cov(new_log_manual_distances_df[,c(1,7,19,20,26,27,34,35)]) 
   isSymmetric.matrix(manual.cov.8w)
   ### AUTOMATIC
   auto.cov.8w = cov(new_log_auto_distances_df[,c(1,7,19,20,26,27,34,35)])
   isSymmetric.matrix(auto.cov.8w) #TRUE
        
        
## 12 traits VERY accurate
   ### MANUAL 
   manual.cov.12 = cov(new_log_manual_distances_df[,c(2,5,6,8,10,14,17,18,21,
                                                      29,30,33)])  
   isSymmetric.matrix(manual.cov.12) #TRUE             
   ### AUTOMATIC                           
   auto.cov.12 = cov(new_log_auto_distances_df[,c(2,5,6,8,10,14,17,18,21,29,
                                                  30,33)])
   isSymmetric.matrix(auto.cov.12) #TRUE
        
## 16 traits VERY accurate
   ### MANUAL 
   manual.cov.16 = cov(new_log_manual_distances_df[,c(2,4,5,6,8,9,10,14,15,
                                                      17,18,21,24,29,30,33)])    
   isSymmetric.matrix(manual.cov.16) #TRUE             
   ### AUTOMATIC                           
   auto.cov.16 = cov(new_log_auto_distances_df[,c(2,4,5,6,8,9,10,14,15,17,18,
                                                  21,24,29,30,33)])
   isSymmetric.matrix(auto.cov.16) #TRUE
        
## 20 traits VERY accurate
   ### MANUAL 
   manual.cov.20 = cov(new_log_manual_distances_df[,c(2,3,4,5,6,8,9,10,12,14,
                                                15,17,18,21,22,24,28,29,30,33)])    
   isSymmetric.matrix(manual.cov.20) #TRUE            
   ### AUTOMATIC                          
   auto.cov.20 = cov(new_log_auto_distances_df[,c(2,3,4,5,6,8,9,10,12,14,
                                              15,17,18,21,22,24,28,29,30,33)])
   isSymmetric.matrix(auto.cov.20) #TRUE
        
## 24 traits VERY accurate
   ### MANUAL 
   manual.cov.24 = cov(new_log_manual_distances_df[,c(2,3,4,5,6,8,9,10,11,
                                12,14,15,16,17,18,21,22,24,25,28,29,30,31,33)])    
   isSymmetric.matrix(manual.cov.24) #TRUE            
   ### AUTOMATIC                          
   auto.cov.24 = cov(new_log_auto_distances_df[,c(2,3,4,5,6,8,9,10,11,12,14,
                                      15,16,17,18,21,22,24,25,28,29,30,31,33)])
   isSymmetric.matrix(auto.cov.24) #TRUE
        
## 28 traits VERY accurate
   ### MANUAL 
   manual.cov.28 = cov(new_log_manual_distances_df[,c(2,3,4,5,6,8,9,10,11,12,13,
                          14,15,16,17,18,21,22,23,24,25,26,28,29,30,31,32,33)])    
   isSymmetric.matrix(manual.cov.28) #TRUE            
   ### AUTOMATIC                          
   auto.cov.28 = cov(new_log_auto_distances_df[,c(2,3,4,5,6,8,9,10,11,12,13,14,15,
                                16,17,18,21,22,23,24,25,26,28,29,30,31,32,33)])
   isSymmetric.matrix(auto.cov.28) #TRUE
        
## 32 traits VERY accurate -> all but 19,20,34,35
   ### MANUAL 
   manual.cov.32 = cov(new_log_manual_distances_df[,c(1,2,3,4,5,6,7,8,9,10,11,12,
                        13,14,15,16,17,18,21,22,23,24,25,26,27,28,29,30,31,32,33)])    
   isSymmetric.matrix(manual.cov.32) #TRUE            
   ### AUTOMATIC                          
   auto.cov.32 = cov(new_log_auto_distances_df[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,
                          15,16,17,18,21,22,23,24,25,26,27,28,29,30,31,32,33)])
   isSymmetric.matrix(auto.cov.32) #TRUE
        
## 34 traits VERY accurate -> all but 19
   ### MANUAL
   manual.cov.34 = cov(new_log_manual_distances_df[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,
                      14,15,16,17,18,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)])    
   isSymmetric.matrix(manual.cov.34) #TRUE            
   ### AUTOMATIC                          
   auto.cov.34 = cov(new_log_auto_distances_df[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
                          16,17,18,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)])
        isSymmetric.matrix(auto.cov.34) #TRUE
       
## 34 bis traits VERY accurate -> all but 20
   ### MANUAL
   manual.cov.34bis = cov(new_log_manual_distances_df[,c(1,2,3,4,5,6,7,8,9,10,11,12,
                    13,14,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)])    
   isSymmetric.matrix(manual.cov.34bis) #TRUE            
   ### AUTOMATIC                          
   auto.cov.34bis = cov(new_log_auto_distances_df[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,
                    14,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)])
   isSymmetric.matrix(auto.cov.34bis) #TRUE
        
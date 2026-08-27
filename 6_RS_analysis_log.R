##########################################################
##   Automatizing phenomic research: ALPHA3D pipeline   ##
##                                                      ##
##             RANDOM SKEWERS ANALYSIS                  ## 
##                                                      ##
##                  Irene Zanandrea                     ##
##########################################################


## load COVARIANCE MATRIX
source("./scripts/5_cov_matrix_log.R")

# auto.cov   -> automated P matrix (log-transform)
# manual.cov -> manual P matrix (log-transform)

#### evolqg package
RandomSkewers(auto.cov, manual.cov, num.vectors = 1000000,  # ALL DATASET
              parallel = FALSE)     

##based on PRECISION
RandomSkewers(auto.cov.4b, manual.cov.4b, num.vectors = 1000000, # 4 BEST traits
              parallel = FALSE) 
RandomSkewers(auto.cov.4w, manual.cov.4w, num.vectors = 1000000, # 4 WORST traits 
              parallel = FALSE)

RandomSkewers(auto.cov.8b, manual.cov.8b, num.vectors = 1000000, # 8 BEST traits
              parallel = FALSE)  
RandomSkewers(auto.cov.8w, manual.cov.8w, num.vectors = 1000000, # 8 WORST traits
              parallel = FALSE)  

RandomSkewers(auto.cov.12, manual.cov.12, num.vectors = 1000000, # 12 HIGH precision
              parallel = FALSE)
RandomSkewers(auto.cov.16, manual.cov.16, num.vectors = 1000000, # 16 HIGH precision
              parallel = FALSE)
RandomSkewers(auto.cov.20, manual.cov.20, num.vectors = 1000000, # 20 HIGH precision
              parallel = FALSE)
RandomSkewers(auto.cov.24, manual.cov.24, num.vectors = 1000000, # 24 HIGH precision
              parallel = FALSE)
RandomSkewers(auto.cov.28, manual.cov.28, num.vectors = 1000000, # 28 HIGH precision
              parallel = FALSE)
RandomSkewers(auto.cov.32, manual.cov.32, num.vectors = 1000000, # 32 HIGH precision
              parallel = FALSE)
RandomSkewers(auto.cov.34, manual.cov.34, num.vectors = 1000000, # 34 HIGH precision
              parallel = FALSE)

RandomSkewers(auto.cov, manual.cov, num.vectors = 1000000,       # ALL DATASET
              parallel = FALSE)   

                         
##########################################################
##   Automatizing phenomic research: ALPHA3D pipeline   ##
##                                                      ##
##                     MESHDIST                         ## 
##                                                      ##
##                  Irene Zanandrea                     ##
##########################################################


# install packages
#install.packages("rgl")
#install.packages("Morpho")
#install.packages("Rvcg")

# load packages
library(rgl)
library(Morpho)
library(Rvcg)


###### foto1 vs foto2
## 1. Read two meshes (PLY)
mesh1 <- vcgImport("./meshdist/42222_foto1.ply")  
mesh2 <- vcgImport("./meshdist/42222_foto2.ply")
#check how big they are
ncol(mesh1$it) 
ncol(mesh2$it)
#downsize them 
mesh1_simpl <- vcgQEdecim(mesh1, percent = 0.3)
mesh2_simpl <- vcgQEdecim(mesh2, percent = 0.3)

## 2. Compute the distance between meshes
## meshDist computes the average minimal distance between triangles
dist_res <- meshDist(mesh1_simpl, mesh2_simpl)


###### foto1 vs fotobis
mesh1bis <- vcgImport("./meshdist/42222_foto1bis.ply")
mesh1bis_simpl <- vcgQEdecim(mesh1bis, percent = 0.3)
dist_res_1 <- meshDist(mesh1_simpl, mesh1bis_simpl)

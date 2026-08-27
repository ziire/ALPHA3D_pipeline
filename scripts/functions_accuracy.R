##########################################################
##   Automatizing phenomic research: ALPHA3D pipeline   ##
##                                                      ##
##      FUNCTIONS FOR PRECISION, BIAS AND ACCURACY      ## 
##                                                      ##
##                  Irene Zanandrea                     ##
##########################################################


#---------------------------------------
# FUNCTION TO PROCESS ALPACA FILES
#---------------------------------------

## Process json file function
process_json <- function(file.json,   #output file after ALPACA -> in the input for the function you have to put the full path                  
                         landmarks_n, #number of landmarks(in our case it`s 34)                               
                         coord_n)     #number of coordinates to extract (3 because we are in 3D)                                
 {data <- fromJSON(file.json)                    #fromJSON = from library(jsonlite)                     
  markups_list <- data$markups                   #extract markups inside jsonfile from list                     
  control_points <- markups_list$controlPoints   #extract controlPoints (inside JSONfile)       
  landmarks <- control_points[[1]]               #inside controlPoints 1^element landmarks        
  landmarks_position <- landmarks$position       #extract position            
  #create matrix
  landmarks_matrix <- do.call(rbind,             #do.call unpacks the list of objects inside landmarks_position
                              landmarks_position)#rbind binds the elements (position) row-wise 
  #create array (contains lots of matrices)
  #we need to create an array because most geomorph functions expect the landmark data as a 3D array.
  landmarks_array <- arrayspecs(            #arrayspecs from geomorph convert a matrix of landmark coordinates into a three-dimensional array
                          landmarks_matrix, #landmarks_matrix,          
                          landmarks_n,      #rows:number of landmarks(landmarks_n)(34)          
                          coord_n,          #columns:number of coordinates(3for3D)
                          sep = NULL)    
  landmarks_array}                                        



#------------------------------
# FUNCTION FOR DISTANCES
#------------------------------

## Euclidean distance
euclidean_distance <- function(x,      #array (array contains more matrices)
                               point1, #index of points from which to 
                               point2)     #calculate euclidean distance
{sum((x[point1, , 1] - x[point2, , 1])^2) ^ 0.5 } 
  #difference between coordinates point1 and point2 
  #^2 = squared. ^ 0.5 = sqaured root. 



#--------------------------------------------------
# FUNCTION to take MEANS between BILATERAL TRAITS
#--------------------------------------------------

## MEANS between columns of BILATERAL TRAITS 
avg_columns <- function(df)            # function takes dataframe as argument
{ id_col <- df[, ncol(df)]             # 1.Save and remove last column (ID) 
  df <- df[,-ncol(df)]                 # Remove the last column 
  avg_list <- list()                   # 2.Initialize an empty list to store the averages 
  for (i in seq(1, ncol(df), by = 2))  # 3.Loop over column in pairs starting from 1st column 
    {mean_values <- rowMeans(df[, c(i, i + 1)], #Mean between rows for 2 columns  
                              na.rm = TRUE) 
      avg_list[[length(avg_list) + 1]] <- mean_values #assign the means to avg_list
    }
  avg_df <- as.data.frame(avg_list)    # 4. Convert average list into a data frame 
  avg_df <- cbind(id_col, avg_df)      # cbind add back last column as ID of avg_df
  return(avg_df)                       # Gives back final data frame
}


#------------------------------
# FUNCTION for IMPRECISION
#------------------------------

## Variance (spread of d1-d2)
variance_formula <- function(d1,  #measure from dataset1
                             d2)  #measures dataset2
{ ((d1 - d2)^2) / 2}              #average squared of d1 and d2

## Calculate imprecision 
calculate_imprecision <- function(data1,      #data1 = 1st dataframe (foto1_distances_df)  
                                  data2,      #data2 = 2nd dataframe (foto2_distances_df)
                                  id_col,     #id_col= name of ID column (ID),
                                  new_suffix) #new suffix for ID (-relative error)
{ # 1.Create a df without ID  
  imprecision_df <- data1 %>%              
    select(-!!sym(id_col))      #we drop a column named in the variable id_col (name of ID column).    
  # 2.start loop for every column of df
  for (col_name in colnames(imprecision_df))
    { impr_value <- variance_formula(                  #apply imprecision formula to each column
                    imprecision_df[[col_name]],        #d1= Current column from data1
                    data2[[col_name]])                 #d2= Corresponding column from data2
      imprecision_df[[paste0("imprecision_",           #add at the column name imprecision_
                             col_name)]] <- impr_value}#assign impr_value in a list
  imprecision_df <- imprecision_df %>%    # 3.Add back the modified ID column                       
      mutate(ID = paste0(data1[[id_col]], # mutate(): modifies existing columns
                         new_suffix))     # ID = paste0() create new ID column with suffix
  total_imprecision <- imprecision_df %>%    # 4. Select new column for final output 
      select(ID, starts_with("imprecision_"))# Select ID + columns starting with imprecision_
  return(total_imprecision)}                  # Gives back final df total_imprecision  
                            


#-------------------------
# FUNCTION for BIAS
#-------------------------

# Bias= d'-d
bias_function <- function(df1, #dataframe1
                          df2, 
                          id_col=1)
{cbind(df1[id_col],                        #unite difference and ID column
       df1[ , -id_col] - df2[ , -id_col])} #calculate differences between the 2dfs
                                           #excluding the first element (ID)


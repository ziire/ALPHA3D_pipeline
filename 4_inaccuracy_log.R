##########################################################
##   Automatizing phenomic research: ALPHA3D pipeline   ##
##                                                      ##
##                    INACCURACY                        ## 
##                                                      ##
##                  Irene Zanandrea                     ##
##########################################################


## load imprecision and bias
source("./scripts/1_imprecision_log.R")
source("./scripts/2_bias_log")



#-------------------------------------
# COMPUTE INACCURACY
#-------------------------------------

### Keep only bias^2
bias_long <- subset(bias_long,select = -c(Difference))

##merge - keep all rows from imprecision_df; the unmatched columns become NA
inaccuracy_df <- merge(log_imprecision_total_long, 
                       bias_long, 
                       by = c ("Trait", "ID"), 
                       all = TRUE) 
##  Inaccuracy = Imprecision + Bias^2
inaccuracy_df_tot <- inaccuracy_df
inaccuracy_df_tot$Inaccuracy <- inaccuracy_df_tot$Bias_2 + 
                                inaccuracy_df_tot$Variance_total

# Create a new column for ordinal numbers from 1 to n
inaccuracy_df_tot$TraitOrdinal <- as.numeric(inaccuracy_df_tot$Trait)
inaccuracy_df_tot$TraitOrdinal <- factor(inaccuracy_df_tot$TraitOrdinal)
   

#-------------------------------------
# MEAN AND MEDIAN of inaccuracy
#-------------------------------------

### TOTAL
mean_inaccuracy_df_tot <- mean(inaccuracy_df_tot$Inaccuracy, na.rm = TRUE)
median_inaccuracy_df_tot <- median(inaccuracy_df_tot$Inaccuracy, na.rm = TRUE)
  ##### Table 3. Mean and median values of inaccuracy
  inaccuracy_tot_mean_median_df  <- data.frame(
                          mean_inaccuracy_total = mean_inaccuracy_df_tot,
                          median_inaccuracy_total= median_inaccuracy_df_tot)
### PER TRAIT
mean_inaccuracy <- inaccuracy_df_tot %>%  
  group_by(Trait) %>%
  summarise(Mean_Inaccuracy = round(mean(Inaccuracy, na.rm = TRUE),4))

median_inaccuracy <- inaccuracy_df_tot %>%  
  group_by(Trait) %>%
  summarise(Median_Inaccuracy = round(median(Inaccuracy, na.rm = TRUE),4))
    
  #### Table S7: Comparison of MEAN for each trait.
  pertrait_inaccuracy <- merge(
                         mean_inaccuracy, median_inaccuracy, 
                          by = c ("Trait"), all = TRUE)
  pertrait_inaccuracy$Trait <- factor(pertrait_inaccuracy$Trait, 
                                      levels = order_traits)  # Chose the trait order
  pertrait_inaccuracy <- pertrait_inaccuracy %>% arrange(Trait)
# Save dataset 
#write.csv(pertrait_inaccuracy, "pertrait_inaccuracy_log.csv", row.names = FALSE)



#-----------------
#   PLOT
#-----------------
  
# 1- imprecision  
ggplot(inaccuracy_df_tot, aes(x = TraitOrdinal, y = Variance_total)) +  
       geom_boxplot() +  
       stat_summary(fun = median, geom = "point", shape = 20, size = 3, 
                    color = "mediumvioletred", fill = "mediumvioletred") +  
        labs(title = "Imprecision of the Linear Measuraments", 
             x = "Linear Measuraments",y = "Logged Variance (%)") +  
        theme_minimal() +  
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.6),
              panel.grid.major.x = element_blank(),  # Remove major vertical grid lines
              axis.text.x = element_text(angle = 45, hjust = 1, size = 24),
              axis.text.y = element_text(hjust = 1, size = 24),
              axis.title.x = element_text(size = 26),
              axis.title.y = element_text(size = 26)) 

# 2- bias squared
ggplot(inaccuracy_df_tot, aes(x = TraitOrdinal, y = Bias_2)) +
        geom_boxplot() +
        stat_summary(fun = median, geom = "point", shape = 20, 
                     size = 3, color = "red", fill = "red") + 
        labs(title = "Bias of AL relative to ML", 
             x = "Linear Measurements", y = "Logged Bias Squared (%)")  +
        theme_minimal() + 
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
              panel.grid.major.x = element_blank(),  # Remove major vertical grid lines
              axis.text.x = element_text(angle = 45, hjust = 1, size = 27),
              axis.text.y = element_text(hjust = 1, size = 25),
              axis.title.x = element_text(size = 23),
              axis.title.y = element_text(size = 23),
              plot.title = element_text(size = 25, face = "bold")) +
       
# 3- Inaccuracy
ggplot(inaccuracy_df_tot, aes(x = TraitOrdinal, y = Inaccuracy)) +
        geom_boxplot(outlier.alpha = 0.5) +
        stat_summary(fun = median, geom = "point", shape = 21, 
                     size = 2,color = "darkmagenta", fill = "darkmagenta") +
        labs(title = "Inaccuracy of the Linear Measuraments",
             x = "Linear Measuraments",
             y = "Logged Inaccuracy (%)")  +
        theme_minimal() + 
        theme(legend.position = "bottom",
              legend.title = element_text(size = 15, face = "bold"),
              legend.text  = element_text(size = 18),
              panel.border = element_rect(color = "black", fill = NA, size = 0.5),
              panel.grid.major.x = element_blank(),
              axis.text.x  = element_text(angle = 45, hjust = 1, size = 20),
              axis.text.y  = element_text(hjust = 1, size = 19),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20),
              plot.title   = element_text(size = 15, face = "bold")) +
        guides(color = "none")  +
       

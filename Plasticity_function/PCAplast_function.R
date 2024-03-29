################################################################################
################### Sample plasticity based on PC distances ####################
############## Written by Colleen B. Bove (colleenbove@gmail.com) ##############
########################## Laste update: 10 Jan 2023 ###########################
################################################################################

## Required packages to run:
# dplyr


## PC distances can be calculated from prcomp() objects or from plotPCA()


## To run the function, enter the following objects:
# PCAplast(pca = XXX, # the PCA dataframe containing the PCA eigenvalues
#          data = XXX, # the condition/treatment data corresponding to samples
#          sample_ID = "XXX", # the name of column that provide unique ID per sample (if blank, will pull rownames for this)
#          num_pca =  "XXX", # the number of PCAs to include in analysis (default is 'all', but you can specify another number with a minimum of 2 PCAs)
#          control_col = "XXX", # what the 'treatment' column is called
#          control_lvl = "XXX", # control level of the treatment. If blank, a control mean per control level is assumed
#          group = "XXX") # the grouping column (i.e., colony). If blank, will assume control level grouping only!



########################################################################################

PCAplast <- function(pca, data, sample_ID = NA, num_pca = "all", control_col, control_lvl = "none", group = NA) {
  
  # rename the user input info
  pca_df <- pca
  data_df <- data
  control_name <- control_col
  control_lvl <- control_lvl
  group_col <- group
  
  # ifelse statement to pull PCs from dataframe or prcomp objects
  if(class(pca_df) == "prcomp"){
    pca_dist <- pca_df$x # grab PC distances from prcomp() object
  } else {
    pca_dist <- pca_df # grab PC distances from data.frame object
  }
  
  
  # check for correct number of PCAs provided
  if(class(num_pca) == "numeric") {
    if(num_pca < 2) { # will throw error if too few PCAs requested
      stop("please select more than 2 PCs to calculate distance")
    } 
  }
  
  if(class(num_pca) == "numeric") {
    if(num_pca > (data.frame(pca_dist) %>% dplyr::select(starts_with("PC")) %>% ncol())) { # will throw error if too many PCAs requested
      stop(paste(num_pca, "PCs requested for calculation, but only", (pca_dist %>% dplyr::select(starts_with("PC")) %>% ncol()), "PCs available. Select appropriate number of PCAs for calculation."))
    } 
  }
  
  
  # oder the dataframe to ensure correct pairing after calculating distance 
  if(sample_ID %in% colnames(data_df)){
    data_df <- data_df[match(row.names(pca_dist), data_df[[sample_ID]]),]
  } else {
    data_df <- data_df[order(as.numeric(row.names(data_df))),]
  }
  
  
  # combine the datasets
  dist_df <- cbind(data_df, pca_dist) 
  
  
  # if there is no control level, modify so function pulls all levels
  if(control_lvl == "none") {
    control_lvl = unique(dist_df[[control_name]])
  } else {
    control_lvl = control_lvl
  }
  
  
  # make dataframe of control grouping only
  if(!is.na(group_col)) { # calculate mean per grouping ID (if provided)
    mean_control <- dist_df %>%
      filter(dist_df[[control_name]] == list(control_lvl)[[1]]) %>% 
      rename_with(tolower) %>% # renames all pc's with lowercase 'PC' (just to differentiate from all sample PCs)
      dplyr::select(colnames((dist_df %>% rename_with(tolower))[tolower(group_col)]), starts_with("pc")) %>%
      group_by_at(vars(tolower(group_col))) %>% 
      summarise_if(is.numeric, mean)
    
    
    # add the control PCA values to treatment samples per grouping
    dist_df2 <- merge(dist_df, mean_control, by.x = group_col, by.y = tolower(group_col))
    
  } else { # calculate mean per control treatment
    mean_control <- dist_df %>%
      filter(dist_df[[control_name]] == list(control_lvl)[[1]]) %>% 
      rename_with(tolower) %>% # renames all pc's with lowercase 'PC' (just to differentiate from all sample PCs)
      dplyr::select(colnames((dist_df %>% rename_with(tolower))[tolower(control_name)]), starts_with("pc")) %>% # select just the PCs 
      #group_by() 
      group_by_at(vars(tolower(control_name))) %>% 
      summarise_if(is.numeric, mean)
    
    # add the control PCA values to all samples 
    #dist_df2 <- merge(dist_df, mean_control, by.x = control_col, by.y = tolower(control_col), all.x = TRUE)
    dist_df2 <- merge(dist_df, mean_control, all = TRUE)
    
  }
  
  
  
  if(sample_ID %in% colnames(data_df)){
    data_df <- data_df[match(row.names(pca_dist), data_df[[sample_ID]]),]
  } else {
    data_df <- data_df[order(as.numeric(row.names(data_df))),]
  }
  
  
  # again, reorder data
  if(sample_ID %in% colnames(data_df)){
    dist_df2 <- dist_df2[order(dist_df2[[sample_ID]]),]
  } else {
    rownames(dist_df2) <- rownames(data_df)
    dist_df2 <- dist_df2[order(as.numeric(row.names(dist_df2))),]
  }
  
  
  ### Calculate sample (PCA) distances from control (pca) using all PCAs
  # make dataframe to populate with pca distances
  full_calc_dist <- data.frame(control_name = dist_df2[control_name])
  
  if(num_pca == "all") {
    ## forloop that will calculate distances between control and sample for all PCs (n will be total number)
    for(n in 1:(dist_df %>% dplyr::select(starts_with("PC")) %>% ncol())){
      # makes the PCA column name for control (lowercase) and sample (uppercase)
      PC_col <- paste0("PC", n)
      pc_col <- paste0("pc", n)
      
      # pulls the PC column for control (lowercase) and sample (uppercase)
      PCx <- dist_df2[PC_col]
      pcx <- dist_df2[pc_col]
      
      pca_calc_dist <- data.frame((PCx - pcx)^2) # calculates the distance between 2 PCs
      full_calc_dist <- cbind(full_calc_dist, pca_calc_dist) # add that distance to running dataframe
    }
  } else {
    ## forloop that will calculate distances between control and sample for SPECIFIED # of PCs (n will be total number)
    for(n in 1:as.numeric(num_pca)){
      # makes the PC column name for control (lowercase) and sample (uppercase)
      PC_col <- paste0("PC", n)
      pc_col <- paste0("pc", n)
      
      # pulls the PC column for control (lowercase) and sample (uppercase)
      PCx <- dist_df2[PC_col]
      pcx <- dist_df2[pc_col]
      
      pca_calc_dist <- data.frame((PCx - pcx)^2) # calculates the distance between 2 PCs
      full_calc_dist <- cbind(full_calc_dist, pca_calc_dist) # add that distance to running dataframe
    }
  }
  
  
  ## final distance calculation (adds all PCA distances and takes squareroot)
  distance <- full_calc_dist %>% 
    mutate(dis_sum = rowSums(across(where(is.numeric)))) %>% 
    mutate(dist = sqrt(dis_sum)) %>% 
    dplyr::select(matches("dist"))
  
  
  ## combine the calculated distance with the metadata and remove controls for final dataframe
  dist_df <- data_df %>% 
    bind_cols(distance) %>% 
    filter(!is.na(dist)) 
  
  ## removes the control levels
  if(length(control_lvl) > 1){
    dist_df <- dist_df
  } else {
    dist_df <- dist_df %>%
      filter(dist_df[[control_name]] != control_lvl)  %>% 
      droplevels()
  }
  
}


#################################################################
###### Session information from last update

# R version 3.6.3 (2020-02-29)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.16
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_3.3.6 dplyr_1.0.8  
# 
# loaded via a namespace (and not attached):
#   [1] rstudioapi_0.13  magrittr_2.0.3   munsell_0.5.0    tidyselect_1.1.1 colorspace_2.0-3 R6_2.5.1        
# [7] rlang_1.0.4      fansi_1.0.3      tools_3.6.3      grid_3.6.3       gtable_0.3.0     xfun_0.29       
# [13] tinytex_0.40     utf8_1.2.2       cli_3.3.0        DBI_1.1.3        withr_2.5.0      ellipsis_0.3.2  
# [19] digest_0.6.29    assertthat_0.2.1 tibble_3.1.7     lifecycle_1.0.1  crayon_1.5.1     farver_2.1.1    
# [25] purrr_0.3.4      vctrs_0.4.1      glue_1.6.2       labeling_0.4.2   compiler_3.6.3   pillar_1.8.0    
# [31] scales_1.2.0     generics_0.1.3   pkgconfig_2.0.3

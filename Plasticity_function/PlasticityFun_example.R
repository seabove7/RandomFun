################################################################################
#################### Example use of PCAplast function use #####################
################### Colleen B. Bove (colleenbove@gmail.com) ####################
################################################################################

## Load packages
# library(devtools) # will need devtools to download ggfortify from GitHub directly 
# install_github('sinhrks/ggfortify')
library(dplyr) # necessary for PCA function
library(ggbiplot) # plotting the PCA
library(ggfortify) # plotting the PCA
library(vegan) # running the PERMANOVA (adonis2())
library(ggpubr) # for arranging multiple plots into single figure
source("Plasticity_function/PCAplast_function.R") # source the plasticity function



## Working with data
data("iris") # Lload the 'iris' dataset
pca_input <- iris[1:4] # select for just the measured parameters to perform PCA 


## Perform principal component analysis (PCA)
pca_df <- prcomp(pca_input, center = TRUE, scale = TRUE)


## plot the PCA coloured by species
pca_plot <- autoplot(pca_df, data = iris, 
         colour = "Species",
         shape = "Species",
         fill = "Species",               
         frame = TRUE, 
         frame.type = "t", # multivariate t-distributions for small n producing heavier tails
         frame.level = 0.95) + # using 95% CI for all ellipses
  theme_classic()
pca_plot


## run the PERMANOVA 
adonis_out <- adonis2(pca_input ~ Species, data = iris, method = 'eu')
adonis_out # view adonis output


## To run the plasticity function, enter the following objects:
plast_out <- PCAplast(pca = pca_df, # the PCA dataframe containing the PCA eigenvalues
         data = iris, # the condition/treatment data corresponding to samples
         #sample_ID = "XXX", # the name of column that provide unique ID per sample (if blank, will pull rownames for this)
         num_pca =  "all", # the number of PCAs to include in analysis (default is 'all', but you can specify another number with a minimum of 2 PCAs)
         #control_lvl = "XXX", # control level of the treatment. If blank, a control mean per control level is assumed
         #group = "XXX", # the grouping column (i.e., colony). If blank, will assume control level grouping only!
         control_col = "Species") # what the 'treatment' column is called


## Plot the plasticity (PC distances): overlay mean and 1 standard deviation
plast_plot <- ggplot(data = plast_out, aes(x = Species, y = dist)) + 
  geom_point(alpha = 0.3, position = position_jitter(width = 0.1)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0, colour = "red") +
  stat_summary(fun = "mean", size = 0.5, colour = "red") +
  theme_classic() +
  ylab("Plasticity")
plast_plot


## Combine PCA and plasticity plot into single figure and save as PNG
ggarrange(pca_plot, plast_plot, labels = "AUTO")
  ggsave("Plasticity_function/sample_plot.png", width = 10, height = 4)

  
  
  
  
  

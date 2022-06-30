# PCAplast function

The following function will calculate principal component distances between samples as a measure of plasticity. This measure was adapted from [Barott et al., 2021](https://www.pnas.org/content/118/22/e2025435118) and was initially developed to calculate physiological plasticity of corals under ocean acidification and warming stressors in a mesocosm experiment ([Bove et al. In Prep](https://www.biorxiv.org/content/10.1101/2021.07.13.452173v1)).

<br/>


### Please cite this function as:
Bove, C. (2022). RandomFun: Plasticity from PC distances (Version 1.0) [Computer software]. DOI

## **Function usage:**

```
plast_out <- PCAplast(pca = XX, # the PCA dataframe containing the PCA eigenvalues
         data = XX, # the condition/treatment data corresponding to samples
         sample_ID = "XXX", # the name of column that provide unique ID per sample (if blank, will pull rownames for this)
         num_pca =  "XXX", # the number of PCAs to include in analysis (default is 'all', but you can specify another number with a minimum of 2 PCAs)
         control_lvl = "XXX", # control level of the treatment. If blank, a control mean per control level is assumed
         group = "XXX", # the grouping column (i.e., colony). If blank, will assume control level grouping only!
         control_col = "XXX") # what the 'treatment' column is called
```         

<br/>


## **Arguments:**

argument | explanation
--- | --- 
pca | PCA dataframe containing the PCA eigenvalues. Can be a prcomp object or data.frame.
data | The condition/treatment data corresponding to samples. This should contain everything you may want to assess these data with downstream.
sample_ID | Name of column that provide unique ID per sample (if this is blank, the function will pull rownames for this)
num_pca | The number of PCs to include in analysis (default is 'all', but you can specify another number with a minimum of 2 PCAs)
control_lvl | Control level of the 'treatment'. If blank, a mean per control level is assumed.
group | The grouping column (i.e., colony). If blank, will assume control level grouping only.
control_col | Column used for the 'treatment' (this can be a grouping by treatment, species, etc.)

<br/>

## **Example:**

An example use of this function can be found on the [PlasticityFun_example.R](https://github.com/seabove7/furry-lamp/blob/main/Plasticity_function/PlasticityFun_example.R) script in this repository using the Iris dataset. 

<p align="center">
<img src="https://github.com/seabove7/furry-lamp/blob/main/Plasticity_function/sample_plot.png" width = "700" />
</p> 

<br/>

---
**Please reach out with questions about the function or if you encounter any bugs!**  
*Email: colleenbove@gmail.com*  
*Website [colleenbove.science](http://colleenbove.science)*  
*Twitter: [@DrSeaBove](https://twitter.com/DrSeaBove)*

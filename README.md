<!-- README.md is generated from README.Rmd. Please edit that file -->
# phenokio
<!-- badges: start -->
<!-- badges: end -->
The goal of the phenokio package is to perform automated grain morphometric analysis from digital images.

□ Code explained: https://agronomy4future.com/archives/25117

## Installation
You can install phenokio() like so:
Before installing, please download Rtools (https://cran.r-project.org/bin/windows/Rtools)

``` r
if(!require(remotes)) install.packages("remotes")
if (!requireNamespace("phenokio", quietly = TRUE)) {
  remotes::install_github("agronomy4future/phenokio", force= TRUE)
}
library(remotes)
library(phenokio)
```
## Example
The following example demonstrates the basic configuration of the phenokio function for analyzing soybean grains.
``` r
phenokio(
  input_folder="C:/Users/Coding",  
  output_folder="C:/Users/Coding/output", 
    image_real_cm=c(20, 20), 
    lower_hsv=c(10, 80, 70),  
    upper_hsv=c(35, 255, 255), 
    extra_hsv=list(
      list(lower=c(0, 0, 10), upper=c(180, 255, 100))),
    margin=10,   
    k_open=c(3, 3),  
    open_iter= 1,    
    k_close=c(5, 5),
    close_iter=2,   
    fill_holes=TRUE, 
    min_component_area_px=50,  
    object_min_area_cm2=0.03, 
    rel_min_frac_of_largest=0.1,  
    max_keep=500,                
    outline_color="green",       
    outline_thickness=5        
)
```
The following example demonstrates the basic configuration of the phenokio function for analyzing black bean grains.
``` r
phenokio(
  input_folder="C:/Users/Coding",  
  output_folder="C:/Users/Coding/output", 
    image_real_cm=c(20, 20), 
    lower_hsv=c(0, 0, 0), 
    upper_hsv=c(180, 255, 60),
    extra_hsv=list(
      list(lower=c(0, 0, 10), upper=c(180, 255, 100))),
    margin=20,          
    k_open=c(3, 3),     
    open_iter=1,
    k_close=c(11, 11),      
    close_iter=2,
    fill_holes=TRUE,     
    min_component_area_px= 300, 
    object_min_area_cm2=0.1,
    rel_min_frac_of_largest=0.1, 
    max_keep=500, 
    outline_color="orange",  
    outline_thickness=5
  )
```
The following example demonstrates the basic configuration of the phenokio function for analyzing wheat grains.
``` r
phenokio(
  input_folder="C:/Users/Coding",  
  output_folder="C:/Users/Coding/output", 
    image_real_cm=c(20, 20), 
    lower_hsv=c(0, 80, 50), 
    upper_hsv=c(179, 255, 255),
    extra_hsv=list(
      list(lower=c(0, 0, 10), upper=c(180, 255, 100))),
    margin=10,          
    k_open=c(3, 3),     
    open_iter=1,
    k_close=c(5, 5),      
    close_iter=1,
    fill_holes=TRUE,     
    min_component_area_px=100, 
    object_min_area_cm2=0.01,
    rel_min_frac_of_largest=0.1, 
    max_keep=500, 
    outline_color="purple",  
    outline_thickness=4
)
```

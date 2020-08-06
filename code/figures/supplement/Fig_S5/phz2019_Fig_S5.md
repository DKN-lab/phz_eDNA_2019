---
title: "Figure S5"
subtitle: 'Extracellular DNA promotes efficient extracellular electron transfer by pyocyanin in *Pseudomonas aeruginosa* biofilms.'
author: 'Scott H. Saunders, Edmund C.M. Tse, Matthew D. Yates, Fernanda Jiménez Otero, Scott A. Trammell, Eric D.A. Stemp, Jacqueline K. Barton, Leonard M. Tender and Dianne K. Newman'
output:
  html_document:
    theme: cosmo
    highlight: tango
    code_folding: show
    toc: yes
    keep_md: true
---

--------

# Notes

Fig S5A - B is a set of images.

----

Setup packages and plotting for the notebook:


```r
# Check packages
source("../../../tools/package_setup.R")

# Load packages
library(tidyverse)
library(cowplot)
library(kableExtra)
library(broom)
library(modelr)

# Code display options
knitr::opts_chunk$set(tidy.opts=list(width.cutoff=60),tidy=FALSE, echo = TRUE, message=FALSE, warning=FALSE, fig.align="center", fig.retina = 2)

# Load plotting tools
source("../../../tools/plotting_tools.R")


#Modify the plot theme
theme_set(theme_notebook())
```

# Fig. S5C 

LC-MS data was manually input from peak integrations in the Empower software. 


```r
ida_phz <- tibble(
  Day = c(1,2,3),
  PYO = c(79.86,115.4,88.6),
  PCA = c(10.4, 18.1,29.3),
  PCN = c(0.254, 0.51, 0.61)
) %>% 
  gather(key = phenazine, value = amount, -Day) %>% 
  group_by(phenazine) %>% 
  mutate(mean = ifelse(Day==1,mean(amount),NA))

plot_ida_phz <- ggplot(ida_phz, aes(x = phenazine, y = amount, shape = factor(Day))) + 
  geom_col(aes(y = mean), fill = 'light gray') + 
  geom_jitter(width = 0.1, height = 0, size = 1) + scale_shape_manual(values = c(21,22,23), guide = F)

plot_ida_phz_styled <- plot_ida_phz+
  labs(x='', y = expression(Concentration ~(mu*M)))

plot_ida_phz_styled
```

<img src="phz2019_Fig_S5_files/figure-html/unnamed-chunk-1-1.png" width="672" style="display: block; margin: auto;" />

# Create figure

Same plot is repeated to preserve legacy formating. Previously Panels 5B-C were panels S5D-E.


```r
theme_figure <- function () {
  theme_classic( ) %+replace%
    theme(
      axis.line = element_line(color = 'black', size = 0.25),
      axis.ticks = element_line(color = 'black', size =0.25),
      axis.text = element_text(color = 'black', size=8),
      axis.title=element_text(color = 'black', size=8),
      strip.text = element_text(color = 'black', size = 8),
      strip.background = element_blank(),
      legend.background = element_blank(),
      legend.title=element_text(color = 'black',size=8),
      legend.text=element_text(color = 'black',size=8),
      legend.text.align=0,
      panel.spacing = unit(0,'cm'),
      plot.margin = margin(t=0.25, b = 0.25, l = 0.25, r = 0.25, unit = 'cm'),
      plot.title = element_text(hjust = 0.5, color = 'black', size = 8)
    )
}

theme_set(theme_figure())

fig_s5 <- plot_grid(plot_ida_phz_styled, plot_ida_phz_styled, plot_ida_phz_styled ,plot_ida_phz_styled , ncol = 3, align = 'hv', axis = 'tblr', labels = c('C','','',''), label_size = 12, scale = 0.95, rel_heights = c(1,1.25))

fig_s5
```

<img src="phz2019_Fig_S5_files/figure-html/unnamed-chunk-2-1.png" width="672" style="display: block; margin: auto;" />

```r
save_plot("../../../../figures/supplement/phz2019_Fig_S5.pdf", fig_s5, base_width = 7, base_height = 4.5)
```

-----


```r
sessionInfo()
```

```
## R version 3.5.3 (2019-03-11)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS  10.15.6
## 
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] lubridate_1.7.4   hms_0.5.3         modelr_0.1.5     
##  [4] broom_0.5.2       kableExtra_1.1.0  cowplot_0.9.4    
##  [7] viridis_0.5.1     viridisLite_0.3.0 knitr_1.23       
## [10] forcats_0.4.0     stringr_1.4.0     dplyr_0.8.3      
## [13] purrr_0.3.3       readr_1.3.1       tidyr_1.0.0      
## [16] tibble_2.1.3      ggplot2_3.3.0     tidyverse_1.3.0  
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_0.2.5 xfun_0.7         haven_2.2.0      lattice_0.20-38 
##  [5] colorspace_1.4-1 vctrs_0.3.1      generics_0.0.2   htmltools_0.4.0 
##  [9] yaml_2.2.0       rlang_0.4.6      pillar_1.4.2     glue_1.3.1      
## [13] withr_2.1.2      DBI_1.0.0        dbplyr_1.4.2     readxl_1.3.1    
## [17] lifecycle_0.1.0  munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0
## [21] rvest_0.3.5      evaluate_0.14    labeling_0.3     Rcpp_1.0.2      
## [25] scales_1.0.0     backports_1.1.4  webshot_0.5.1    jsonlite_1.6    
## [29] fs_1.3.1         gridExtra_2.3    digest_0.6.21    stringi_1.4.3   
## [33] grid_3.5.3       cli_1.1.0        tools_3.5.3      magrittr_1.5    
## [37] crayon_1.3.4     pkgconfig_2.0.3  ellipsis_0.3.0   xml2_1.2.2      
## [41] reprex_0.3.0     assertthat_0.2.1 rmarkdown_1.13   httr_1.4.1      
## [45] rstudioapi_0.10  R6_2.4.0         nlme_3.1-137     compiler_3.5.3
```

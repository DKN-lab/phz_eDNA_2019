---
title: "Figure 6"
subtitle: 'Extracellular DNA promotes efficient extracellular electron transfer by pyocyanin in *Pseudomonas aeruginosa* biofilms.'
author: 'Scott H. Saunders, Edmund C.M. Tse, Matthew D. Yates, Fernanda Jiménez-Otero, Jacqueline K. Barton, Leonard M. Tender and Dianne K. Newman'
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

Panel A is a diagram.

To see how we got from the raw electrochemical scans to the datasets used here, please see the following processing notebooks:

* [IDA ∆phz biofilm processing](https://scott-saunders.github.io/phz_eDNA_2019/code/processing/IDA_dPHZ/IDA_dPHZ_processing.html)
* [IDA WT biofilm processing](https://scott-saunders.github.io/phz_eDNA_2019/code/processing/IDA_WT/IDA_WT_processing.html)
* [IDA blank processing](https://scott-saunders.github.io/phz_eDNA_2019/code/processing/IDA_blank/IDA_blank_processing.html)

Then see how those data were analyzed in these notebooks for supplemental figures S6 and S7:

* [Fig. S6](https://scott-saunders.github.io/phz_eDNA_2019/code/figures/supplement/Fig_S6/phz2019_Fig_S6.html)
* [Fig. S7](https://scott-saunders.github.io/phz_eDNA_2019/code/figures/supplement/Fig_S7/phz2019_Fig_S7.html)

These supplemental figure notebooks produced model coefficients that are used in this notebook.

----

Setup packages and plotting for the notebook:


```r
# Load packages
library(tidyverse)
library(cowplot)
library(kableExtra)

# Code display options
knitr::opts_chunk$set(tidy.opts=list(width.cutoff=60),tidy=FALSE, echo = TRUE, message=FALSE, warning=FALSE, fig.align="center", fig.retina = 2)

# Load plotting tools
source("../../tools/plotting_tools.R")


#Modify the plot theme

theme_set(theme_notebook())
```

# Fig. 6B

need to read in the model coefficients and grid predictions for these two datasets.


```r
df_dphz_nls_1 <- read_csv("../supplement/Fig_S6/phz2019_dPHZ_Dphys_nls_coefs.csv") %>% 
  filter(exp == 1 & run == 1 & term == 'a')

df_blank_nls_75 <- read_csv("../supplement/Fig_S6/phz2019_blank_Dphys_nls_coefs.csv") %>% 
  filter(PHZadded == '75uM' & term == 'a')

df_dphz_swv_decay <- read_csv("../../processing/processed_data/phz_eDNA_2019_signals_long.csv") %>% 
  filter(echem == 'SWV' & exp == 1 & run == 1 & reactor == 'transfer' & electrode == 'i1') %>% 
  mutate(IDA = 'biofilm', signal = signal - df_dphz_nls_1$estimate) %>% 
  select(time, signal, IDA)

df_blank_decay <- read_csv("../../processing/processed_data/phz_eDNA_2019_swv_blank_tran_time_signals.csv") %>% 
  filter(PHZadded == '75uM') %>% 
  mutate(IDA = 'blank', signal = signal - df_blank_nls_75$estimate) %>%
  select(time, signal, IDA)

df_dphz_preds <- read_csv("../supplement/Fig_S6/phz2019_dPHZ_Dphys_preds.csv") %>% 
  filter(exp == 1 & run == 1) %>% 
  select(time, pred, pred_low, pred_high) %>% 
  mutate(pred = pred - df_dphz_nls_1$estimate,
         pred_low = pred_low - df_dphz_nls_1$estimate, 
         pred_high = pred_high - df_dphz_nls_1$estimate)%>% 
  mutate(IDA = 'biofilm')

df_blank_preds <- read_csv("../supplement/Fig_S6/phz2019_blank_Dphys_preds.csv") %>% 
  filter(PHZadded == '75uM') %>% 
  select(time, pred, pred_low, pred_high) %>% 
  mutate(pred = pred - df_blank_nls_75$estimate,
         pred_low = pred_low - df_blank_nls_75$estimate, 
         pred_high = pred_high - df_blank_nls_75$estimate) %>% 
  mutate(IDA = 'blank')

df_preds <- bind_rows(df_dphz_preds, df_blank_preds)

df_decays <- bind_rows(df_dphz_swv_decay, df_blank_decay)

plot_blank_dphz_decay <- ggplot(df_preds, aes(x = time, y = pred, group = IDA, fill = time)) + geom_ribbon(aes(ymin = pred_low, ymax = pred_high), fill = 'light gray') +
  geom_path(linetype = 2)+
  geom_point(data =df_decays, aes(x = time, y = signal) , shape = 21) + guides(fill = 'none')

plot_blank_dphz_decay_styled <- plot_blank_dphz_decay+
  labs(x = 'Time (min)', y = expression(I[swv]~(nA)), fill = 'Time (min)') +
  scale_fill_viridis(guide = F) +
  scale_y_continuous(labels = nA_label)


plot_blank_dphz_decay_styled
```

<img src="phz2019_Fig_6_files/figure-html/unnamed-chunk-1-1.png" width="672" style="display: block; margin: auto;" />

Let's subtract the a intercept value of these models so that we can compare the decay.



# Fig. 6C


```r
df_swv <- read_csv("../../processing/processed_data/phz_eDNA_2019_swv_raw_1_1.csv")
df_swv_sig <- read_csv("../../processing/processed_data/phz_eDNA_2019_signals_long.csv") %>% 
  filter(echem == 'SWV' & exp == 1 & run == 1 & electrode == 'i1' & reactor == 'transfer')

# Plot Layout
plot_swv <- ggplot(df_swv %>% filter(E<=0), aes(x = E , y = current )) +
  geom_vline(xintercept = -0.265, linetype=2, color = 'gray') + 
  geom_path(aes(group = rep, color = rep)) + 
  geom_point(data = df_swv_sig, aes(x = E_from_maxs , y = current_from_maxs , fill = rep), shape = 21,size = 2) 

# Plot Styling
plot_swv_styled <- plot_swv + 
  annotate('text', x = -0.2, y = 3.5e-6,label = expression({E^0}[pyo]) )+
  scale_x_reverse(labels = mV_label)+
  scale_y_continuous(labels = nA_label)+
  scale_fill_viridis(guide = F) + 
  scale_color_viridis(guide = F) + 
  labs(x = "E (mV vs. Ag/AgCl)", y = expression(I[swv]~(n*A)), color = "Scan #") + 
  theme(legend.position = c(0.15,0.75), legend.background = element_blank())

plot_swv_styled
```

<img src="phz2019_Fig_6_files/figure-html/unnamed-chunk-2-1.png" width="672" style="display: block; margin: auto;" />

## Inset

convert to norm time??


```r
# Plot Layout
plot_swv_sig <- ggplot(data = df_swv_sig, aes(x = time, y = signal))+
  geom_point(shape = 21, size = 1, aes(fill = time))

# Plot Styling
plot_swv_sig_styled <- plot_swv_sig +
  labs(x = 'Time (min)', y = expression(I[swv]~(nA)), fill = 'Time (min)') +
  scale_fill_viridis(guide = F) +
  scale_y_continuous(labels = nA_label)

plot_swv_sig_styled
```

<img src="phz2019_Fig_6_files/figure-html/unnamed-chunk-3-1.png" width="672" style="display: block; margin: auto;" />

# Fig. 6D

need to mod

```r
df_gc <- read_csv("../../processing/processed_data/phz_eDNA_2019_gc_raw_1_1.csv")
df_gc_sig <- read_csv("../../processing/processed_data/phz_eDNA_2019_signals_long.csv") %>% 
  filter(echem == 'GC' & exp == 1 & run == 1 & reactor == 'transfer')


# Plot Layout
plot_gc <- ggplot() + 
  geom_vline(xintercept = -0.265, linetype=2, color = 'gray') + 
  geom_path(data=df_gc %>% filter(electrode=='i1'), 
            aes(x = E , y = current , group = rep, color = rep)) + 
  geom_path(data=df_gc %>% filter(electrode=='i2'), 
            aes(x = E , y = current, group = rep, color = rep)) +
  geom_point(data = df_gc_sig, 
             aes(x = E_from_maxs , y = -current_from_maxs , fill = rep), shape = 21, size = 2) 

# Plot styling
plot_gc_styled <- plot_gc +
  annotate('text', x = -0.2, y = 3.5e-7,label = expression({E^0}[pyo]) )+
  scale_x_reverse(labels = mV_label)+
  scale_y_continuous(labels = nA_label)+
  scale_fill_viridis(guide=F) + 
  scale_color_viridis(guide=F) + 
  labs(x = "E (mV vs. Ag/AgCl)", y = expression(I[gc]~(nA)), color = "Scan #") + 
  theme(legend.position = c(0.15,0.75), 
        legend.background = element_blank())

plot_gc_styled
```

<img src="phz2019_Fig_6_files/figure-html/unnamed-chunk-4-1.png" width="672" style="display: block; margin: auto;" />


## Inset

need to mod

```r
# Plot Layout
plot_gc_sig <- ggplot(data = df_gc_sig, aes(x = time, y = signal))+
  geom_point(shape = 21, size = 1, aes(fill = time))

# Plot Styling
plot_gc_sig_styled <- plot_gc_sig +
  labs(x = 'Time (min)', y = expression(I[gc]~(nA)), fill = 'Time (min)') +
  scale_fill_viridis(guide = F) +
  scale_y_continuous(labels = nA_label)

plot_gc_sig_styled
```

<img src="phz2019_Fig_6_files/figure-html/unnamed-chunk-5-1.png" width="672" style="display: block; margin: auto;" />

# Fig. 6E

need to mod

```r
df_swv_gc_sig <- read_csv("../../processing/processed_data/phz_eDNA_2019_swv_gc_signals.csv") %>% filter(exp == 1 & run == 1)

# Plot Layout
plot_swv_gc <- ggplot(df_swv_gc_sig %>% filter(rep>0), 
                           aes(x = signal_SWV, y = signal_GC)) +
  geom_smooth(method = 'lm', se = T, color = 'black', linetype = 2) +
  geom_point(shape = 21, size = 2, aes(fill = time_SWV))

# Plot Styling
plot_swv_gc_styled <- plot_swv_gc + 
  scale_x_continuous(breaks=c(0,2.5e-7,5.0e-7, 7.5e-7, 1.0e-6), labels = nA_label)+
  scale_y_continuous(labels = nA_label)+
  labs(x = expression(I[swv]~(nA)), y = expression(I[gc]~(nA))) +
  scale_fill_viridis(guide=F)


plot_swv_gc_styled
```

<img src="phz2019_Fig_6_files/figure-html/unnamed-chunk-6-1.png" width="672" style="display: block; margin: auto;" />

# Fig. 6F

1. Take biofilm linear models from fig. S6 and convert into Dap.
2. Take blank linear model from Fig. S7 and convert into Dap
3. Take nonlinear models from fig. S6 and join with Dap data. Then calculate Dphys from nonlinear fit coefficient and Dap.

## $D_{ap}$


```r
df_dphz_lm <- read_csv("../supplement/Fig_S6/phz2019_dPHZ_Dap_lm_coefs.csv") %>% mutate(IDA = 'biofilm')

df_blank_lm <- read_csv("../supplement/Fig_S7/phz2019_blank_Dap_lm_coefs.csv") %>% mutate(IDA = 'blank')

df_all_lm <- bind_rows(df_dphz_lm, df_blank_lm)
```

Let's redefine the function we used in [Fig. S7]() to calculate $D_{ap}$ from the slope of the SWV vs. GC line.


```r
dap_from_swvGC <- function(m, t_p=1/(2*300)){
  
  psi <-  0.7
  A <-  0.025 #cm^2
  S <-  18.4 #cm
  
  d_ap <- (m*A*psi)^2 / (S^2 * pi * t_p)
  
  d_ap
}
```

Now we'll convert slope into $D_{ap}$ estimates (with confidence intervals).


```r
df_all_dap <- df_all_lm %>% 
  filter(term=='signal_SWV') %>% 
  mutate(dap=dap_from_swvGC(m = estimate)) %>% 
  mutate(dap_high = dap_from_swvGC(m = conf.high)) %>% 
  mutate(dap_low = dap_from_swvGC(m = conf.low)) %>% 
# Add back id variables for matching below  
  mutate(run = case_when(
    run_id == 'Rep 1' ~ 1,
    run_id == 'Rep 2' ~ 2,
    run_id == 'Rep 3' ~ 3
  )) %>% 
    mutate(exp = case_when(
    exp_id == 'Biofilm 1' ~ 1,
    exp_id == 'Biofilm 2' ~ 2,
    exp_id == 'Blank 1' ~ 3
  )) %>% 
  mutate(run = ifelse(exp_id == 'Blank 1',4,run))


ggplot(df_all_dap, aes(x = exp_id, y = dap)) + geom_pointrange(aes(ymin = dap_low, ymax = dap_high))
```

<img src="phz2019_Fig_6_files/figure-html/unnamed-chunk-9-1.png" width="672" style="display: block; margin: auto;" />

## $D_{phys}$

Combine nls fits for blank and dphz. Add Dap estimates for each. Add I0 estimates for each. Define Dphys calculator function then apply. 


```r
df_dphz_nls <- read_csv("../supplement/Fig_S6/phz2019_dPHZ_Dphys_nls_coefs.csv")
df_dphz_i0 <- read_csv("../../processing/processed_data/phz_eDNA_2019_swv_signals.csv") %>% 
  filter(reactor == 'soak' & electrode == 'i1') %>% 
  mutate(run = as.double(run)) %>% 
  select(exp, run, i0 = signal)

df_dphz_nls_i0 <- left_join(df_dphz_nls, df_dphz_i0, by = c('exp','run')) %>% mutate(IDA = 'biofilm')



df_blank_nls <- read_csv("../supplement/Fig_S6/phz2019_blank_Dphys_nls_coefs.csv")
df_blank_i0 <- read_csv("../../processing/processed_data/phz_eDNA_2019_swv_blank_signals.csv") %>% 
  filter(reactor == 'soak') %>% 
  select(PHZadded, i0 = signal)

df_blank_nls_i0 <- left_join(df_blank_nls, df_blank_i0, by = c('PHZadded')) %>% 
  mutate(exp = 3, IDA = 'blank') %>% 
# Add an ID variable for matching below
  mutate(run = case_when(
    PHZadded == '10uM' ~ 4,
    PHZadded == '25uM' ~ 4,
    PHZadded == '50uM' ~ 4,
    PHZadded == '75uM' ~ 4,
    PHZadded == '100uM' ~ 4
  ))

df_all_nls_i0 <- bind_rows(df_dphz_nls_i0, df_blank_nls_i0) %>% filter(term == 'b') %>% 
  select(exp, run,b_estimate = estimate, b_low = conf.low, b_high = conf.high, i0, IDA, PHZadded)

df_all_dap_i0 <- left_join(df_all_nls_i0, df_all_dap %>% select(exp, run, IDA, dap, dap_high, dap_low),
                           by = c('exp','run','IDA'))
```


```r
dphys_from_nls <- function(estimate, i0, dap, t_s=0.1){
  
  dphys <- ( i0^2 * dap * t_s ) / (pi* estimate^2)
  
  dphys
}

df_dphys <- df_all_dap_i0 %>% 
  mutate(dphys = dphys_from_nls(estimate = b_estimate, i0 = i0, dap = dap, t_s = 0.02)) %>% 
  mutate(dphys_high = dphys_from_nls(estimate = b_high, i0=i0, dap = dap_high, t_s = 0.02)) %>%
  mutate(dphys_low = dphys_from_nls(estimate = b_low, i0=i0, dap = dap_low, t_s = 0.02))

ggplot(df_dphys, aes(x = IDA, y = dphys)) + geom_pointrange(aes(ymin = dphys_low, ymax = dphys_high, shape = factor(exp))) + scale_y_log10()
```

<img src="phz2019_Fig_6_files/figure-html/unnamed-chunk-11-1.png" width="672" style="display: block; margin: auto;" />



```r
df_plot_dap <- df_dphys %>% 
  select(exp, run, IDA, estimate = dap, estimate_high = dap_high, estimate_low = dap_low ) %>% 
  mutate(coef = 'Dap') %>% 
  distinct()

df_plot_dphys <- df_dphys %>% 
  select(exp, run, IDA, estimate = dphys, estimate_high = dphys_high, estimate_low = dphys_low ) %>% 
  mutate(coef = 'Dphys')

df_plot_dap_dphys <- bind_rows(df_plot_dap, df_plot_dphys)

plot_dap_dphys <- ggplot(df_plot_dap_dphys, aes(x = coef, y = estimate, shape = factor(exp))) + 
  geom_pointrange(aes(ymin = estimate_low, ymax = estimate_high), position = position_jitter(width =0.1, height = 0), fatten = 3,stroke = 0.5) + facet_wrap(~IDA, scales = 'free') + 
  scale_y_log10(limits = c(1e-8, 2e-5), labels = scales::trans_format("log10", scales::math_format(10^.x))) 

plot_dap_dphys_styled <- plot_dap_dphys +
  labs(x = NULL, y = expression(D~(cm^2 / sec)))+ 
  scale_shape_manual(values = c(21,22,23), guide = F)

plot_dap_dphys_styled
```

<img src="phz2019_Fig_6_files/figure-html/unnamed-chunk-12-1.png" width="672" style="display: block; margin: auto;" />

# Create figure


```r
theme_set(theme_figure())

fig_6_insets <- plot_grid(NULL, NULL,plot_swv_sig_styled, plot_gc_sig_styled, ncol = 2, rel_heights = c(1,2))

fig_6 <- plot_grid(fig_6_insets, plot_blank_dphz_decay_styled, 
                   plot_swv_styled, plot_gc_styled, 
                   plot_swv_gc_styled, plot_dap_dphys_styled, 
                   ncol = 2, labels = 'AUTO', label_size = 12, scale = 0.95, align = 'hv', axis = 'tblr')

fig_6 <- plot_grid(fig_6_insets, plot_swv_styled, plot_gc_styled, 
                   plot_swv_gc_styled, plot_blank_dphz_decay_styled, plot_dap_dphys_styled, 
                   ncol = 3, labels = 'AUTO', label_size = 12, scale = 0.95, align = 'hv', axis = 'tblr')

#fig_6

save_plot("../../../figures/phz2019_Fig_6.pdf", fig_6, base_height = 4, base_width = 7)
```

-----


```r
sessionInfo()
```

```
## R version 3.5.2 (2018-12-20)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS Mojave 10.14.6
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
##  [1] viridis_0.5.1     viridisLite_0.3.0 kableExtra_1.0.1 
##  [4] cowplot_0.9.4     forcats_0.3.0     stringr_1.3.1    
##  [7] dplyr_0.8.1       purrr_0.2.5       readr_1.3.1      
## [10] tidyr_0.8.2       tibble_2.1.3      ggplot2_3.2.0    
## [13] tidyverse_1.2.1  
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_0.2.5 xfun_0.7         haven_2.0.0      lattice_0.20-38 
##  [5] colorspace_1.4-0 generics_0.0.2   htmltools_0.3.6  yaml_2.2.0      
##  [9] rlang_0.4.0      pillar_1.3.1     glue_1.3.1       withr_2.1.2     
## [13] modelr_0.1.2     readxl_1.2.0     munsell_0.5.0    gtable_0.2.0    
## [17] cellranger_1.1.0 rvest_0.3.2      evaluate_0.14    labeling_0.3    
## [21] knitr_1.23       broom_0.5.1      Rcpp_1.0.1       scales_1.0.0    
## [25] backports_1.1.3  webshot_0.5.1    jsonlite_1.6     gridExtra_2.3   
## [29] hms_0.4.2        digest_0.6.18    stringi_1.2.4    grid_3.5.2      
## [33] cli_1.1.0        tools_3.5.2      magrittr_1.5     lazyeval_0.2.1  
## [37] crayon_1.3.4     pkgconfig_2.0.2  xml2_1.2.0       lubridate_1.7.4 
## [41] assertthat_0.2.1 rmarkdown_1.13   httr_1.4.0       rstudioapi_0.9.0
## [45] R6_2.4.0         nlme_3.1-140     compiler_3.5.2
```

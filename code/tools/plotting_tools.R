library(ggplot2)
library(viridis)

# Theme

theme_notebook <- function () { 
  theme_classic( ) %+replace% 
    theme(
      axis.line = element_line(color = 'black'),
      axis.ticks = element_line(color = 'black'),
      axis.text = element_text(color = 'black',),
      axis.title=element_text(color = 'black'),
      strip.text = element_text(color = 'black'),
      strip.background = element_blank(),
      legend.background = element_blank(),
      legend.title=element_text(color = 'black'),
      legend.text=element_text(color = 'black'),
      legend.text.align=0, 
      plot.title = element_text(hjust = 0.5)
    )
}

theme_figure <- function () {
  theme_classic( ) %+replace%
    theme(
      axis.line = element_line(color = 'black', size = 0.25),
      axis.ticks = element_line(color = 'black', size =0.25),
      axis.text = element_text(color = 'black', size=6),
      axis.title=element_text(color = 'black',size=6),
      strip.text = element_text(color = 'black',size = 6),
      strip.background = element_blank(),
      legend.background = element_blank(),
      legend.title=element_text(color = 'black',size=6),
      legend.text=element_text(color = 'black',size=6),
      legend.text.align=0,
      panel.spacing = unit(0,'cm'),
      plot.margin = margin(t=0.25, b = 0.25, l = 0.25, r = 0.25, unit = 'cm'),
      plot.title = element_text(hjust = 0.5)
    )
}

# Color palette

colorblind_palette <- c("#E69F00", "#56B4E9", "#000000","#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

scale_colour_discrete <- function(...) scale_colour_manual(values = colorblind_palette)
scale_fill_discrete <- function(...) scale_fill_manual(values = colorblind_palette)

scale_colour_continuous <- scale_colour_viridis
scale_fill_continuous <- scale_fill_viridis

# Custom axis labels

nA_label <- function(x){
# from A to nA
  lab <- x / 1e-9
}

ns_label <- function(x){
  # from s to ns
  lab <-  paste(x / 1e-9, "ns", sep = '')
}

mV_label <- function(x){
# From V to mV
  lab <- x * 1000
}

fold_label <- function(x){
#e.g. 2 to 2x
  lab <- paste(x, "x", sep = '')
}
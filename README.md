---
title: "Extracellular DNA promotes efficient extracellular electron transfer by pyocyanin in *Pseudomonas aeruginosa* biofilms."
author: 'Scott H. Saunders, Edmund C.M. Tse, Matthew D. Yates, Fernanda Jiménez Otero, Scott A. Trammell, Eric D.A. Stemp, Jacqueline K. Barton, Leonard M. Tender and Dianne K. Newman'
fontsize: 12pt
output:
  html_document:
    theme: cosmo
    code_folding: hide
    keep_md: true
---

--------

# Abstract

Extracellular electron transfer (EET), the process whereby cells access electron acceptors or donors that reside many cell lengths away, enables metabolic activity by microorganisms, particularly under oxidant-limited conditions that occur in multicellular bacterial biofilms. Although different mechanisms underpin this process in select organisms, a widespread strategy involves extracellular electron shuttles, redox-active metabolites that are secreted and recycled by diverse bacteria. How these shuttles catalyze electron transfer within biofilms without being lost to the environment has been a long-standing question. Here, we show that phenazine electron shuttles mediate efficient EET through interactions with extracellular DNA (eDNA) in *Pseudomonas aeruginosa* biofilms, which are important in nature and disease. Retention of pyocyanin (PYO) and phenazine carboxamide in the biofilm matrix is facilitated by binding to eDNA. In vitro, different phenazines can exchange electrons in the presence or absence of DNA and phenazines can participate directly in redox reactions through DNA; the biofilm eDNA can also support rapid ET between intercalators. Electrochemical measurements of biofilms indicate that retained PYO supports an efficient redox cycle with rapid EET and slow loss from the biofilm. Together, these results establish that eDNA plays a previously unrecognized role facilitating phenazine metabolic processes in P. aeruginosa biofilms, suggesting a model for how extracellular electron shuttles achieve retention and efficient EET in biofilms. 

--------

<br>

If you would like to see the website with rendered code notebooks that generate the figures in this paper go to [this website](https://dkn-lab.github.io/phz_eDNA_2019/) and look at the section "Computational notebooks."

<br>

-----

# How to use this GitHub repository

Hello! This is a GitHub repository for the paper described above. In the window above, you can see that this repository contains a file structure with three main folders: `code`, `data`, and `figures`. To put it simply then, this repository contains the data and code to create the figures in the paper!

### Step 1. Don't panic

You may not be familiar with GitHub repositories and this may look a little scary to you. Don't panic - you can still understand this repository without understanding the underlying code or every single file. First, don't worry about the commit messages, which are the short pieces of text next to the files (with time stamps like "2 months ago"), and focus on the three folders `code`, `data`, and `figures`. 

### `/figures`

To orient you, let's start with the final product - the figures. If you open the `figures` folder you will see .pdf files for each of the main figures 1-5, as well as a `supplement` folder containing .pdf files for figures S1 - S7. Clicking on any of these files will show you a part of that figure that was generated from the data and code within this respository. All of these figures had final modfications done in Adobe Illustrator, so some of the legends and things may not look perfect, but the point is that anyone could see how the plots in each figure were generated from the underlying data. 

### `/data`

Now that you've seen the final product, go ahead and explore the underlying data. In the `data` folder you will find subfolders labelled `Electrochemistry`, `LC-MS`, and `Spectroscopy`, because those are the three general types of data that we collected for this paper. Going into each of the subfolders you will find mostly raw data saved as .csv files. Clicking on those files (if they're not too big), will show you the spreadsheet like data that underlies everything in the paper. If you feel so inclined, you could download any of these files and plot the data yourself in your favorite program (e.g. Excel). 

### `/code`

How do we get from data to figures in a well documented manner? Code. In the `code` folder you will find subfolders `figures`, `processing`, and `tools`. `figures` contains the code specifically used to make each main and supplemental figure. You can see that each figure has it's own folder and code notebook file (.Rmd). Each one of these notebooks outputs the corresponding .pdf file for the figure. If you are familiar with R code, take a look at some of the notebooks! Even if you are familiar, it is somewhat easier to look at the pre-rendered .html versions of these notebooks. **Therefore, I recommend checking out [this website](https://dkn-lab.github.io/phz_eDNA_2019/) and clicking on the links in the "Computational notebooks" section.**

The `processing` folder contains notebooks that do not specifically make figures, but help "process" complicated raw data files into forms that are later used in the figure notebooks. These are also pre-rendered on the website. 

Lastly, there's the `tools` folder. This contains a few .R scripts that have made my life easier. Specifically, there are `echem_processing_tools.R`, which are used in the processing notebooks to read in large numbers of files and extract simple metadata from the filenames and file headers. `plotting_tools.R` has the ggplot2 themes, color palettes and scale labels used throughout the notebooks. 

If you feel so inclined, you could download this whole repository and open any of the .Rmd notebook files in R, and recreate or modify any of the figures. The version numbers of R and associated packages are documented at the bottom of each notebook. 


--------

## License
![](https://licensebuttons.net/l/by/3.0/88x31.png)

All creative works (writing, figures, etc) are licensed under the [Creative
Commons CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) license. All software is distributed under the standard MIT license as follows

```
Copyright 2019 Scott H Saunders 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
`

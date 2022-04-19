# SFPSA: Self Fulfilling Prophecy Sensitivity Analysis
An R-package that allows users to implement the methods of **DO FORECASTS OF BANKRUPTCY CAUSE BANKRUPTCY?A MACHINE LEARNING SENSITIVITY ANALYSIS.** The github repo for that paper can be found [here](https://github.com/demetrios1/bankruptcy_sensitivity), with the arxiv link [here](https://arxiv.org/abs/2106.04503).  Go [here](https://demetrios1.github.io/SFPSA/) for the vignette.

How to get started.  
<!--The first thing you want to  do is 

```install.packages('https://github.com/demetrios1/bankruptcy_sensitivity/blob/main/monotone_bart/fastbart_2.0.tar.gz', repos = NULL, type="source")```

This will install one of the necessary packages from [Jared Murray](https://jaredsmurray.github.io/), called [monbart](https://github.com/jaredsmurray/monbart).  After this is done, 
-->
```
library(devtools)
install_github("jaredsmurray/monbart", ref='main')
install_github("demetrios1/SFPSA", ref="main") 
```
This will install the package and the requisite packages attached. See the attached [github pages](https://demetrios1.github.io/SFPSA/) for a tutorial!

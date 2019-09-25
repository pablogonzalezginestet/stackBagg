# EnsBagg: Ensemble IPCW Bagging

**EnsBagg** is an `R` package that uses inverse probability of censoring weighted
(IPCW) bagging approach as a pre-processing step to allow for all existing and any
newly developed machine learning methods for classification to be applied to right-censored data with or
without competing risk.

**EnsBagg** provides the files to reproduce the simulation studies and the real data application in the forthcoming paper:
* Gonzalez Ginestet, P. et al. (2019). "Ensemble IPCW bagging: a case study in the HIV care
registry".


## Important Note

The real data application in the paper relies on Swedish HIV care register (InfCareHIV) that cannot be openly shared. The data available on the folder [missing] is based on a fake, simulated  data that mimics the data used for the analysis. The R scripts in this repository are the same as those used for the analysis, but analysis results cannot be exactly reproduced. 

## Installation

You can install the development version of `EnsBagg` from [GitHub](https://github.com/pablogonzalezginestet/EnsBagg) with:


```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github("pablogonzalezginestet/EnsBagg")
```

## Example

```R
library(EnsBagg)
train<- as.data.frame(EnsBagg::exampleData[[1]])
test<- as.data.frame(EnsBagg::exampleData[[2]])
pred_auc <- EnsBagg::ensBagg(train.data = train, test.data = test, xnam=c("X9","X18"),tao = 26.5,weighting = "CoxPH",folds = 5)
```

## Vignette


See the vignette for details about examples usage of the functions found in  `EnsBagg` .


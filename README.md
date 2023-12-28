## Install

``` r
install.packages("devtools")
devtools::install_github("Biores410/FateND")
library(FateND)
``` 

## Example
``` r
# Example 1
library(RColorBrewer)

data(committed_data)

fit_committed <- get_fit_committed(committed_data)

genepair <- get_genepair(fit_committed[[1]],fit_committed[[2]])

TF_committed <- FateND_committed(fit_committed[[3]],genepair[[1]],genepair[[2]])
```

``` r
#Example 2
library(RColorBrewer)

data(embryos_data)

fit_embryos <- get_fit_embryos(epi_data,hypo_data)

genepair <- get_genepair(fit_embryos[[1]],fit_embryos[[2]])

TF_embryos <- FateND_embryos(fit_embryos[[3]],fit_embryos[[4]],genepair[[1]],genepair[[2]])
```

``` r
# Example 3
library(RColorBrewer)

data(HSC_data)

fit_HSC <- get_fit_HSC(Ery_data,DC_data,Mono_data)

genepair <- get_genepair(fit_HSC[[1]],fit_HSC[[2]])

TF_HSC <- FateND_HSC(fit_HSC[[3]],fit_HSC[[4]],fit_HSC[[5]],genepair[[1]],genepair[[2]])
```

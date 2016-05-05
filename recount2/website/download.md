
```{r}
## If needed
install.packages('devtools')

## Install recount from GitHub
devtools::install_github('leekgroup/recount')
## Will be available from Bioconductor in the future
# source('http://bioconductor.org/biocLite.R')
# biocLite('recount')

## Browse the vignetets for a quick description of how to use the package
library('recount')
browseVignettes('recount')

## Run the example for downloading a study
example('download_study', package = 'recount')
```

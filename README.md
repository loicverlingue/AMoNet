# AMoNet package V1
Author: Loic Verlingue

## Installation

You can install AMoNet R package directly from Rstudio with the following steps.

Be sure to have a Curie intranet connection to acess gitlab.curie.fr.

### From gitlab

Within the ```install_git()``` function you can:

* Change the ```user``` and ```password``` with your own.

* Add the argument ```build_opts=c("--no-resave-data", "--no-manual")``` to ```install_git()``` and build vignettes before installation, but should take additional installation time.

```{r install}
library(devtools)
install_git('https://user:password@gitlab.curie.fr/lverling/amonet.git')
library(AMoNet)
```

## Usage

A step by step introduction to basic AMoNet usage are provided in the vignettes.

To show vignettes, you have several options:

* Build vignettes yourself with
```{r}
build_vignettes("AMoNet")
browseVignettes("AMoNet")
```

* Dowload vignettes in Rmarkdown format at https://gitlab.curie.fr/lverling/amonet/tree/master/vignettes and open it with R studio.

* If installation was made with building vignettes (or installation with source package) just do:
```{r}
browseVignettes("AMoNet")
```

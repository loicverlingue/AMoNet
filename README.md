# AMoNet package V1
Author: Loic Verlingue
Date: 2019-10-30

# Installation

Welcome to AMoNet R package.

You can install AMoNet R package directly from Rstudio with the following steps.

To download it from github or gitlab.curie.fr (be sure to have a Curie intranet connection), you may need to install the packages with dependencies in AMoNet first: check warnings.

Otherwise, you can install AMoNet package from source file.

## From gitlab or github

Within the ```install_git()``` function you can:

*with gitlab
Change the ```user``` and ```password``` with your own.

```{r, install, eval=FALSE}
library(devtools)
install_git("https://user:password@gitlab.curie.fr/lverling/amonet.git", dependencies = TRUE)
library(AMoNet)
```

*with github
```{r, install, eval=FALSE}
library(devtools)
install_git("https://github.com/loicverlingue/AMoNet.git", dependencies = TRUE)
library(AMoNet)
```

## From source file

```{r, install from source, eval=FALSE}
# Download source package: AMoNet_0.1.0.tar.gz 
install.packages("path_to_file/AMoNet_0.1.0.tar.gz", repos = NULL, type="source")
```

## From CRAN
under testing

# Usage

A step by step introduction to basic AMoNet usage are provided in the vignettes. To show vignettes:

* If you have installed AMoNet from source just run
```{r, eval=FALSE}
browseVignettes("AMoNet")
```

* Or build vignettes yourself with
```{r, eval=FALSE}
build_vignettes("AMoNet")
browseVignettes("AMoNet")
```

* Dowload vignettes in Rmarkdown format at https://gitlab.curie.fr/lverling/amonet/tree/master/vignettes or https://github.com/loicverlingue/AMoNet/tree/master/vignettes and open it with R studio.


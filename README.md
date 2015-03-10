# ***addmtoolbox***
---

The ***addmtoolbox*** package provides a set of functions that facilitate the process of fitting (attentional) drift diffusion models to experimental data.

The basic drift diffusion process is implemented in a highly optimized ***c++*** function, which is complemented with a host of ***R*** functions that help with *data preparation*, facilitation of *parameter optimization* and *parallelization*.

The package is in an early development stage, but basic functionality is in working condition. Make sure to read the **Documentation** section for info how to use the package. Content that is mentioned in the **package vignettes** is currently the safest to run. 

Having that said, the package will improve dramatically in the next weeks/month so don't be put off by a few bugs in the meanwhile.


## Installation
---
The following guide is specifically for Mac OS systems. Windows systems have not yet been tested. On Linux/Ubuntu system, directly follow up with step: *"Installing addmtoolbox from github"*,  as gsl is often/always included in the distribution as well as all utilities that you gain on Macs by installing [***Xcode***](https://developer.apple.com/xcode/). This guide will carry you through the complete installation process, starting from the installation of the ***R progamming language*** on your system.

### Installation relevant dependencies
The ***addmtoolbox*** package depends notably on the ***RcppZiggurat*** package which in turn depends on the ***RcppGSL*** package. For the ***RcppGSL*** package to work, we need to install the [***GNU GSL scientific library***](http://www.gnu.org/software/gsl/). 

First, you need to have [***Xcode***](https://developer.apple.com/xcode/) installed on your system, which if not present can be downloaded for free from the [***Apple Appstore***](https://itunes.apple.com/de/app/xcode/id497799835?mt=12). After installing 
[***Xcode***](https://developer.apple.com/xcode/), you also need the [***command line tools***](https://idmsa.apple.com/IDMSWebAuth/login?&appIdKey=891bd3417a7776362562d2197f89480a8547b108fd934911bcbea0110d07f757&path=%2F%2Fdownloads%2Findex.action) from Apple (Apple wants you to sign in to download those).

### Installing GSL Library
Below is a step by step guide to installing the [***GNU GSL scientific library***](http://www.gnu.org/software/gsl/).

1. First [download](http://mirrors.ibiblio.org/gnu/ftp/gnu/gsl/) the current version and unpack it if not done automatically. (choose the file which has the highest version number: ***gsl-1.16.tar.gz*** at the time of writing)
2. Go to the terminal and navigate to the directory to which you downloaded the ***gsl library***. (useful terminal commands: cd, ls)
3. Type ```./configure``` wait until the the execution finished and then type ```make``` to build the library. (you can find these instructions in the ***INSTALL*** file inside of the downloaded ***gsl*** folder). The build process takes a while now.
4. Lastly type ```make install```, which will finally install the library on your system. It is likely that you have to use ```sudo make install```, as the default folder which the library is installed in requires superuser access.

### R
Please proceed only if you have successfully installed [***Xcode***](https://developer.apple.com/xcode/) and the [***command line tools***](https://idmsa.apple.com/IDMSWebAuth/login?&appIdKey=891bd3417a7776362562d2197f89480a8547b108fd934911bcbea0110d07f757&path=%2F%2Fdownloads%2Findex.action). If you install ***R*** without these things present on you system, there is some risk that the ***addmtoolbox*** will not run/install. 

Given that you followed the instructions so far, we are now ready to  get the ***R*** programming language ready for usage on our system (if not present). The package is built under ***R Version 3.1.2.*** and I strongly suggest to [install this](http://cran.cnr.berkeley.edu) or (if present) a later version.

### RStudio
The de-facto standard ***IDE*** for ***R*** programming. I strongly suggest that you download [this version](http://www.rstudio.com/products/rstudio/download/preview/). In the following I assume that you use [***RStudio***](http://www.rstudio.com/) as your ***IDE***.

We can proceed by opening [***RStudio***](http://www.rstudio.com/).

### Installing addmtoolbox from github

To install the ***addmtoolbox*** package from github, first the ***devtools*** package must be loaded in R.
This can be easily done by typing ```install.packages("devtools")``` in the ***R Console*** to install the package and the typing ```library(devtools)``` to load it.

Now type ```install_github('AlexanderFengler/addmtoolbox', build_vignettes=TRUE)``` to install the ***addmtoolbox*** package and make it available in your ***R Session*** with ```library(addmtoolbox)```.

*You should be able to use the provided functions now.*

## Documentation

Besides the standard documentation you will find in the help window in [***RStudio***](http://www.rstudio.com/), the package come with vignettes, one of which carries you through a complete model-fit. 

Be sure to include the ```build_vignettes = TRUE``` parameter when installing the package, otherwise the following instruction will not work.

To list the vignettes of the ***addmtoolbox*** package, type ```vignette(package = "addmtoolbox")```. This lets you see the names of the vignette-documents that come with package. Most importantly you are going to see **addmtoolbox_modelfit_walkthrough** vignette.

To open the vignette, type ```vignette("addmtoolbox_general_info")```. This vignette will carry you through a model-fit to get you started with the package.

## Useful Background
---
Below a list of ***R packages*** that are used in the source code. In case you would like to really understand and contribute to the codebase, I think it will make it considerably easier if you first familiarize yourself with the basic utilities these packages offer (in case you have not already), as they are used quite heavily.

1. [dplyr](http://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html) (for summarising data.frames)
2. [data.table](http://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.pdf) (for high performance versions of data.frames, especially the ***setkey()*** function) 
3. [foreach](http://cran.r-project.org/web/packages/foreach/vignettes/foreach.pdf) (parallel for loops)
4. [Rcpp](http://adv-r.had.co.nz/Rcpp.html) (linking ***c++*** with ***R***) 
5. [RcppZiggurat](http://cran.r-project.org/web/packages/RcppZiggurat/vignettes/RcppZiggurat.pdf) (faster version of ```rnorm()``` )



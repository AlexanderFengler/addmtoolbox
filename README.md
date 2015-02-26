# addmtoolbox

The ***addmtoolbox*** package provides a set of functions that facilitate the process of fitting (attentional) drift diffusion models to epxerimental data.

The basic drift diffusion process is implemented in a highly optimized ***c++*** function, which is complemented with a host of ***R*** functions that help with *data preparation*, facilitation of *parameter optimization* and *parallelization*.

The package is in an early development stage.


# Installing the package

The following guide is specifically for Mac OS systems. Windows systems have not yet been tested. On Linux/Ubuntu system, directly follow up with step ***UNKNOWN*** as gsl is often/always included in the distribution. 

## Dependencies
The ***addmtoolbox*** package depends notably on the ***RcppZiggurat*** package which in turn depends on the ***RcppGSL*** package. For the ***RcppGSL*** package to work, we need to install the [GNU GSL scientific library](http://www.gnu.org/software/gsl/). 

Moreover, you need to have [***Xcode***](https://developer.apple.com/xcode/) installed on your system, which if not present can be downloaded for free from the [***Apple Appstore***](https://itunes.apple.com/de/app/xcode/id497799835?mt=12).

## Installing GSL Library
Below is a step by step guide to installing the [GNU GSL scientific library](http://www.gnu.org/software/gsl/).

1. First [download](http://mirrors.ibiblio.org/gnu/ftp/gnu/gsl/) the current version and unpack it if not done automatically. (choose the file which has the highest version number: ***gsl-1.16.tar.gz*** at the time of writing)
2. Go to the terminal and navigate to the directory to which you downloaded the ***gsl library***. (useful terminal commands: cd, ls)
3. Type ```./configure``` wait until the the execution finished and then type ```make``` to build the library. (you can find these instructions in the ***INSTALL*** file inside of the downloaded ***gsl*** folder). The build process takes a while now.
4. Lastly type ```make install```, which will finally install the library on your system. It is likely that you have to use ```sudo make install```, as the default folder which the library is installed in requires superuser access.

By now we can proceed by opening [***RStudio***](http://www.rstudio.com/).

## Installing addmtoolbox from github

To install the ***addmtoolbox*** package from github, first the ***devtools*** package must be loaded in R.
This can be easily done by typing ```install.packages("devtools")``` in the ***R Console*** to install the package and the typing ```library(devtools)``` to load it.

Now type ```install_github('AlexanderFengler/addmtoolbox')``` to install the ***addmtoolbox*** package and make it available in your ***R Session*** with ```library(addmtoolbox)```.

You should be able to use the provided functions now.




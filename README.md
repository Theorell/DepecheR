# DepecheR on GitHub (2019-02-27 edit: see the last section below if you have trouble installing)

DepecheR is an R package for clustering cytometry and single-cell RNA sequencing data according to a minimal set of biological markers. 
This is done using repeated application of a penalized k-means algorithm, which penalizes partitions that depend on a large number of markers. DepecheR tunes the level of penalization by testing a range of penalization levels, thus minimizing the need for assumptions and arbitrary decisions. Technical documentation is found in /man/, but are most easily accessible by installing the package and accessing the help files. Please see the DepecheR wiki pages (https://github.com/Theorell/DepecheR/wiki/Downloading-DepecheR,-getting-it-to-run,-and-accessing-help-files) for additional information on installation, etc.

Very recently, DepecheR has been accepted by BioConductor. Due to BioConductor-internal reasons, the current master branch version is therefore only compatible with the current development version of R, namely 3.6, which means that it will not install on most systems. For this reason, we have created a separate branch that works with R 3.5. The only difference between these two versions is that the dependency for R 3.6 is changed to 3.5 in the DESCRIPTION file, so this hack is easy to obtain for anyone. But if you just want to install the current version, use the following commands: 

Pre-install devtools if not already in the system

```
install.packages("devtools")
```
Now, install DepecheR
```
library(devtools)
install_github("Theorell/DepecheR", ref = "Rversion3.5")
```

Soon, this hack will be unnecessary, as R 3.6 is expected to come out very soon. 

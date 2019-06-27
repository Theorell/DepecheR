# DepecheR is now on Bioconductor! Please go to https://bioconductor.org/packages/release/bioc/html/DepecheR.html to get access to vinjettes, and to install the latest release or development versions. 

DepecheR is an R package for clustering cytometry and single-cell RNA sequencing data according to a minimal set of biological markers. 
This is done using repeated application of a penalized k-means algorithm, which penalizes partitions that depend on a large number of markers. DepecheR tunes the level of penalization by testing a range of penalization levels, thus minimizing the need for assumptions and arbitrary decisions. 
DepecheR however also contains a function suite for statistical interpretation of group differences, and for visualization of said group differences, in addition to population densities, marker expressions, et cetera. 


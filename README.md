# cyclone
## What is cyclone
cyclone is a set of computational methods to assign cell-cycle stage from single-cell transcriptome data. They include a PCA-based method, a random forest based method and a custom-built predictor.

By Antonio Scialdone, Oliver Stegle, John Marioni and Florian Buettner.
Find more details in the accompanying publication [1].


## How to use cyclone
Cyclone is implemented in R and python. If you would like to use the PCA-based Gaussian Naive Bayes Classifier (GNB), have a look at the [ipython notebook](https://github.com/PMBio/cyclone/blob/master/py/demo/demo_cyclone.ipynb) in the demo folder. For the custom-built predictor (pairs method), have a look in the R folder.

## Prerequisites
The python-based predictors (includng GNB) require a working installation of python 2.7 with scipy, matplotlib, h5py and sklearn.

## References
[1] Scialdone A, Natarajan KN, Saraiva LR, Proserpio V, Teichmann SA, Stegle O, Marioni JC and Buettner F, "Computational assignment of cell-cycle stage from single-cell transcriptome data", Methods, 2015, doi:10.1016/j.ymeth.2015.06.021

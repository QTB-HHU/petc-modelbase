# Model of the photosynthetic electron transport chain
## Model
This repository contains the mathematical model described in the paper published by [Ebenhoeh *et al.* 2014](http://dx.doi.org/10.1098/rstb.2013.0223), that has been originally developed in MATLAB. The model comprises a set of seven coupled ordinary differential equations and captures the temporal evolution of major protein complexes and key products of photosynthetic light reaction for model organism alga *Chlamydomonas reinhardtii*. The model has been written using [modelbase](https://github.com/QTB-HHU/modelbase) package.

## Install requirements

* NumPy
* matplotlib
* [modelbase](https://github.com/QTB-HHU/modelbase)

## Sample analysis
A sample analysis is given in the accompanying notebook. We provide a set of prepared experimental protocols that can be easily amended by the user:
* anoxia - simulates fluorescence trace in anoxia induced state transitions.
* dld - simulates fluorescence trace in Dark-Light-Dark experiment.



# Turing patterns on expanding domain
Simulations for Turring patterns on an apically expanding domain.
The details about the models and numerical implementations can be found in our manuscript:
Y. Liu, P. Maini, R. Baker, "Organisation of diffusion-driven stripe formation in expanding domains" (2021), to be submitted to arXiv.
Two models are implemented: the Schnackenberg model, and a model for the CDIMA chemical reaction. 
Written in Matlab. The numerical scheme used is a finite difference discretization for space, implicit-explicit time stepping with Crank-Nicolson.

## Branches
The 'public' branch contains cleaned up codes that simulate the models. If you are not the authors, you should look at this branch.
The 'main' branch additionally contain codes for analysis, auxillary scripts, as well as work-in-progress, it is unlikely to be interesting to anyone else.
Ignore the other branches

## Files

* cdima\_1d.m, cdima\_2d.m: for simulating the CDIMA model
* schnackenberg\_1d.m, schnackenberg\_2d.m: for simulating the Schnackenberg model
* biggerFont.m, plot\_kymograph.m, tightEdge.m: utilities for plotting, called by the other scripts

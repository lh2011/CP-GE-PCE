# CP-GE-PCE
Conformalized gradient-enhanced polynomial chaos expansion for electromagnetic simulations.

This repository contains the code for the masters thesis "Machine Learning with Polynomial Chaos Expansion for Electromagnetic Simulations" conducted at TU Darmstadt in coorporation with the Centrum Wiskunde & Informatica, Amsterdam. Furthermore, it also contains the relevant data that were computed and used in the thesis for figures and tables.

To reproduce the experiments carried out in the thesis, the following files can be run as is:
 - fieldPlots.m to generate the field plots
 - conformalPrediction.py to recreate the conformalization experiments
 - compute_PCE.py for compression and surrogate accuracy experiments


To generate your own data set, the following steps are required:

1. generate_samples.py: specifify the parameter subset, angle-mode and training/validation/test partition sizes
2. generate_samples.m: generates the data (inlcuding sensitivities) by running the motor simulation for every sampled parameter combination

All experiments can then be run using your new dataset. To run the autoencoder-based experiments, autoencoders first have to be trained on the new datasets.

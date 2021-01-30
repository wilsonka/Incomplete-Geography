This repository contains the code for "Estimation of Health and Demographic Indicators with Incomplete Geographic Information" by Wilson and Wakefield.

Setting up the simulation analysis is done through [`cluster-runDA.R`](cluster-runDA.R). Depending on the simulation set-up (jitterered/masking and spatial covariates or no spatial covariates):

| Model Number (Table 2) | jitter | spcov  |
|  :---: | ------------- | ------------- |
| 3a  | True  | False |
| 6a  | False  | False |
| 3b  | True | True |
| 6b  | False | True|

The main function is [`run-DA()`](run-DA.R#L110), which will run the entire INLA within MCMC algorithm. This function calls [`run_one()`](run-DA.R#L46), which runs a single iteration of the MCMC. First, [`determinep.c()`](run-DA.R#L28) is used to calculate the conditional probability of the potential locations (eqn 8 in the paper). Next, a location is sampled using [`updateloc.c()`](run-DA.R#L14). There is some preprocessing before fitting the INLA model and obtaining samples from the approximate posterior for the model parameters using functions found in [`binomial-model.R`](binomial-model.R).

There are some other useful functions for setting up the simulation:
- [`coverings-on-cluster.R`](coverings-on-cluster.R) includes code for obtaining the normalization factors for each enumeration area (eqn 2 in the paper)
- [`cluster-prior.R`](cluster-prior.R) includes code for creating what is called `Alist` which contains the potential locations for each jittered or masked points. It is used in `determinep.c()`.

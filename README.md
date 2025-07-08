# MathematicalModelsOfCoagulation-AUnifiedApproach
 Supplemental code for the paper "Mathematical Models of Coagulation - A Unified Approach""

# Compiling the Unified Model ODE
In order to speed up the time to solve the model, we recommend compiling the ODE right hand side function (unifiedODE). This can be done with the following MATLAB command:
```Matlab
codegen unifiedODE -args {rand(79,1), rand(109,1)}
```
This generates the compiled code in a mex function named `unifiedODE_mex`. To point the ODE solver to this, edit line 5 in unifiedODESolver.m to point to `unifiedODE_mex` rather than `unifiedODE`

# Generating the data
The raw data, and the functions for importing it into Matlab are in the `data` directory. Running the `dataToIC` script or the `dataToOD` script will load the raw data, handle any missing values and save the initial conditions and OD curves to `IC.mat` and `od.mat`.

# Running all the code
The code can be run, from start to finish by using the `runall.m` script.

# Cross validation
While the code that was used to generate the cross validation, and the posteriors for the cross validation are included, the data files used require patient level demographic information which is not included in the repo.
For this reason, this repo does not generate test and training scores for the different cross validation folds.

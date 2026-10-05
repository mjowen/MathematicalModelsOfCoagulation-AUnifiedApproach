# MathematicalModelsOfCoagulation-AUnifiedApproach
 Supplemental code for the paper "Mathematical Models of Coagulation - A Unified Approach""

# Add to path
All folders should be added to the MATLAB path. This can be done by selecting each folder, right-clicking, clicking "Add to Path" and "Selected Folders".

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

# Make changes to the code
The following descriptions are to inform users what sections of the code need to be updated and how in order to swap out parts of the fitting process.

## Using different prior distributions
The prior distributions are stored in `parameters/priors.mat` as two vectors `priorMean` and `priorSTD` reporting mu and sigma from the distribution `10^N(mu, sigma)` where `N` is a normal distribution.
Each reaction is indexed and ordered acording to the `Prior Form Name` column of `parameters/Draft Parameters.xlsx`. The names in `Reduced Prior Form Name` match the forrest plot figure of the paper.

To change the values, simply change the variables `priorMean` and `priorSTD` in the `parameters/priors.mat` file.
Once this is complete, the files `parameters.mat` (which stored the median of the prior distributions in both prior and mass action form) and `percentiles.mat` (which stores the 5th and 95th percentiles of the prior distributions in prior form) need to be updated.

## Using different data
### ...that is OD data
You can change the OD data (or the initial conditions corresponding to the OD data) in the file `data/pramis_data.xlsx`.
This file contains two sheets, the first for OD data, where each column represents an individual (IDs in the first row) measuring OD at each timepoint.
The second sheet contains initial condition data (ICs). Here, each row represents an individual (IDs in the first column), measuring each of the initial conditions with the labelled units.
This data is then processed into `data/IC.mat` by `data/dataToIC.m` and `data/od.mat` by `data/dataToOD.m`. 

If you are changing the meaning or number of the columns for the initial conditions, you will also need to update `setIC.m` which generates the initial condition vector to be used by the model.

If you change the number of individuals, then some hardcoded values in `evaluateUnifiedModel.m` need to be updated.

### ...that is not OD data
If you want to use data that is not OD data, you must update as above (ideally still storing your data in the `od.mat` file) and then update the `scoreParam.m` function which calculates the cost of the parameters.

For example, if you wanted to fit to ETP data, you could remove all OD values apart from the last one, and then change `scoreParam.m` to only consider the final timepoint.

## Using a different model
The model is stored in `unifiedODE.m` and solved by calling `unifiedODESolver.m`.
If you are changing the model, you likely also need to change the parameters/prior distributions, as discussed about above. You will also need to change the function `parameters/prior2param.m` to convert sampled prior form parameters to mass action parameters.
You will also need to change the function `setIC.m` which specifies how the `IC.mat` data is used to produce a full model IC vector.

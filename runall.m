% Add this directory and all subdirectories to the path
addpath(genpath("."))

% Generate the IC and od data
run("data/dataToIC")
run("data/dataToOD")

% Perform fitting
% Warning: this is very expensive to run, and will not output the same
% posteriors as we generated (different seed). However, the script that
% calls the fitting is shown here.
%run("ABCSMCFitting")

% Evaluate the unified model after fitting
% This will load the posteriors we found, not the files generated in the
% above step.
run("evaluateUnifiedModel")

% Evaluate the previous models (Table 1)
run("previousModels/evaluatePreviousModels")

% Uncertainty quantification
run("uncertainty/densityPlotOverPriors")

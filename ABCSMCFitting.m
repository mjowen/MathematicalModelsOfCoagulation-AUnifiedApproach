parpool
load('IC.mat')
load('od.mat')
load('parameters.mat')
load('priors.mat')

% Used to control which parameters to fit. Here, we fit all parameters
% except those for substrate activation. Indices are from the prior form of
% the parameters
ki = 1:84;

modelSolver = @unifiedODESolver;

N = 2000;
tolEnd = 10;
tolFirst = 250;
U = 200;

id = 1:333; % Option to control cross validation (select which individuals to use for training)
[pop,w,scores] = SMCABCearlyrejection(ki,kPrior,priorMean,priorSTD,IC(id,:),od(id,:),N,U,tolFirst,tolEnd,modelSolver);

save('SMC-Out.mat','pop','w','scores')

poolobj = gcp('nocreate');
delete(poolobj);

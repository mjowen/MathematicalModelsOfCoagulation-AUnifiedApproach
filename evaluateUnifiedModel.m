load('posterior.mat')
load('parameters.mat')
k = kPrior;
k(1:84) = mode(nextPop);
modelSolver = @unifiedODESolver;
evaluateFit(prior2param(k),1:333,modelSolver,true)

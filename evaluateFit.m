function evaluateFit(k,ids,modelSolver,etpCheck)

load('IC.mat','IC')
load('od.mat','od')
trainingIC = IC(ids,:);
trainingOD = od(ids,:);
testIC = IC(setdiff(1:333,ids),:);
testOD = od(setdiff(1:333,ids),:);

scoreTrain = scoreParam(k,trainingIC,trainingOD,modelSolver);
scoreTest = scoreParam(k,testIC,testOD,modelSolver);

disp(['IDs for training: ',num2str(ids)])
disp(['Training Score: ',num2str(scoreTrain)])
disp(['Test Score: ',num2str(scoreTest)])

if etpCheck
    etp = od(:,end);
    modETP = zeros(size(etp));
    
    [c0, odIdx] = setIC([15e-12, 1.4e-6, 2e-8, 1e-8, 7e-10, 9e-8, 1.6e-7, 3e-8, 3.4e-6, 2.5e-9, 8.73e-6]);
    
    [~,~,~,sol] = modelSolver(k,c0,20*60,10^(-14));
    if (sol.x(end)<1200)
        etpPooled = NaN;
    else
        etpPooled = deval(sol,1200,odIdx);
    end
    
    for i = 1:length(od)
        c0 = setIC(IC(i,:));
        [t,~,ODcurve,~] = modelSolver(k,c0,20*60,10^(-14));
        
        if t(end)==1200
            modETP(i) = ODcurve(1,end)/etpPooled*100;
        else
            modETP(i) = NaN;
        end
    end
    disp("ETP Results")
    disp("RMSE")
    a = fitlm(modETP, etp, intercept=false);
    disp(a.RMSE)
    disp("R^2")
    a = fitlm(modETP, etp);
    disp(a.Rsquared.Ordinary)
    disp("p-value")
    disp(a.Coefficients.pValue(2))
end

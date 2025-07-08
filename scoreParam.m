function score = scoreParam(k,IC,OD,modelSolver)

numPatients = size(IC,1);
%For patient in IC
timePoints = 0:30:1200; %Timepoints (in seconds) where the measurements are taken
score = 0;


[c0, odIdx] = setIC([15e-12, 1.4e-6, 2e-8, 1e-8, 7e-10, 9e-8, 1.6e-7, 3e-8, 3.4e-6, 2.5e-9, 8.73e-6]);

% Find pooled ETP to rescale model predicted OD curves
[~,~,~,sol] = modelSolver(k,c0,20*60,10^(-14));
if (sol.x(end)<1200)
    score = NaN;
else
    pooledOD = deval(sol,1200,odIdx);
end

for i = 1:numPatients
    [c0, odIdx] = setIC(IC(i,:));
    
    [~,~,~,sol] = modelSolver(k,c0,20*60,10^(-14));
    % If model failed to solve, abort parameters
    if (sol.x(end)<1200 || isnan(score))
        score = NaN;
        break
    else
        modelOD = deval(sol,timePoints,odIdx);
        modelOD = modelOD/pooledOD*100;
        score = score + sum((modelOD-OD(i,:)).^2);
    end

end
% Rescale to get final score
score = sqrt(score/numPatients);

end

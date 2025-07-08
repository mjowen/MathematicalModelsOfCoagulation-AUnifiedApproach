function [pop,w,scores] = SMCABCearlyrejection(ki,k50,priorMean,priorSTD,IC,OD,N,U,tolFirst,tolEnd,modelSolver)

pop = zeros(1,N,length(ki));
w = zeros(1,N);
scores = zeros(1,N);

if length(OD) == 333
    filename = 'saved_pop.mat';
    disp('Using full data')
else
    r = num2str(rand());
    filename = ['saved_pop_CV_', r, '.mat'];
    disp(['Using cross validation data. Random Name: ', r])
end

counter = 0;
parfor i = 1:N
    warning('off','MATLAB:ode23tb:IntegrationTolNotMet')
    warning('off','MATLAB:nearlySingularMatrix')
    warning('off','MATLAB:singularMatrix')
    warning('off','MATLAB:illConditionedMatrix')
    while true
        counter = counter + 1;
        params = prior(priorMean,priorSTD); %First sample from prior
        k = k50;
        k(ki) = params(ki);
        score = scoreParam(prior2param(k), IC, OD, modelSolver);
        if score<tolFirst
            pop(1,i,:) = k(ki);
            scores(1,i) = score;
            w(1,i) = 1/N;
            disp(num2str(i))
            break
        end
    end
end

disp(['Primary Acceptance Rate: ',num2str(N/counter)])
tol = max(scores); %Current Tolerance
tols = tol;
t = 1;
disp(['Unique Particles at Start: ', num2str(length(unique(scores(t,:))))])
while tol>tolEnd
    disp(['Beginning population ', num2str(t)])
    v = rand(1,N);
    if U > resampling(tol,w(t,:),scores(t,:),v)
        disp('Unable to meet U requirements')
    end
    [tol,idx] = UBisection(tol,w(t,:),scores(t,:),v,U);
    tols = [tols tol];
    disp(['Current Tolerance Bound: ', num2str(tol)])
    disp(['Unique Particles: ', num2str(length(unique(scores(t,idx))))])
    disp(['Desired U: ',num2str(U)])
    
    currentPop = reshape(pop(t,idx,:),[N,length(ki)]);
    currentWeights = w(t,idx)/sum(w(t,idx));
    currentScores = scores(t,idx);
    nextPop = zeros(size(currentPop));
    nextScores = zeros(size(currentScores));
    moveCounter = 0;
    acceptCounter = 0;
    parfor i = 1:N
        particle = currentPop(i,:);
        k = k50;
        k(ki) = particle;
        trialParticle = K(particle,currentPop);
        trialK = k50;
        trialK(ki) = trialParticle;
        %disp(['Move Ratio',num2str(prior(trialK,ki,priorMean,priorSTD)/prior(k,ki,priorMean,priorSTD))])
        if rand < prior(trialK,ki,priorMean,priorSTD)/prior(k,ki,priorMean,priorSTD)
            moveCounter = moveCounter + 1;
            score = scoreParam(prior2param(trialK),IC,OD,modelSolver);
            if score<tol
                particle = trialParticle;
                acceptCounter = acceptCounter + 1;
            else
                score = currentScores(i);
            end
            nextPop(i,:) = particle;
            nextScores(i) = score;
        else
            nextPop(i,:) = particle;
            nextScores(i) = currentScores(i);
        end
    end
    disp(['Move Rate (%): ',num2str(moveCounter/N*100)])
    disp(['Accept Rate (%): ',num2str(acceptCounter/N*100)])
    pop = cat(1,pop, reshape(nextPop,[1,N,length(ki)]));
    scores = cat(1,scores, nextScores);
    w = cat(1,w, currentWeights);
    t = t+1;
    disp('---------')
    disp('---------')
    save(filename)
end
end

function out = K(particle,currentPop)

covMat = 0.5*cov(log10(currentPop));
out = 10.^(mvnrnd(log10(particle),nearestSPD(covMat)));

end

function [tol,idx] = UBisection(tolOld,weights,scores,v,U)
% Find the best tolerance threshold that will give U (or close to U) unique
% particles/scores
tolLower = max(mink(scores,2)); %Tolerance lower bound
tolUpper = tolOld; %Tolerance upper bound
tol = (tolLower+tolUpper)/2; %Test value of tolerance

while (tolUpper-tolLower)/2 > 1e-10 %Run until value of tolerance if found in range of +- error=1e-10
    if sign(resampling(tol,weights,scores,v)-U) == sign(resampling(tolLower,weights,scores,v)-U)
        tolLower = tol;
    else
        tolUpper = tol;
    end
    tol = (tolLower+tolUpper)/2;
    if resampling(tol,weights,scores,v)==U
        break
    end
end

disp(['UBisection: Unique Scores: ',num2str(resampling(tol,weights,scores,v))])
disp(['UBisection: Desired U: ',num2str(U)])
disp('Scores:')
disp(sort(scores))

[~,idx] = resampling(tol,weights,scores,v);

end

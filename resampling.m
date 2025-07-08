function [U,idx] = resampling(tol,weights,scores,v)
% Sample the scores that are less than tol randomly (using the random
% vector v), weighted by the weights
weights(scores>tol) = 0;
weights = weights/sum(weights);
cumWeight = cumsum(weights);
idx = zeros(size(v));
for i = 1:length(v)
    idx(i) = find(cumWeight>v(i),1);
end

U = length(unique(scores(idx))); %Ensure U unique scores in place of U unique particles

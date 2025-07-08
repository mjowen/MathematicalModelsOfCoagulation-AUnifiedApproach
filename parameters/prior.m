function out = prior(varargin)
if nargin == 2
    % Sample from parameters
    mean = varargin{1};
    std = varargin{2};
    out = 10.^normrnd(mean,std);
else
    % Evaluate pdf of prior distribution
    k  = varargin{1};
    ki = varargin{2};
    priorMean = varargin{3};
    priorSTD = varargin{4};
    if sum(k<=0)>0
        out = 0;
    else
        out = prod(normpdf(log10(k(ki)),priorMean(ki),priorSTD(ki)));
    end
end

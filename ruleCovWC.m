function [probF unWeight] = ruleCovWC(feature, candidate, w)
    n = size(feature, 1); % all examples
    tW = sum(w); % total weight
    temp = repmat(candidate, n, 1); % repeat matrix
    res = abs(feature-temp); % feature difference, same feature will be zero after minus
    unWeight = find(sum(res, 2) == 0); % matched feature pairs
    inter = sum(w(unWeight)); % sum of weighted sample
    K = 2; % two class problem
    % probF = (inter+1)/(n+class); % probability of this feature pair using Laplace estimate
    if tW == 0
        probF = (inter+1)/(tW+K); % probability of this feature pair using Laplace estimate
    else
        probF = inter/tW;
    end
end

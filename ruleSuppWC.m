function [probF inter2] = ruleSuppWC(feaLabel, candidate, w)
    n = size(feaLabel, 1); % all examples
    temp = repmat(candidate, n, 1); % repeat matrix
    res = abs(feaLabel-temp); % feature difference, same feature will be zero after minus
    unWeight = find(sum(res, 2) == 0); % matched feature pairs
    inter = sum(w(unWeight)); % sum of weighted sample
    feature = feaLabel(:, 1:end-1); % only features
    cond = repmat(candidate(1, 1:end-1), n, 1); % feature condition
    res2 = feature-cond;
    unWeight2 = find(sum(res2, 2) == 0);
    inter2 = sum(w(unWeight2)); % no. of times feature appears, n(cond)
    K = 2; 
    % probF = (inter+1)/(inter2+class); % probability of this feature pair using Laplace estimate
    if inter2 == 0
        probF = (inter+1)/(inter2+K); % probability of this feature pair using Laplace estimate
    else
        probF = inter/inter2;
    end
    
end

function [probF inter] = ruleSP(feaLabel, candidate)

n = size(feaLabel, 1); % all examples
temp = repmat(candidate, n, 1); % repeat matrix
res = abs(feaLabel-temp); % feature difference, same feature will be zero after minus
inter = length(find(sum(res, 2) == 0)); % no. of times feature and this lable appear, n(cond,class)
feature = feaLabel(:, 1:end-1); % only features
cond = repmat(candidate(1, 1:end-1), n, 1); % feature condition
res2 = abs(feature-cond);
inter2 = length(find(sum(res2, 2) == 0)); % no. of times feature appears, n(cond)
K = 2; 
% probF = (inter+1)/(inter2+class); % probability of this feature pair using Laplace estimate
if inter2 == 0
    probF = (inter+1)/(inter2+K); % probability of this feature pair using Laplace estimate
else
    probF = inter/inter2;
end

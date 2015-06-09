function [probF index] = ruleCov(feature, candidate)

n = size(feature, 1); % all examples
temp = repmat(candidate, n, 1); % repeat matrix
res = abs(feature-temp); % feature difference, same feature will be zero after minus
index = find(sum(res, 2) == 0); % matched feature pairs
inter = length(find(sum(res, 2) == 0)); % how many common paris found
K = 2; % two class problem
% probF = (inter+1)/(n+class); % probability of this feature pair using Laplace estimate
if n == 0
    probF = (inter+1)/(n+K); % probability of this feature pair using Laplace estimate
else
    probF = inter/n;
end

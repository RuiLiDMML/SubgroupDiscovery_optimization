function [ruleDataCov ruleDataSupp] = basRule(data, label, class)
%% compute basic rule coverage and precision

if nargin < 3
    class = 1;
end

uniLabel = unique(label);
ruleDataCov = zeros(size(data, 1), size(data, 2));
ruleDataSupp = ruleDataCov;

for i = 1:size(data, 2)
    fea = data(:, i);
    uni = unique(fea);
    for j = 1:length(uni)
        candiFea = uni(j);
        [gF fIdx] = ruleCov(fea, candiFea); % coverage
        [pF count] = ruleSP([fea label], [candiFea uniLabel(class)]); % support
        ruleDataCov(fIdx, i) = gF;
        ruleDataSupp(fIdx, i) = pF;
    end
end











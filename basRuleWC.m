function [ruleDataCov ruleDataSupp] = basRuleWC(data, label, class, weight)
%% compute basic rule coverage and precision using weighted covering
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
        [gF fIdx] = ruleCovWC(fea, candiFea, weight); % coverage
        pF = ruleSuppWC([fea label], [candiFea uniLabel(class)], weight); % support
        ruleDataCov(fIdx, i) = gF;
        ruleDataSupp(fIdx, i) = pF;
    end
end


end
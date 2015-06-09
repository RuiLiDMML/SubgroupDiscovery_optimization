function [linear] = covSuppWC(data, label, class, method, weight, feaWeight)
%%% coverage and support computation using weighted covering, i.e. sample
%%% weigths, adding a feature weights computed from SVM

[n dim] = size(data);

if strcmp(method, 'equal')
    weight = ones(n, 1); % samples are equal
end

uniLabel = unique(label);
idx1 = find(label == uniLabel(1));
nPos = sum(weight(idx1));

idx2 = find(label == uniLabel(2));
nNeg = sum(weight(idx2));

p0F1 = nPos/sum(weight); % fraction of positive targets
p0F2 = nNeg/sum(weight);
p0FPN = [p0F1; p0F2];

p0F = p0FPN(class);
p0Fmat = p0F.*ones(n, dim);

[cover supp] = basRuleWC(data, label, class, weight);

p0F = p0FPN(class);
p0Fmat = p0F.*ones(n, dim);
pgMat = (cover.*(supp>p0FPN(class))).^feaWeight.*(supp-p0Fmat);
pgMat(isinf(pgMat)) = 0;
linear = zeros(dim, 1);
for i = 1:dim
    uni = unique(pgMat(:, i));
    quali = uni(find(uni)); % qualified feature values
    linear(i, 1) = sum(quali);
    clear uni quali
end











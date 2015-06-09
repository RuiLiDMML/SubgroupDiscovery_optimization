function [dependency NMI] = muInfoFea(data)
% mutual information between feature pair
% dependency: symmetric matrix represents feature denpendency

[n p] = size(data);
dependency = zeros(p, p);

for i = 1:p
    for j = i:p
        fea1 = data(:, i);
        fea2 = data(:, j);
        depen = muInfo(fea1, fea2);
        dependency(i, j) = depen;
        dependency(j, i) = depen;
    end
end

% normalised mutual information
NMI = zeros(p, p);
for i = 1:p
    for j = i:p
        NMI(i, j) = dependency(i, j)/sqrt(dependency(i, i)*dependency(j, j));
        NMI(j, i) = NMI(i, j);
    end
end

NMI(isnan(NMI)) = 0;
NMI(isinf(NMI)) = 0;
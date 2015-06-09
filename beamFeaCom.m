function [comb] = beamFeaCom(beamFea, beamValue, dimension, data)
%%% given beam features, output next level candidate features, just by
%%% adding one feature

[m dim] = size(beamFea);
feaComb = [];
%% feature combination
ct = 0;
for i = 1:size(beamFea, 1)
    fea = beamFea(i, :);
    index = 1:dimension;
    index(fea) = []; % remove features in beam
    for j = 1:length(index)
        ct = ct + 1;
        feaComb(ct, :) = [fea index(j)];
        feaCombValue(ct, :) = beamValue(i, :);
    end
end

%% feature value combination
combs = [];
uniFea = cell(dim+1, 1);
if ~isempty(feaComb)
    for i = 1:size(feaComb, 1)
        feaC = feaComb(i, :);
        targetF = feaC(end); % the feature to be added
        targetFV = data(:, targetF); % feature values
        uniTFV = unique(targetFV);
        for j = 1:dim
            uniFea{j, 1} = feaCombValue(i, j);
        end
        uniFea{dim+1, 1} = uniTFV;
        combsF = feaUniComb(uniFea, dim+1);
        guide = repmat(feaC, size(combsF, 1), 1); % feature pairs
        combs = [combs; [combsF guide]]; % first half col. are feature values and the last half col. are the resp. feature indices

        clear uniTFV combsF guide
    end
end
%% check repetation
[allCombF, allCombFIndex] = repRemove(combs);
comb = [allCombF allCombFIndex];


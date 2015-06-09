function [quad linear] = dataInfoQP(data, label, class, method)
%%% compute H and f for QP
% quad: quadratic term, a symmetric matrix
% linear: linear term, a vector
[n dim] = size(data);
[cover supp] = basRule(data, label, class);
a = 1; % pg score parameter

uniLabel = unique(label);
nPos = length(find(label == uniLabel(1)));
nNeg = length(find(label == uniLabel(2)));

p0F1 = nPos/n; % fraction of positive targets
p0F2 = nNeg/n;
p0FPN = [p0F1; p0F2];

switch method
    case 'MIpg'
        %% mutual information between feature pair
        quad = muInfoFea(data);
%         quad(logical(eye(size(quad))))=0; % set diagnoal to 0, modified on 22.05.2012
        %% linear term: pg score
        [ruleDataCov ruleDataSupp] = basRule(data, label, class);
        p0F = p0FPN(class);
        p0Fmat = p0F.*ones(n, dim);
        pgMat = (ruleDataCov.*(ruleDataSupp>p0F)).^a.*(ruleDataSupp-p0Fmat);
        linear = zeros(dim, 1);
        for i = 1:dim
            uni = unique(pgMat(:, i));
            quali = uni(find(uni)); % qualified feature values
            linear(i, 1) = sum(quali);
            clear uni quali
        end
        
        
    case 'entropy' % use mutual information based method
        %% mutual information between feature pair
        quad = muInfoFea(data);
        %% feature value importance (support)
        feaImp = zeros(dim, 1);
        for i = 1:dim
            imp = unique(supp(:, i));
            impI = find(imp > 0.5); % feature values predict the studied class
            impV = imp(impI);
            feaImp(i, 1) = sum(impV);
            quad(i, i) = -sum(impV); 
            coverV = cover(impI, i); % qualified coverage
            feaCov(i, 1) = sum(coverV);
        end
    %     linear = feaImp./sum(feaImp);
        linear = feaCov./sum(feaCov);
    
    case 'CS' % coverage and support method
        feaImp = zeros(dim, 1);
        quad = zeros(dim, dim); % support, coverage matrix info.
        linear = zeros(dim, 1); % support info.
        coverSum = zeros(dim, 1); % coverage sum
        values = cell(dim, 1);
        for i = 1:dim
            imp = unique(supp(:, i));
            impI = find(imp > 0.5); % indices of feature values predict the studied class
            impV = imp(impI); % feature values predict the studied class
            impVdata = zeros(length(impV), 1);
            if ~isempty(impV)
                for j = 1:length(impV)
                    temp = find(supp(:, i) == impV(j));
                    impVdata(j, 1) = data(temp(1), i); % need only the first index pointing to the value
                end
                values{i, 1} = impVdata; % save qualified values in a cell
                clear temp impVIdx
            end
            
            feaImp(i, 1) = sum(impV); % vector of support values
            coverV = cover(impI, i); % qualified coverage
            coverSum(i, 1) = sum(coverV);
            linear(i, 1) = sum(impV);
        end      

        
        %% pair-wise coverage intersection
        for i = 1:dim-1
            tar = values{i, 1};
            fea = data(:, i);
            if ~isempty(tar)
                for j = i+1:dim
                    tar2 = values{j};
                    fea2 = data(:, j);
                    if ~isempty(tar2)
                        len = length(tar);
                        len2 = length(tar2);
                        %%% all possible value combinations
                        count = 0;
                        for p = 1:len
                            for q = 1:len2
                                count = count + 1;
                                idx = (fea == tar(p)); % logical index of target value
                                idxSum = sum(idx); % how many are found
                                idx2 = (fea2 == tar2(q));
                                idxSum2 = sum(idx2);
                                inter = idx+idx2;
                                index = find(inter == 2); % if appear together then must be 2
                                frac(count, 1) = length(index)/n; % fraction of intersected value pairs
                            end 
                        end
                    else
                        frac = 0;
                        %%% end all possible value combinations
                    end
                    quad(i, j) = sum(frac);
                    quad(j, i) = sum(frac);
                end
            end
        end
        quad(logical(eye(size(quad)))) = coverSum; % assign diagonal    
    end

end













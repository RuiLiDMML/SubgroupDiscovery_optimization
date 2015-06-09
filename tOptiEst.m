function SDrule = tOptiEst(data, label, thres)
%%% subgroup discovery by optimistic estimate
SDrule = [];
if nargin < 3
    thres = 0.1;
end
[n dim] = size(data);
uniLabel = unique(label);
nPos = length(find(label == uniLabel(1)));
nNeg = length(find(label == uniLabel(2)));
maxDepth = 4; % maximal feature number for combination

a = 1; % pg score parameter
p0F1 = nPos/n; % fraction of positive targets
p0F2 = nNeg/n;
p0FPN = [p0F1; p0F2];

for c = 1:2 %% for pos and neg class
    sampCovCount = zeros(n, 1); % sample cover counts
    ct = 0;
    ct2 = 0;
    p0F = p0FPN(c);
    counter = 0;
    SDrules = [];
    pgS = [];
    pgOE = [];
    pgOEIndex = [];
    val = [];
    feas = [];
    %% a single feature, 1st round
    flag = 'false';
    for i = 1:dim
        fea = data(:, i); 
        uni = unique(fea);
        for j = 1:length(uni)
            candiFea = uni(j);
            [gF idx]= featureIntersect(fea, candiFea); % coverage
            [pF count] = featureLabIntersect([fea label], [candiFea uniLabel(c)]); % support
            pgTemp = gF^a*(pF-p0F);

           %%% compute rule significance, it considers both classes
           if c == 1
               [pFS countS] = featureLabIntersect([fea label], [candiFea uniLabel(2)]); % support of class2
               sig1 = count*log2((count+1)/((nPos*gF)+2)); % class positive part
               sig2 = countS*log2((countS+1)/((nNeg*gF)+2)); % class negative part
           else
               [pFS countS] = featureLabIntersect([fea label], [candiFea uniLabel(1)]); % support of class1
               sig1 = count*log2((count+1)/((nNeg*gF)+2)); % class positive part
               sig2 = countS*log2((countS+1)/((nPos*gF)+2)); % class negative part
           end
           sig = 2*(sig1+sig2);
            
            %%% qualified rules
            if pgTemp > thres
                ct = ct + 1;
                counter = counter + 1; % global counter
                pgS(counter, 1) = pgTemp;
                rule{counter, 1} = candiFea; % feature values
                ruleFea{counter, 1} = i; % features used
                ruleCov(counter, 1) = gF; % coverage
                ruleSupp(counter, 1) = pF; % support
                ruleSig(counter, 1) = sig;
                sampleCover(counter, 1) = length(idx);
                sampCovCount(idx, 1) = sampCovCount(idx, 1) + 1;
            end
            %%% tight optimistic estimate
            pgOETemp = (gF*pF)^a*(1-p0F);
            if pgOETemp > thres
                ct2 = ct2 + 1;
                pgOE(ct2, 1) = pgOETemp;
                pgOEIndex(ct2, 1) = i; % which feature
                pgOEIndex(ct2, 2) = candiFea; % feature value
                flag = 'true';
            end
        end
    end
    %% feature pairs
    feaNum = 2;
    pgOEidx = pgOEIndex;
    allCombF = [];
    while strcmp(flag, 'true')
        allCombF = [];
        allCombFIndex = [];
        feas = unique(pgOEIndex(:, 1:feaNum-1));        
        %% generate all possible feature pair combinations
        tarFea = pgOEIndex(:, 1:end/2); % features from last round
        tarFeaUni = unique(tarFea); % features from last round
        tarFeaV = pgOEIndex(:, end/2+1:end); % feature values from last round
        %% all possible feature pairs values, Apriori_gen style
            %%% useful feature and their qualified values
        uniFea = [];
        for j = 1:length(tarFeaUni)
            idx = find(tarFea==tarFeaUni(j));
            uni = unique(tarFeaV(idx));
            temp = repmat(tarFeaUni(j), length(uni), 1);
            add = [temp uni];
            uniFea = [uniFea; add];
        end

        allCombF = [];
        allCombFIndex = [];
        for j = 1:size(pgOEIndex, 1)
            candi = tarFea(j, :);
            candiV = tarFeaV(j, :);
            diff = setdiff(tarFeaUni, candi);
            for z = 1:length(diff)
               new = [candi diff(z)]; % feature pair
               idx = find(uniFea(:, 1) == diff(z));
               val = uniFea(idx, 2); % distinct values
               for p = 1:length(val)
                   all(p, 1:feaNum) = [candiV val(p)];
               end
               allCombF = [allCombF; all]; %% newly formed feature pairs values
               allCombFIndex = [allCombFIndex; repmat(new, length(val), 1)]; % newly formed feature pairs
               clear all new
            end
        end
        %% check repeation
        [allCombF, allCombFIndex] = repRemove([allCombF allCombFIndex] );

        %%% pg score
        ct = 0;
        ct2 = 0;
        pgIndex = [];
        pgOEIndex = [];
        for i = 1:size(allCombF, 1)      
            len = size(allCombFIndex, 2); % feature length
            candiFea = allCombF(i, :); % feature values
            fea = data(:, allCombFIndex(i, :)); % features
            [gF idx] = featureIntersect(fea, candiFea); % coverage
            [pF count] = featureLabIntersect([fea label], [candiFea uniLabel(c)]); % support
            pgTemp = gF^a*(pF-p0F);
           %%% compute rule significance, it considers both classes
           if c == 1
               [pFS countS] = featureLabIntersect([fea label], [candiFea uniLabel(2)]); % support of class2
               sig1 = count*log2((count+1)/((nPos*gF)+2)); % class positive part
               sig2 = countS*log2((countS+1)/((nNeg*gF)+2)); % class negative part
           else
               [pFS countS] = featureLabIntersect([fea label], [candiFea uniLabel(1)]); % support of class1
               sig1 = count*log2((count+1)/((nNeg*gF)+2)); % class positive part
               sig2 = countS*log2((countS+1)/((nPos*gF)+2)); % class negative part
           end
           sig = 2*(sig1+sig2);
            
            if pgTemp > thres                          
                counter = counter + 1; % global counter
                pgS(counter, 1) = pgTemp;
                rule{counter, 1} = candiFea; % feature values
                ruleFea{counter, 1} = allCombFIndex(i, :); % features used
                ruleCov(counter, 1) = gF; % coverage
                ruleSupp(counter, 1) = pF; % support
                ruleSig(counter, 1) = sig;  
                sampleCover(counter, 1) = length(idx);
                sampCovCount(idx, 1) = sampCovCount(idx, 1) + 1;

            end
            %%% optimistic estimate
            pgOETemp = (gF*pF)^a*(1-p0F);
            if pgOETemp > thres
                ct2 = ct2 + 1;
                pgOE(ct2, 1) = pgOETemp;
                pgOEIndex(ct2, 1:len) = allCombFIndex(i, :); % which feature
                pgOEIndex(ct2, len+1:2*len) = candiFea; % feature value
                flag = 'true';
            end
        end
        feaNum = feaNum + 1;
        if isempty(pgOEIndex) || size(allCombFIndex, 2) >= maxDepth
            flag = 'false';
        end
        
    end

    %% sort subgroup description by quality
    if ~isempty(pgS)
        remove = find(pgS < thres); % quality smaller than 0.1, then remove
        pgS(remove) = [];
        [pgS index] = sort(pgS, 'descend');
        
        rule(remove) = [];
        ruleFea(remove) = [];
        ruleCov(remove) = [];
        ruleSupp(remove) = [];
        ruleSig(remove) = [];
        rules = cell(size(pgS, 1), 1); % ordered new rules
        rulesFea = cell(size(pgS, 1), 1); 
        
        len = size(pgS, 1);
        rules = cell(len, 1); % ordered new rules
        rulesFea = cell(len, 1); 
        ruleCover = zeros(len, 1);
        ruleSupport = zeros(len, 1);
        ruleSigf = zeros(len, 1);
        sampleFreq = zeros(len, 1);
        
        for i = 1:size(pgS, 1)
            rules{i, 1} = rule{index(i)}';
            rulesFea{i, 1} = ruleFea{index(i)}';
            ruleCover(i, 1) = ruleCov(index(i));
            ruleSupport(i, 1) = ruleSupp(index(i));
            ruleSigf(i, 1) = ruleSig(index(i));
            sampleFreq(i, 1) = sampleCover(index(i));
        end
        SDrules.pg = pgS;
        SDrules.rule = rules;
        SDrules.feature = rulesFea;  
        SDrules.coverage = ruleCover;
        SDrules.support = ruleSupport;
        SDrules.ruleSignificance = ruleSigf;
        SDrules.sampleFreq = sampleFreq;
        SDrules.sampCovCount = sampCovCount;
        %% remove irrelevant rules
        SDrules = ruleFilter(SDrules); 
        SDrule{1, c} = SDrules;
    end
    clear rule ruleFea ruleCov ruleSupp ruleSig rules rulesFea

end 


%% sub-functions

%%% the probability of feature pairs appear, i.e. how often [1 3] appear in
%%% whole feature set, i.e. coverage cov(R) = n(cond)/n_s
function [probF indexR] = featureIntersect(feature, candidate)
    n = size(feature, 1); % all examples
    temp = repmat(candidate, n, 1); % repeat matrix
    res = abs(feature-temp); % feature difference, same feature will be zero after minus
    indexR = find(sum(res, 2)==0); % matched indices
    inter = length(indexR); % how many common paris found
    nK = 2; % two class problem
    if n == 0
        probF = (inter+1)/(n+nK); % probability of this feature pair using Laplace estimate
    else
        probF = inter/n;
    end
end

%%% the prob. of feature, label pairs appear, i.e. support
function [probF inter] = featureLabIntersect(feaLabel, candidate)
    n = size(feaLabel, 1); % all examples
    temp = repmat(candidate, n, 1); % repeat matrix
    res = abs(feaLabel-temp); % feature difference, same feature will be zero after minus
    inter = length(find(sum(res, 2) == 0)); % no. of times feature and this label appear, n(cond,class)
    feature = feaLabel(:, 1:end-1); % only features
    cond = repmat(candidate(1, 1:end-1), n, 1); % feature condition
    res2 = abs(feature-cond);
    inter2 = length(find(sum(res2, 2) == 0)); % no. of times feature appears, n(cond)
    nK = 2;
    if inter2 == 0
        probF = (inter+1)/(inter2+nK); % probability of this feature pair using Laplace estimate
    else
        probF = inter/inter2;
    end
end

%%% check repeations of next round feature combinations, e.g. features [1 2 3] with
%%% value [2 0 4] is same as features [3 2 1] with value [4 0 2]
function [L R] = repRemove(pg)
    idx = [];
    indexAll = [];
    pghalfL = pg(:, 1:end/2);
    pghalfR = pg(:, end/2+1:end);
    [pghalfR id] = sort(pghalfR, 2);
    %%% reorder pghalfL
    for w = 1:size(id, 1)
        temp = pghalfL(w, :);
        pghalfL(w, :) = temp(id(w, :));
    end
    
    [newmat,index] = unique(pghalfR,'rows','first');  % Finds indices of unique rows
    repeatedIndex = setdiff(1:size(pghalfR,1),index);  % Finds indices of repeats
    repeat = unique(pghalfR(repeatedIndex, :), 'rows');
    [val index] = intersect(pghalfR, repeat, 'rows');
    
    for w = 1:size(index, 1)
        rep = repmat(val(w, :), size(pghalfR, 1), 1);
        res = sum(abs(pghalfR - rep), 2);
        ind = find(res == 0);
        values = pghalfL(ind, :);
        [valT ix] = unique(values, 'rows');
        if length(ix) ~= size(values, 1)
            reIdx = setdiff(1:size(values, 1), ix);
            reduIdx = ind(reIdx); % redudant indicies
            indexAll = [indexAll; reduIdx];
            
        end
    end
    %%% remove redudant ones
    pg(indexAll, :) = [];
    L = pg(:, 1:end/2);
    R = pg(:, end/2+1:end);
end
%% end sub-functions

end
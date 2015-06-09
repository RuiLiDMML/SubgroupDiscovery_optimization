function [SDrule targetIdx] = SDexhauSA(data, label, target, thres, weight)
%%% subgroup discovery by exhaustive search on a single target
%%% evalution quality function: gF^a*(pF-p0F), where a = 1
% a is adaptively computed from the svm weights and each feature has an 'a'

if nargin < 4
    thres = 0.1;
end

method = 'counting sum';

[n p] = size(data);
targetIdx = zeros(n, 1);
uniLabel = unique(label);
nPos = length(find(label == uniLabel(1)));
nNeg = length(find(label == uniLabel(2)));
%% all possible feature combinations
combs = cell(p, 1);
for i = 1:p
   combination = nchoosek(1:p, i);
   combs{i, 1} = combination;
end

p0F1 = nPos/n; % fraction of positive targets
p0F2 = nNeg/n;
p0FPN = [p0F1; p0F2];


if target == 1
    tar = 1;
else
    tar = 2;
end

pgThres = thres;
for c = tar:tar %% for pos and neg class
    sampCovCount = zeros(n, 1);
    ct = 0;
    p0F = p0FPN(c);
    SDrules = [];
    pgS = [];
    for i = 1:p
        combination = combs{i, 1};
        for j = 1:size(combination, 1)
           a = mean(weight(combination(j, :))); % adaptively computed from svm
           fea = data(:, combination(j, :));% feature subsets to study
           dim = length(combination(j, :)); % number of features
           %%% find distinctive values in each feature
           uniFea = cell(dim, 1); % save unique value for each feature
           for p = 1:dim
               singleFea = fea(:, p);
               uniFea{p, 1} = unique(singleFea);
           end % match k = 1:dim
           if dim ~= 1
               combsF = feaUniComb(uniFea, dim); % all combinations of this feature subsets
           else
               combsF = uniFea{1, 1};
           end
           %% subgroup discovery evaluation based on the features
           for k = 1:size(combsF, 1)
               ct = ct + 1;
               candiFea = combsF(k, :); % feature paris
               [gF idx] = featureIntersect(fea, candiFea); % coverage
               [pF count] = featureLabIntersect([fea label], [candiFea uniLabel(c)]); % support
               pgTemp = gF^a*(pF-p0F);

               %%% compute rule significance, it considers both classes
               if c == 1
                   [pFS countS] = featureLabIntersect([fea label], [candiFea uniLabel(2)]); % support of class2
                   sig1 = count*log2((count)/((nPos*gF))); % class positive part
                   sig2 = countS*log2((countS)/((nNeg*gF))); % class negative part
               else
                   [pFS countS] = featureLabIntersect([fea label], [candiFea uniLabel(1)]); % support of class1
                   sig1 = count*log2((count)/((nNeg*gF))); % class positive part
                   sig2 = countS*log2((countS)/((nPos*gF))); % class negative part
               end
               sig = 2*(sig1+sig2);
               
               if pgTemp > pgThres
                   if strcmp(method, 'normal')
                       targetIdx(idx, 1) = 1; % mark found indices as 1
                   else
                       % counting sum, weighted covering style
                       targetIdx(idx, 1) = targetIdx(idx, 1) + 1; % mark found indices as 1
                   end
                   pgS(ct, 1) = gF^a*(pF-p0F); % quality score
                   rule{ct, 1} = candiFea'; % feature values
                   ruleFea{ct, 1} = combination(j, :)'; % features used
                   ruleCov(ct, 1) = gF; % coverage
                   ruleSupp(ct, 1) = pF; % support
                   ruleSig(ct, 1) = sig; 
                   sampleCover(ct, 1) = length(idx);
                   sampCovCount(idx, 1) = sampCovCount(idx, 1) + 1;
               else
                   ct = ct - 1;
               end
           end
        end % match j = 1:size(combination, 1)
    end % match i = 1:p
    %% sort subgroup description by quality
    if ~isempty(pgS)
        remove = find(pgS < thres); % quality smaller than thres, then remove
        pgS(remove) = [];
        [pgS index] = sort(pgS, 'descend');
        rule(remove) = [];
        ruleFea(remove) = [];
        ruleCov(remove) = [];
        ruleSupp(remove) = [];
        ruleSig(remove) = [];
        len = size(pgS, 1);
        rules = cell(len, 1); % ordered new rules
        rulesFea = cell(len, 1); 
        ruleCover = zeros(len, 1);
        ruleSupport = zeros(len, 1);
        ruleSigf = zeros(len, 1);
        sampleFreq = zeros(len, 1);
        
        for i = 1:size(pgS, 1)
            rules{i, 1} = rule{index(i)};
            rulesFea{i, 1} = ruleFea{index(i)};
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
        SDrule{1, 1} = SDrules;
    else
        SDrule{1, 1} = 'empty';
    end
    clear rule ruleFea ruleCov ruleSupp ruleSig rules rulesFea sampleFreq
end


%% sub-functions

%%% the probability of feature pairs appear, i.e. how often [1 3] appear in
%%% whole feature set, i.e. coverage cov(R) = n(cond)/n_s
function [probF indexR] = featureIntersect(feature, candidate)
    n = size(feature, 1); % all examples
    temp = repmat(candidate, n, 1); % repeat matrix
    res = abs(feature-temp); % feature difference, same feature will be zero after minus
    indexR = find(sum(res, 2)==0); % matched indices
    inter = length(find(sum(res, 2) == 0)); % how many common paris found
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
%% end sub-functions



end



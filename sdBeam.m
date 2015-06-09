function SDrules = sdBeam(data, label, thres, beam)
%%% subgroup discovery by beam search

if nargin < 3
    thres = 0.05;
end
 
if nargin < 4
    beam = 15;
end

[n dim] = size(data);
dimension = size(data, 2);
uniLabel = unique(label);
nPos = length(find(label == uniLabel(1)));
nNeg = length(find(label == uniLabel(2)));
 
p0F1 = nPos/n; % fraction of positive targets
p0F2 = nNeg/n;
p0FPN = [p0F1; p0F2];

for i = 1:dim
   combination = nchoosek(1:dim, i);
   combs{i, 1} = combination;
end
 
depth = 4;
 

tol = 10^(-12);
a = 1;
 
for c = 1:2
    sampCovCount = zeros(n, 1); % sample cover counts
    class = c;
    changeFlag = 'true';
    iter = 0;
    ct = 0;
    p0F = p0FPN(c);
    SDrule = [];
    pgS = [];
    
    while iter < depth && strcmp(changeFlag, 'true') == 1
        iter = iter + 1;
        if iter == 1
            %% work for single feature, 1st level
            combination = combs{1, 1};
            for i = 1:dim
                fea = data(:, i);
                uniFea = unique(fea);
                for j = 1:length(uniFea)
                   candiFea = uniFea(j);
                   ct = ct + 1;
                   %%% compute rule quality
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
                       pgS(ct, 1) = pgTemp; % quality score
                       rule{ct, 1} = candiFea; % feature values
                       ruleFea{ct, 1} = i; % features used
                       ruleCov(ct, 1) = gF; % coverage
                       ruleSupp(ct, 1) = pF; % support
                       ruleSig(ct, 1) = sig; 
                       sampleCover(ct, 1) = length(idx);
                       sampCovCount(idx, 1) = sampCovCount(idx, 1) + 1;
                   else
                       ct = ct - 1;
                   end

                end % match j = 1:length(uniFea)
                clear fea uniFea
            end % match i = 1:dim
            
            if ct > beam
               [val idx] = sort(pgS, 'descend');
               value = cell2mat(rule);
               beamValue = value(idx(1:beam)); % keep the best rules feature values
               features = cell2mat(ruleFea);
               beamFea = features(idx(1:beam)); % kepp the best rules feature
               [feaCombBm] = beamFeaCom(beamFea, beamValue, dimension, data);
               clear value features
            else
                if ct ~= 0
                    [feaCombBm] = beamFeaCom(cell2mat(ruleFea), cell2mat(rule), dimension, data);
                end
            end
            if ct == 0
                changeFlag = 'false';
            end
            
        else % match iter == 1
            pgSB = [];
            ruleB = [];
            ruleFeaB = [];
            ctB = 0;
            %% work for multiple feature >= 2nd level
            combination = feaCombBm(:, end/2+1:end); % feature indices
            combsF = feaCombBm(:, 1:end/2); % feature value
            for j = 1:size(combsF, 1)
               candiFea = combsF(j, :);
               candiFeaIdx = combination(j, :);
               fea = data(:, candiFeaIdx);% feature subsets to study
               ct = ct + 1;
               ctB = ctB + 1;
               %%% compute rule quality
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
                   pgS(ct, 1) = pgTemp; % quality score
                   rule{ct, 1} = candiFea'; % feature values
                   ruleFea{ct, 1} = candiFeaIdx'; % features used
                   ruleCov(ct, 1) = gF; % coverage
                   ruleSupp(ct, 1) = pF; % support
                   ruleSig(ct, 1) = sig;
                   sampleCover(ct, 1) = length(idx);
                   sampCovCount(idx, 1) = sampCovCount(idx, 1) + 1;
                   %%% info. for beam tracking
                   pgSB(ctB, 1) = pgTemp; % quality score
                   ruleB{ctB, 1} = candiFea; % feature values
                   ruleFeaB{ctB, 1} = candiFeaIdx; % features used
               else
                   ct = ct - 1;
                   ctB = ctB -1;
               end

            end % match j = 1:size(combsF, 1)
            
            if ctB > beam
               [val idx] = sort(pgSB, 'descend');
               value = cell2mat(ruleB);
               beamValue = value(idx(1:beam), :); % keep the best rules feature values
               features = cell2mat(ruleFeaB);
               beamFea = features(idx(1:beam), :); % kepp the best rules feature
               [feaCombBm] = beamFeaCom(beamFea, beamValue, dimension, data);
               clear value features
            else
                if ctB ~= 0
                    [feaCombBm] = beamFeaCom(cell2mat(ruleFeaB), cell2mat(ruleB), dimension, data);
                end
            end
            if ctB == 0
                changeFlag = 'false';
            end
            
        end % match iter == 1
    end % match iter <= depth && strcmp(changeFlag, 'true') == 1
    
    %% reorganize subgroup
    if ~isempty(pgS)
        [pgS index] = sort(pgS, 'descend');
        len = size(pgS, 1);
        rules = cell(len, 1); % ordered new rules
        rulesFea = cell(len, 1); 
        ruleCover = zeros(len, 1);
        ruleSupport = zeros(len, 1);
        ruleSigf = zeros(len, 1);
        sampleFreq = zeros(len, 1);
        
        for i = 1:size(pgS, 1)
            rules{i, 1} = rule{index(i)}; % ruleFilter needs row vector, see following
            rulesFea{i, 1} = ruleFea{index(i)};
            ruleCover(i, 1) = ruleCov(index(i));
            ruleSupport(i, 1) = ruleSupp(index(i));
            ruleSigf(i, 1) = ruleSig(index(i));
            sampleFreq(i, 1) = sampleCover(index(i));
        end
        SDrule.pg = pgS;
        SDrule.rule = rules;
        SDrule.feature = rulesFea;  
        SDrule.coverage = ruleCover;
        SDrule.support = ruleSupport;
        SDrule.ruleSignificance = ruleSigf;
        SDrule.sampleFreq = sampleFreq;
        SDrule.sampCovCount = sampCovCount;
        %% remove irrelevant rules
        SDrule = ruleFilter(SDrule); 
        SDrules{1, c} = SDrule;
    else
        SDrules{1, c} = 'empty';
    end
    clear rule ruleFea ruleCov ruleSupp ruleSig rules rulesFea    
    
    
end % match c = 1:2
 
 
 
%% sub-functions
 
%%% the probability of feature pairs appear, i.e. how often [1 3] appear in
%%% whole feature set, i.e. coverage cov(R) = n(cond)/n_s
function [probF indexR] = featureIntersect(feature, candidate)
    n = size(feature, 1); % all examples
    temp = repmat(candidate, n, 1); % repeat matrix
    res = abs(feature-temp); % feature difference, same feature will be zero after minus
    indexR = find(sum(abs(res), 2)==0); % matched indices
    inter = length(find(sum(abs(res), 2) == 0)); % how many common paris found
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
    inter = length(find(sum(abs(res), 2) == 0)); % no. of times feature and this label appear, n(cond,class)
    feature = feaLabel(:, 1:end-1); % only features
    cond = repmat(candidate(1, 1:end-1), n, 1); % feature condition
    res2 = feature-cond;
    inter2 = length(find(sum(abs(res2), 2) == 0)); % no. of times feature appears, n(cond)
    nK = 2;
    if inter2 == 0
        probF = (inter+1)/(inter2+nK); % probability of this feature pair using Laplace estimate
    else
        probF = inter/inter2;
    end
end
%% end sub-functions
 
end
 
 
 
 
 
 
 
 
 
 
 
 
 



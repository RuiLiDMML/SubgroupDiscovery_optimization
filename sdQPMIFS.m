function SDRules = sdQPMIFS(data, label, thres, alpha)
%%% subgroup discovery by quadratic programming using mutual information

[n dim] = size(data);
method = 'MIpg';
indexCheck = zeros(n, 1);
tol = 10^(-3);
if nargin < 4
    alpha = 0.5;
end

maxWidth = 4; % maximal feature number for combination

delta = 0.1;

[n p] = size(data);
uniLabel = unique(label);
nPos = length(find(label == uniLabel(1)));
nNeg = length(find(label == uniLabel(2)));
p0F1 = nPos/n; % fraction of positive targets
p0F2 = nNeg/n;
p0FPN = [p0F1; p0F2];
diffPN = p0F1-p0F2;

[score] = infoGainID3(data, label);
out.weight = abs(score);
% weighting function
newA = abs(exp(-out.weight)); % function 1

newAMat = repmat(newA', size(data, 1), 1);

for c = 1:2
    [ruleDataCov{c} ruleDataSupp{c}] = basRule(data, label, c);
    p0F = p0FPN(c);
    p0Fmat = p0F.*ones(n, dim);
    pgMat = (ruleDataCov{c}.*(ruleDataSupp{c}>p0FPN(c))).^newAMat.*(ruleDataSupp{c}-p0Fmat);
    linear = zeros(dim, 1);
    for i = 1:dim
        uni = unique(pgMat(:, i));
        quali = uni(find(uni)); % qualified feature values
        suppWeight{c}(i, 1) = sum(quali);
        clear uni quali
    end
end

for c = 1:2 % two class problem
    feaWeight = suppWeight{c};
    feaSDCheck = zeros(1, dim); % check usage of features
    idxCount = zeros(n, 1);
    iter = 0;
    class = c; % studied class
    [depen temp] = dataInfoQP(data, label, class, method);

    stopFlag = 0;
    exhaustiveCheck = 1;

    if strcmp(method, 'CS')
        H = -depen.*alpha;  
    else
        H = depen.*alpha; 
    end
    
    f = -feaWeight.*(1-alpha);
    
    feaCmp = ones(1, dim);
    options.LargeScale = 'off';
    classIdx = find(label == class);
    nClass = length(classIdx);
    fold = 0;
    fnew = 0;
    fval = 0;
    while stopFlag == 0  
        iter = iter + 1;
        %% quadratic programming, minimization approach
%         A = [];
%         b = [];
        A = ones(1, dim);
        b = 0.1; % sparsity control, Lasso like, inequality constraint
        
        Aeq = [];
        beq = [];
        lb = zeros(dim, 1);

        [x,fval,exitflag,output,lambda] = quadprog(H,f,A,b,Aeq,beq,lb,[],[],options);
        objValue(iter+1) = fval;
        
        if exitflag ~= 1
            disp('an optimal solution is not found!!!!!!!!!!!!!!!!')
        end
        x(find(x < tol), 1) = 0; % round off x
        [sortX indexX] = sort(x, 'descend');
        
        if length(find(x)) > maxWidth
            % if too many features, then select top ones
            candFeaIdx = indexX(1:maxWidth);
        else
            candFeaIdx = find(x);
        end
        [candFeaIdx] = sort(candFeaIdx, 'ascend');

        %%% check if exhaustive search is necessary, if done before!
        featuresUsed{iter} = candFeaIdx;
        if iter > 1
            [idx exhaustiveCheck] = sdCheck(featuresUsed, sampleCover);
        end
        
        if exhaustiveCheck == 1 || iter == 1
            [sd idx] = SDexhauSA(data(:, candFeaIdx), label, class, thres, out.weight(candFeaIdx));    
            %% correct the feature index by ture (ind)
            if ~strcmp(sd{1, 1}, 'empty')
                for k = 1:size(sd{1, 1}.feature, 1)
                    features = sd{1, 1}.feature{k, 1};
                    sd{1, 1}.feature{k, 1} = candFeaIdx(features);
                end     
            end
        end
        
        sdRules{iter, 1} = sd;
        sampleCover{iter, 1} = idx;
        %% reweight linear term
        idxCount = idxCount + idx;
%         weight = 1./(idxCount+1); % multiplicative descreasing
%         weight = 0.5.^(idxCount);

        weight = exp(-idxCount); % exponential descreasing
        feaWeight = covSuppWC(data, label, class, method, weight, newAMat);
        f = -feaWeight.*(1-alpha);
        
        % stop condition check
        if iter == 30 || abs((objValue(end)-objValue(end-1))/objValue(end)) < tol ...
                || sum(abs(f)) == 0
            stopFlag = 1;
        end
        clear idx
    end
    %% reorganize the sd rules
    SDRule = ruleOrgan(sdRules);
    SDRules{1, c} = SDRule;
    SDRules{1, c}.fval = objValue;
    clear SDRule sdRules sd
end








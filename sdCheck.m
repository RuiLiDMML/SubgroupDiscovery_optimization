function [idx check] = sdCheck(features, sample)
% check if exhaustive sd search necessary
n = length(features);

compared = features{end}; % the last one is the candidate features ought to be checked
idx = 0;
check = 1; % check is necessary
for i = 1:n-1
    tar = features{i};
    if length(tar)==length(compared)
        res = sum(abs(compared-tar));
        if res == 0 % found checked
            check = 0;
            idx = sample{i};
            break
        end
    end
end


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
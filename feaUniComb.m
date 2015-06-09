function combs = feaUniComb(F, n)
%%% unique feature combinations
if n > 20
    error('cannot hand feature more than 10')
end

switch n
    case 2
        combs = allcomb(F{1,1}, F{2,1});
    case 3
        combs = allcomb(F{1,1},F{2,1},F{3,1});
    case 4
        combs = allcomb(F{1,1},F{2,1},F{3,1},F{4,1});
    case 5
        combs = allcomb(F{1,1},F{2,1},F{3,1},F{4,1},F{5,1});
    case 6
        combs = allcomb(F{1,1},F{2,1},F{3,1},F{4,1},F{5,1},F{6,1});
    case 7
        combs = allcomb(F{1,1},F{2,1},F{3,1},F{4,1},F{5,1},F{6,1},F{7,1});
    case 8
        combs = allcomb(F{1,1},F{2,1},F{3,1},F{4,1},F{5,1},F{6,1},F{7,1},F{8,1});
    case 9
        combs = allcomb(F{1,1},F{2,1},F{3,1},F{4,1},F{5,1},F{6,1},F{7,1},F{8,1},F{9,1});
    case 10
        combs = allcomb(F{1,1},F{2,1},F{3,1},F{4,1},F{5,1},F{6,1},F{7,1},F{8,1},F{9,1},F{10,1});
    case 11
        combs = allcomb(F{1,1},F{2,1},F{3,1},F{4,1},F{5,1},F{6,1},F{7,1},F{8,1},F{9,1},F{10,1},F{11,1});
    case 12
        combs = allcomb(F{1,1},F{2,1},F{3,1},F{4,1},F{5,1},F{6,1},F{7,1},F{8,1},F{9,1},F{10,1},F{11,1}...
            ,F{12,1});   
    case 13
        combs = allcomb(F{1,1},F{2,1},F{3,1},F{4,1},F{5,1},F{6,1},F{7,1},F{8,1},F{9,1},F{10,1},F{11,1}...
            ,F{12,1},F{13,1});  
    case 14
        combs = allcomb(F{1,1},F{2,1},F{3,1},F{4,1},F{5,1},F{6,1},F{7,1},F{8,1},F{9,1},F{10,1},F{11,1}...
            ,F{12,1},F{13,1},F{14,1});  
    case 15
        combs = allcomb(F{1,1},F{2,1},F{3,1},F{4,1},F{5,1},F{6,1},F{7,1},F{8,1},F{9,1},F{10,1},F{11,1}...
            ,F{12,1},F{13,1},F{14,1},F{15,1});          
end

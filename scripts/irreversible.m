function [bModelIrrev] = irreversible(bModel)
%
%
lb = bModel.lb;ub = bModel.ub;
[nm nr] = size(bModel.S);
%
bModel.lb = -200*ones(nr,1);bModel.ub = 200*ones(nr,1);
bModel.rev = ones(nr,1);
%
[bModelIrrev x y z] = convertToIrreversible(bModel);
%
bModelIrev.lb = zeros(2*nr,1);bModelIrrev.ub = zeros(2*nr,1);
for i=1:nr
    %
    if lb(i) < 0 && ub(i) <=0
        bModelIrrev.ub(i+nr)   = -lb(i);
        bModelIrrev.lb(i+nr)   = -ub(i); 
    end
    %
    if lb(i) == 0 && ub(i) > 0
        bModelIrrev.ub(i) = ub(i);
    end
    %
    if lb(i) > 0 && ub(i) > 0
        bModelIrrev.lb(i) = lb(i);
        bModelIrrev.ub(i) = ub(i);
    end
    %
    if lb(i) < 0 && ub(i) > 0
        bModelIrrev.ub(i) = ub(i);
        bModelIrrev.ub(i+nr)   = -lb(i);
    end
end

end


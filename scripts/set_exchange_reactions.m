function [outmodel] = set_exchange_reactions(inmodel)

model = inmodel;

%Set unit growth rate
model.lb(3741) = 1;model.ub(3741) = 1+1e-5;
model.lb(3744) = 0; %No maintenance ATP constraint
model.ub(3744) = 200; %No maintenance ATP constraint

%Free all exchanges
onrxns = setdiff(1:3744,find(model.lb == 0 & model.ub == 0));
all_ex_rxns = find(~cellfun(@isempty,strfind(model.rxns,'EX_')));
ex_rxns = intersect(onrxns,all_ex_rxns);

%Free all exhanges
model.lb(ex_rxns) = -200; model.ub(ex_rxns) = 200;

model.lb(ex_rxns) = 0;
%Allow these reactions to be exported
in_rxns = [1353 1354 1385 1386 1408 1409 1412 1465 1472 1495 1524 1526];
model.lb(in_rxns) = -200;

outmodel = model;

end


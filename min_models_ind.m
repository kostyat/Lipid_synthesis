clear
rng(xSEED)
%%
load results/min_model_params_new
load data/zeromodel.mat
model0 = model;
load data/small_model
%%
flux_all = samp(nutN,:);
% Shut down all exchanges
all_exs = find(strncmp('EX_',model0.rxns,3));
model0.lb(all_exs) = 0;model0.ub(all_exs) = 0;
%% Allow favorite exchanges
model0.lb(exInds) = -1000;model0.ub(exInds) = 1000;
% nutrients exchange reactions:
nutI = [1264 1277 1279 1281 1385 1386 1387 1389 1427 1437 1439];
nutI = [nutI 1449 1457 1490 1499 1524 1551 1560 1567 1576];

%% For every cancer type
rxnMatInd = {};
for ctype = xN1:xN2
    mu = gRates(ctype);
    alph = alpha0s(ctype);
    gModel = convert_to_gurobi(model0);
    gModel.lb(3741) = mu;gModel.ub(3741) = mu+1e-10;
    gModel.obj = zeros(3744,1);gModel.obj(nutI) = -2*alph*flux_all(:,ctype);
    q = sparse(3744,1);q(nutI) = 1;Q = diag(q);gModel.Q = Q;
    f = gurobi(gModel,params);
    fm = alph*flux_all(:,ctype);tt = (abs((f.x(nutI)-fm)./fm));
    fsl = f.x(nutI);

    %%
    gModel = convert_to_gurobi(model0);
    gModel.lb(3741) = mu;gModel.ub(3741) = mu+1e-10;
    gModel.lb(nutI) = fsl-1e-10;gModel.ub(nutI) = fsl-1e-10;
    %%
    noremove  = [nutI 3741 3744];
    rxnMatrix = [];
    nModels = 1000; % Number of small models
    for iM=1:nModels
        rxnx = 1:3744;rxnx = setdiff(rxnx,noremove);
        tModel = gModel;tModel.obj = zeros(3744,1);
        nonr = [];
        for iter=1:10000
            fModel = tModel;
            rxnInd = randsample(rxnx,1);
            if ~ismember(rxnInd,nonr)     
                fModel.lb(rxnInd) = 0;fModel.ub(rxnInd) = 0;
                sol = gurobi(fModel,params);
                if strcmp(sol.status,'OPTIMAL')
                    tModel = fModel;rxnx(rxnx==rxnInd) = [];
                else
                    nonr = [nonr rxnInd];
                end
            end
            %
            if mod(iter,2500) == 0
                [iM length(nonr)]
            end
        end
        tmp = sparse(3744,1);tmp(nonr) = 1;tmp(noremove) = 1;
        rxnMatrix = [rxnMatrix tmp];
    end
    rxnMatInd(ctype) = {rxnMatrix};
    %save xFILENAME rxnMatInd
    save(["xFILENAME"+"_"+ctype], 'rxnMatrix')
end


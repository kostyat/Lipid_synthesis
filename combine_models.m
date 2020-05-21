fold = "run_files/";
%%
for run=0:11
    for i=1:5
        fname = "rxnMatrixInd_run"+run+"_"+(i+run*5)+".mat";
        if(isfile(fold+fname))
            load(fold+fname)
            min_mod_conv = [];
            for(j=1:1000)
                min_mod_conv = [min_mod_conv find(rxnMatrix(:,j))'];
            end
            models_probs(i+run*5).tab = tabulate(min_mod_conv);
        else
            models_probs(i+run*5).tab = [];
        end
    end
end

%%
load results/min_model_params_new

load data/zeromodel
model0 = model;
load small_model

% reactions that cannot be removed:
remove_idx = 1:length(model.rxns);
remove_idx([biomass_rxn,main_atp_rxn,modelNut]) = [];

%%
flux_all = samp(nutN,:);
% Shut down all exchanges
all_exs = find(strncmp('EX_',model0.rxns,3));
model0.lb(all_exs) = 0;model0.ub(all_exs) = 0;
%% Allow favorite exchanges
model0.lb(exInds) = -1000;model0.ub(exInds) = 1000;
nutI = modelNut;

%% Loop over cancer types
fsl_mat = [];
for ctype = 1:length(models_probs)
    if isempty(models_probs(ctype).tab)
        fsl_mat = [fsl_mat zeros(length(nutI),1)];
        continue
    end
    mu = gRates(ctype);
    alph = alpha0s(ctype);
    %% constrain exchange reactions using the alph that was already optimized
    gModel = convert_to_gurobi(model0);
    gModel.lb(3741) = mu;gModel.ub(3741) = mu+1e-10;
    gModel.obj = zeros(3744,1);gModel.obj(nutI) = -2*alph*flux_all(:,ctype);
    q = sparse(3744,1);q(nutI) = 1;Q = diag(q);gModel.Q = Q;
    f = gurobi(gModel,params);
    fm = alph*flux_all(:,ctype);tt = (abs((f.x(nutI)-fm)./fm));
    fsl = f.x(nutI);
    fsl_list(ctype).fsl = fsl; fsl_list(ctype).tt = tt;
    fsl_list(ctype).tt_mean = mean(tt); fsl_list(ctype).fm = fm;
    fsl_mat = [fsl_mat fsl];
    %% set lb and ub of exchange reactions to the optimized fluxes
    gModel = convert_to_gurobi(model0);
    gModel.lb(3741) = mu;gModel.ub(3741) = mu+1e-10;
    gModel.lb(nutI) = fsl-1e-10;gModel.ub(nutI) = fsl-1e-10;
    %% find the minimal model
    rxn_count = models_probs(ctype).tab;
    rxn_count_sorted = sort(unique(rxn_count(:,2)), 'descend');
    %%
    for i=1:length(rxn_count_sorted)
        %%
        cutoff = rxn_count_sorted(i);
        which_remove = find(rxn_count(:,2) < cutoff);
        rxns_remove = rxn_count(which_remove,1);
        rxnx = intersect(rxns_remove, remove_idx);
        tModel = gModel;tModel.obj = zeros(3744,1);
        cutoff
        tModel.lb(rxnx) = 0;tModel.ub(rxnx) = 0;
        sol = gurobi(tModel,params);
        %%
        if(strcmp(sol.status, 'OPTIMAL'))
            minmodel = tModel;
            sol_minmodel = sol;
            cutoff_minmodel=cutoff;
            break
        end
    end
    min_models(ctype).model = minmodel;
end

%% get model sizes
model_sizes = [];
for(ctype=1:length(min_models))
    if isempty(min_models(ctype).model)
        model_sizes = [model_sizes 0];
    else
       model_sizes = [model_sizes sum(min_models(ctype).model.lb~=0 | min_models(ctype).model.ub~=0)]; 
    end
end

%% save
save results/FinalModels.mat min_models model_sizes

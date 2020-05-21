clear
%
addpath('scripts/')
load data/zeromodel.mat

params=struct();
params.outputflag=0;

load FinalModels

%%
costs = [];
for i=1:length(min_models)
    if isempty(min_models(i).model)
        costs = [costs NaN];
        nad_rxns(i).cost = NaN;
        nad_rxns(i).flux = NaN;
        continue
    end
    %%
    model = convert_to_cobra(min_models(i).model, model);

    % Set NAD reactions
    [model] = set_nad_reactions(model);

    % Set up model with only a few input reactions
    ac_model = set_exchange_reactions(model);

    %% Don't let any other carbon sources except glucose and acetate
    ac_model.lb([1524 1386]) = 0;

    % shut down GLUDym, ALR, PEPCK, ACACT1r, ALCD21_L, and LCADi
    badrxs = [2072 383 2942 167 362 2375];
    ac_model.lb(badrxs) = 0;ac_model.ub(badrxs) = 0; 

    % Add accoa->lipids reactions
    ac_model = add_accoa_reaction(ac_model);

    %Allow acetate
    ac_model.lb(1238) = -200;ac_model.ub(1238) = 200;
    ac_model.lb(253) = -200;ac_model.ub(253) = 200;

    %% Convert to irreversible model
    imodel = irreversible(ac_model);
    [iM,iN] = size(imodel.S);

    %Conver to gurobi model
    gModel = convert_to_gurobi(imodel);

    %% Minimize NAD cost
    nad_c = 2043;nad_m = 2045;
    nad_c_rxns_inds = find(gModel.A(nad_c,:));
    nad_c_dirs = gModel.A(nad_c,find(gModel.A(nad_c,:)));
    nad_m_rxns_inds = find(gModel.A(nad_m,:));
    nad_m_dirs = gModel.A(nad_m,find(gModel.A(nad_m,:)));
    nad_c_prod_inds = nad_c_rxns_inds(find(nad_c_dirs > 0));
    nad_m_prod_inds = nad_m_rxns_inds(find(nad_m_dirs > 0));
    gModel.obj([nad_c_prod_inds nad_m_prod_inds]) = 1;

    %% Solve initial model
    sol1 = gurobi(gModel, params);

    %% Update new model with this as a constraint and minimize 1 norm
    gModel.A(end+1,:) = gModel.obj; gModel.sense(end+1) = '<'; 
    gModel.rhs(end+1) = sol1.objval + 0.0001;
    gModel.obj = ones(iN,1);

    %Solve new model and get realistic flux distribution
    sol2 = gurobi(gModel, params);
    [iM,iN] = size(gModel.A);iN = iN/2;
    flx = sol2.x(1:iN) - sol2.x(iN+1:end);

    %Lipid NAD+ cost from glutamine
    nad_accoa_ac = sol1.objval;
    costs = [costs nad_accoa_ac];
    % 
    nadcon1 = find(ac_model.S(2043,:).*flx' < 0);
    nadcon2 = find(ac_model.S(2045,:).*flx' < 0);
    nadcon = [nadcon1 nadcon2];

    flxn = flx(nadcon);[a b] = sort(flxn,'descend');
    nad_rxns(i).cost = model.rxnNames(nadcon(b));

    nad_rxns(i).flux = flx(nadcon(b));
end

%% Save
save results/acetate_accoa_costs.mat costs nad_rxns



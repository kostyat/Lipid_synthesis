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
    ser_model = set_exchange_reactions(model);

    %% Don't let any other carbon sources except glucose and glutamine
    ser_model.lb([1524]) = 0;

    % shut down GLUDym, ALR, and PEPCK, ACACT1r, ALCD21_L, and LCADi
    badrxs = [2072 383 2942 167 362 2375];
    ser_model.lb(badrxs) = 0;ser_model.ub(badrxs) = 0;
    
    % Constrain aspartate transaminase
    %ser_model.lb([501])  = -0.5104;
    %ser_model.lb([501])  = -0.7470;
    ser_model.lb([501])  = -0.8235;

    %Change biomass reaction of lipid model
    biomass_mets = find(model.S(:,3741));
    ser_mets = [2436];
    nonser_mets = setdiff(biomass_mets,ser_mets);
    ser_model.S(nonser_mets,3741) = 0;

    %% Convert to irreversible model
    imodel = irreversible(ser_model);
    [iM,iN] = size(imodel.S);

    %Convert to gurobi model
    gModel = convert_to_gurobi(imodel);

    %%Minimize total NAD+ production in cytosol and mitochondria
    nad_c = 2043;nad_m = 2045;
    nad_c_rxns_inds = find(gModel.A(nad_c,:));
    nad_c_dirs = gModel.A(nad_c,find(gModel.A(nad_c,:)));
    nad_m_rxns_inds = find(gModel.A(nad_m,:));
    nad_m_dirs = gModel.A(nad_m,find(gModel.A(nad_m,:)));
    nad_c_prod_inds = nad_c_rxns_inds(find(nad_c_dirs > 0));
    nad_m_prod_inds = nad_m_rxns_inds(find(nad_m_dirs > 0));
    gModel.obj([nad_c_prod_inds nad_m_prod_inds]) = 1;

    %Solve initial model
    sol1 = gurobi(gModel, params);

    %% Update new model with this as a constraint and minimize 1 norm
    gModel.A(end+1,:) = gModel.obj; gModel.sense(end+1) = '<';
    gModel.rhs(end+1) = sol1.objval + .0001;gModel.obj = ones(iN,1);

    %Solve new model and get realistic flux distribution
    sol2 = gurobi(gModel, params);
    [iM iN] = size(gModel.A);iN = iN/2;
    flx = sol2.x(1:iN) - sol2.x(iN+1:end);

    %Lipid NAD+ cost from serine
    nad_ser = sol1.objval;
    costs = [costs nad_ser];
    
    nadcon1 = find(ser_model.S(2043,:).*flx' < 0);
    nadcon2 = find(ser_model.S(2045,:).*flx' < 0);
    nadcon = [nadcon1 nadcon2];

    flxn = flx(nadcon);[a b] = sort(flxn,'descend');
    nad_rxns(i).cost = model.rxnNames(nadcon(b));
    
    nad_rxns(i).flux = flx(nadcon(b));
end

%% Save
save results/serine_costs.mat costs nad_rxns






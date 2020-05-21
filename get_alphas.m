%initCobraToolbox;
changeCobraSolver('gurobi','all');
load('data/data.mat')

%% First, set all exchange fluxes to zero:
all_ex_fluxes = find(not(cellfun('isempty',strfind(model.rxns, "EX_"))));
model.ub(all_ex_fluxes) = 0;
model.lb(all_ex_fluxes) = 0;
%Then, set the bounds on fluxes of exchange reactions that we want to keep to +-1000:
model.ub(exInds) = 1000;
model.lb(exInds) = -1000;

%% find the biomass reaction
biomass_rxn = find(not(cellfun('isempty',strfind(model.rxns, "biomass"))));
%set the ATP maintenance reaction
main_atp_rxn = find(not(cellfun('isempty',strfind(model.rxns, "maintenance_ATP"))));
model.ub(main_atp_rxn) = 3.84+1e-10;
model.lb(main_atp_rxn) = 3.84-1e-10;
%Set the objective function as the biomass reaction:
model.c = model.c.*0;
model.c(biomass_rxn) = 1;

FBAsolution = optimizeCbModel(model,'max');
%This is supposed to be infeasible, so make a new model that's the same but isn't constrained by the fluxed in the dataset:
%because we want them to be constrained by something proportional to them (w/ alpha)

%% %%%%%
%% Switch to Gurobi solver
model_gur.varnames = model.rxns;
%set up the constraints matrix
model_gur.A = model.S;
% set constraints on the right hand side
model_gur.rhs = zeros(size(model.S,1),1);
model_gur.sense = '=';
%set up the quadratic objective function:
model_gur.Q = zeros(length(model.rxns));
Q_diag = diag(model_gur.Q);
Q_diag(modelNut) = 1;
model_gur.Q = sparse(model_gur.Q + diag(Q_diag));

%set up the linear objective function for various values of alpha
model_gur.obj = Q_diag;
params.outputflag = 0;
alphas = 0:0.00001:0.005;

%% Cycle over the 60 cell lines

for j = 1:length(gRates)
    %Set the growth rate for the cell line
    model.ub(biomass_rxn) = gRates(j)+1e-10;
    model.lb(biomass_rxn) = gRates(j)-1e-10;
          
    %add the lb and ub rows:
    model_gur.lb = model.lb;
    model_gur.ub = model.ub;

    %set up exchange flux for each cancer type
    flux_137 = samp(:,j);
    flux_20 = flux_137(nutN);

    objvals = [];
    for i = 1:length(alphas)
        model_gur.obj(modelNut) = -2*alphas(i)*flux_20;
        results  = gurobi(model_gur, params);
        objvals(i) = results.objval;
        sol = results.x;
        dX = abs(sol(modelNut)-(alphas(i)*(flux_20)));
        dX = dX./flux_20;
        erx(i)=mean(abs(dX));
    end

    alpha0 = alphas(find(erx == min(erx)));
    alpha0s(j) = alpha0;
    j
end

%%
save("results/min_model_params_new.mat");


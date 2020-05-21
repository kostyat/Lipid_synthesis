function [outmodel] = set_nad_reactions(inmodel)
model = inmodel;

%% Cystosolic NAD+

%Reactions/metabolites of interest
nad_c = 2043;
nad_c_rxns_inds = find(model.S(nad_c,:));
nad_c_dirs = model.S(nad_c,find(model.S(nad_c,:)));

%Reactions to leave on
desat_inds = find(~cellfun(@isempty,strfind(model.rxns,'DESAT')));
mdh_ind = 2511; 
glydh_ind = 1861;
ldh_ind = 2386;
onrxns_c = [mdh_ind; glydh_ind; ldh_ind; desat_inds];

%Turn of all sources of NAD+ in the cytosol except MA and GP shuttles
for i = 1:length(nad_c_rxns_inds)
    
    rxn_ind = nad_c_rxns_inds(i);
    dir = nad_c_dirs(i);
    
    %If reaction consumes NAD+
    if dir < 0
        model.lb(rxn_ind) = 0;
    end
    
    %If reaction produces NAD+
    if dir > 0
        model.ub(rxn_ind) = 0;
    end
    
    %If reaction is allowable (free both directions)
    if ismember(rxn_ind,onrxns_c)
        model.lb(rxn_ind) = -200;
        model.ub(rxn_ind) = 200;
    end
   
end

%% Mitochondrial NAD+

%Reactions/metabolites of interest
nad_m = 2045;
c1_ind = 2664;
nad_m_rxns_inds = find(model.S(nad_m,:));
nad_m_dirs = model.S(nad_m,find(model.S(nad_m,:)));

%Reactions to leave on
onrxns_m = [c1_ind];


%Turn of all sources of NAD+ in the cytosol except MA and GP shuttles
for i = 1:length(nad_m_rxns_inds)
    
    rxn_ind = nad_m_rxns_inds(i);
    dir = nad_m_dirs(i);
    
    %If reaction consumes NAD+
    if dir < 0
        model.lb(rxn_ind) = 0;
    end
    
    %If reaction produces NAD+
    if dir > 0
        model.ub(rxn_ind) = 0;
    end
    
    %If reaction is allowable (free both directions)
    if ismember(rxn_ind,onrxns_m)
        model.lb(rxn_ind) = 0;
        model.ub(rxn_ind) = 200;
    end
    
  
end

outmodel = model;

end


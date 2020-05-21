function [outmodel] = add_accoa_reaction(inmodel)

glc_model = inmodel;
% Add a accoa -> "lipid" + coa reaction
[iM iN] = size(glc_model.S);
glc_model.S(393,iN+1) = -1;
glc_model.S(780,iN+1) = 1;
glc_model.S(iM+1,iN+1) = 1;
glc_model.lb(iN+1) = -200;glc_model.ub(iN+1) = 200;
glc_model.mets{iM+1} = 'Lipids';
glc_model.rxns{iN+1} = 'Lipid synthesis';
glc_model.c(iN+1) = 0;glc_model.b(iM+1) = 0;
%
lipid_mets = [733 751 2191 2200 2214 2231 2333 2463];
nAcUnit = [18 32 16 16 16 16 16 18]; % Number of accoa units
accoaUnits = sum(full(glc_model.S(lipid_mets,3741)).*nAcUnit');
glc_model.S(:,3741) = 0;glc_model.S(iM+1,3741) = accoaUnits;

outmodel = glc_model;

end


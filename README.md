# Lipid_synthesis
FBA model of lipid synthesis in hypoxia.

First run `get_alphas.m` to obtain the scaling parameters alpha for each tumor type.

Then run `run-min-model_ind.sh` to create 1000 minimal models in the `run_files/` folder.

Then run `combine_models.m` to obtain a consensus minimal model (saved as `resuls/FinalModels.m`).

Then, to get the NAD+ consumption cost for each nutrient source, run each file in the `cost-calcs/` folder.

All code written by Brian W. Ji, Purushottam D. Dixit, and Konstantine Shapiro-Tchourine.

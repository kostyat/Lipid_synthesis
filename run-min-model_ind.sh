#!/usr/bin/

ar1=($(seq 1 5 60))
ar2=($(seq 5 5 60))

mkdir run_files

for i in $(seq 0 1 $((${#ar1[@]}-1)))
do
  cp min_models_ind.m run_files/min_models_ind_run$i.m
  perl -pi -e "s{xFILENAME}{run_files/rxnMatrixInd_run$i}g" run_files/min_models_ind_run$i.m
  perl -pi -e "s{xSEED}{42}g" run_files/min_models_ind_run$i.m
  perl -pi -e "s{xN1}{${ar1[i]}}g" run_files/min_models_ind_run$i.m
  perl -pi -e "s{xN2}{${ar2[i]}}g" run_files/min_models_ind_run$i.m
  /usr/local/MATLAB/R2019a/bin/matlab -nodisplay -nosplash -nodesktop -r "run('run_files/min_models_ind_run$i.m');exit;" > run_files/out$i.txt&
done


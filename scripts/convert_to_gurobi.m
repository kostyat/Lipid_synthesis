function gModel = convert_to_gurobi(model)
[M,N] = size(model.S);
gModel.A = model.S;
gModel.lb = model.lb;
gModel.ub = model.ub;
gModel.rhs = model.b;
gModel.sense(1:M)  = '=';
gModel.obj = model.c;

gModel.obj = zeros(N,1);
gModel.modelsense = 'min';
gModel.sense(1:M) = '=';

end

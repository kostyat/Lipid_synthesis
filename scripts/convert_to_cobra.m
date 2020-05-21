function model = convert_to_cobra(gModel, zero_model)
[M,N] = size(gModel.A);
model = zero_model;
model.S = gModel.A;
model.lb = gModel.lb;
model.ub = gModel.ub;
model.b = gModel.rhs;
model.sense(1:M)  = '=';
model.c = gModel.obj;

model.c = zeros(N,1);
model.modelsense = 'min';
model.csense(1:M) = '=';

end

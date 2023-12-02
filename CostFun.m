function Cost = CostFun(par)
global model
N1 = model.N1;
N2 = model.N2;

N = Num(par,N1);

s = N-N2;
Cost = norm(s);
function N = Num(par,N1)

P = gtm_optimization(par); 
% P = gtm(par); 
N = P*N1;


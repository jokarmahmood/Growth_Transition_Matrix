clc; clear;
%% This code finds matrix G that satisfy N_k+1=G*N_k using proposed method:

D1 = load('R3.mat');
N1 = D1.R3'*10^9;

DD = load('X3.mat');
P = DD.X3;

N2 = P*N1;
length_class = length(N1); 
step = 1;
degree = 2;
x = 5.5:step:7.5;
nx = length(x);
xnew = 5:step:8;
nxnew = length(xnew);

global model
model.length_class = length_class;
model.N1 = N1;
model.N2 = N2;
model.x = x;
model.step = step;

% Construct B-Spline coeeficients of degree p 
% index = find(N2);
model.degree = degree;
model.knot = knots(xnew,nxnew-1,degree);
% par = [c0 ... ci]
par = ones(1,nxnew);                  % vector of model parameters
% Optimization 
lb = zeros(1,nxnew);
ub = max(N2/10^9)*ones(1,nxnew);
disp('Optimization started to estimate the model parameters...')
fprintf(1,'%d parameters should be estimated\n', nxnew);
%%%%%% method PSO
options = optimoptions('particleswarm','SwarmSize',50,'HybridFcn',@fmincon);
[Epar,fval,exitflag,output] = particleswarm(@CostFun,length(par),lb,ub,options);
% result
fprintf(1,'fval: %e\n', fval);
for k = 1:nxnew
    fprintf('c%d: %d\n', k-1,Epar(k));
end
G = gtm_optimization(Epar);  % Estimated growth transition matrix
%G = gtm(par); 

M = G*N1;              % Estimated N_t+1

figure(1)
plot(x,N2,'b','LineWidth',2)
hold on
plot(x,M,'r--','LineWidth',2)
legend('True','Simulation','Location','northwest')
xlabel('length')
h = ylabel('$N_{t+1}$','fontweight','bold','fontsize',14);
set(h,'Interpreter','latex')
%title('Proposed method')

figure(2)
subplot(1,2,1)
plot(x,N1,'LineWidth',2)
hold on
plot(x,N2,'LineWidth',2)
xlabel('length')
ylabel('number')
title('Transition')
h = legend('$N_t$ ','$N_{t+1}$ ');
set(h,'Interpreter','latex')
subplot(1,2,2)
imagesc(P)
title('True growth transition matrix')

median_e = sum(sum(abs(P-G)))/(nx*(nx+1)/2);
fprintf(1,'Matrix median difference error: %d\n', median_e*100);

figure(3)
subplot(1,2,1), imagesc(P)
title('True P')
subplot(1,2,2), imagesc(G)
title('Proposed  approach')

%%
%close all;
function P = gtm_optimization(par)
global model
x = model.x;
step = model.step;
nx = length(x);
coeff = par;
P = zeros(nx,nx);
I = zeros(1,nx);
f = @ (y) BSinterp2(y,coeff);
for j = 1:nx   
    I(j) = simpson(f,x(j)-step/2,x(end)+step/2); % Simpson
end
for i = 1:nx % index(1):index(end) 
    for j = 1:i    
         P(i,j) = simpson(f,x(i)-step/2,x(i)+step/2)/I(j); % Simpson
    end
end
for j = 1:nx % index(1):index(end)    
    P(:,j) = P(:,j)/sum(P(:,j)); 
end
function P = gtm(par)
global model
x = model.x;
index = model.index;
nx = length(x);
P = zeros(nx,nx);
coeff = zeros(nx,1);
coeff(index) = par;
for i = index(1):index(end)
    for j = 1:i    
        a = BSinterp2(x(i)-0.25,coeff);
        a(a<0) = 0;
        b = BSinterp2(x(i),coeff);
        b(b<0) = 0;
        c = BSinterp2(x(i)+0.25,coeff);
        c(c<0) = 0;
        P(i,j) = (1/18)*(a+4*b+c); % Simpson
    end
end

for j = 1:nx
    if sum(P(:,j))>0
        P(:,j) = P(:,j)/sum(P(:,j));
    end 
end
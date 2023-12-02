function U = knots(t,n,p)

m = n+p+1;
U = zeros(1,m+1);
U(1:p+1) = t(1);
for j = 2:n-p+1
    U(j+p) = (1/p)*sum(t(j:j+p-1));
end
U(m-p+1:m+1) = t(end);
    

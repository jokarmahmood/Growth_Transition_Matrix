function y = BSinterp2(x,c)
% This function interpolates using B-splines of degree p 
% to find the function values at vector x

global model
p = model.degree;
n = length(x);
if n == 1
    phi = B_spline(p,x);
    y = c*phi;
else
    y = zeros(n,1);
    for i = 1:n
        phi = B_spline(p,x(i));
        y(i) = c*phi;
    end
end
function N = B_spline(p,u)
% This function computes the coefficients of B-spline curves of degree p

global model
U = model.knot;
m = length(U)-1;
n = m-p-1;
N = zeros(n+1,1);
if u == U(1) 
    N(1) = 1;
elseif u == U(end)
    N(n+1) = 1;
else
    for j = 1:m
        if u>=U(j) && u<U(j+1)
            k = j;
            break
        end
        if u>=U(j) && U(j+1)==U(end)
            k = j;
            break
        end
    end
    N(k) = 1;
    for d = 1:p
        N(k-d) = (U(k+1)-u)/(U(k+1)-U(k-d+1))*N(k-d+1);
        for i = k-d+1:k-1
             N(i) = (u-U(i))/(U(i+d)-U(i))*N(i)+(U(i+d+1)-u)/(U(i+d+1)-U(i+1))*N(i+1);
        end
        N(k) = (u-U(k))/(U(k+d)-U(k))*N(k);
    end
end
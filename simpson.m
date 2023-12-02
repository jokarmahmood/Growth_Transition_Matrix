function I = simpson(f,a,b) 
x = linspace(a,b,100);
h = (b-a)/100;
I = h/3*(f(x(1))+2*sum(f(x(3:2:end-2)))+4*sum(f(x(2:2:end)))+f(x(end)));
% n = length(x);
% 
% a = f(x(1)); a(a<0) = 0;
% S = a;
% for i = 2:n-1
%     a = f(x(i)); a(a<0) = 0;
%     if rem(i,2) == 0  
%         S = S + 4*a;
%     else
%         S = S + 2*a;
%     end
% end
% a = f(x(n)); a(a<0) = 0;
% S = S + a;  
% I = h/3*S;
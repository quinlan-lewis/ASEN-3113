function [Temp_func] = T_x_t(alpha,L,H,T_0,x,t)
Tmp = zeros(length(t),1);

for j = 1:length(t)
    sum = 0;
    for i = 1:10
        b_n = (-8*H*L/(pi^2))*((-1)^(i+1)/(2*i-1)^2);
        lambda = (2*i-1)*pi/(2*L);
        
        sum = sum + b_n*sin(lambda*x)*exp(-lambda^2*alpha*t(j));
    end
    Tmp(j) = T_0 + H*x + sum;
end

Temp_func = Tmp;

end


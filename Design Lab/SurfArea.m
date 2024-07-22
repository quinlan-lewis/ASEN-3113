function [Nadir, Radiator, Others] = SurfArea(T,As,rho_s)
sigma = 5.67e-8;
alpha0 = 0.54;
alpha1 = 0.17;
eps0 = 0.75;
eps1 = 0.92;
rho_e = 33;

syms A_N A_R A_O

eq_N = T == ((A_N*As(1)*alpha0*(rho_e + rho_s) + 20)/(eps0*sigma*A_N))^(1/4);
Nadir = vpasolve(eq_N,A_N);

eq_R = T == ((A_R*As(2)*alpha1*(rho_s) + 20)/(eps1*sigma*A_R))^(1/4);
Radiator = vpasolve(eq_R,A_R);

Others = zeros(1,4);
for i = 3:6
eq_O = T == ((A_O*As(i)*alpha0*(rho_s) + 20)/(eps0*sigma*A_O))^(1/4);
Others(i-2) = vpasolve(eq_O,A_O);
end

end


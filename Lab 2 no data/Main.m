clear all; clc; close all;

exp_data_Al2 = load('Spring_data\Spring 2021 Data\Aluminum30V_290mA');
exp_data_Al1 = load('Spring_data\Spring 2021 Data\Aluminum25V_240mA');
exp_data_Brass = load('Spring_data\Spring 2021 Data\Brass_30V_285mA');
exp_data_Steel = load('Spring_data\Spring 2021 Data\Steel22V_203mA');

%below linear distance is assumed the same for all tests
L = 1.375 + .5*7 + 1; %[in]
alphas = [(4.82e-5); (4.05e-6); (3.56e-5)]; %[m^2/s]
alphas = 1550.*alphas; %[in^2/s]
d1 = 1.375; %[in]
Linear_dist = [0; d1; d1+(.5); d1+(.5*2); d1+(.5*3); d1+(.5*4); d1+(.5*5); d1+(.5*6); d1+(.5*7)]; %[in]

%interpolations
%AL1:
Temps_Al1 = exp_data_Al1(length(exp_data_Al1),2:10);
[T_0_Al1, H_Al1, Temps_out_Al1] = Lin_interp(Temps_Al1,Linear_dist);
figure(); hold on;
plot(Linear_dist,Temps_Al1);
plot(Linear_dist,Temps_out_Al1);
title('Alumium, V = 25'); xlabel('Linear Distance [in]'); ylabel('Temperature [deg C]')
legend('Experimental data','Fitted data','location','SouthEast')

%AL2:
Temps_Al2 = exp_data_Al2(length(exp_data_Al2),2:10);
[T_0_Al2, H_Al2, Temps_out_Al2] = Lin_interp(Temps_Al2,Linear_dist);
figure(); hold on;
plot(Linear_dist,Temps_Al2);
plot(Linear_dist,Temps_out_Al2);
title('Alumium, V = 30'); xlabel('Linear Distance [in]'); ylabel('Temperature [deg C]')
legend('Experimental data','Fitted data','location','SouthEast')

%Brass:
Temps_Brass = exp_data_Brass(length(exp_data_Brass),2:10);
[T_0_Brass, H_Brass, Temps_out_Brass] = Lin_interp(Temps_Brass,Linear_dist);
figure(); hold on;
plot(Linear_dist,Temps_Brass);
plot(Linear_dist,Temps_out_Brass);
title('Brass, V = 30'); xlabel('Linear Distance [in]'); ylabel('Temperature [deg C]')
legend('Experimental data','Fitted data','location','SouthEast')

%Steel:
Temps_Steel = exp_data_Steel(length(exp_data_Steel),2:10);
[T_0_Steel, H_Steel, Temps_out_Steel] = Lin_interp(Temps_Steel,Linear_dist);
figure(); hold on;
plot(Linear_dist,Temps_Steel);
plot(Linear_dist,Temps_out_Steel);
title('Steel, V = 22'); xlabel('Linear Distance [in]'); ylabel('Temperature [deg C]')
legend('Experimental data','Fitted data','location','SouthEast')

%Question 2
%AL1:
x = Linear_dist(2:end);
t_Al1 = exp_data_Al2(:,1);
Temp_func_Al1 = zeros(length(t_Al1),length(x));
for i = 1:length(x)
    tmp = T_x_t(alphas(1),L,H_Al1,T_0_Al1,x(i),t_Al1);
    for j = 1:length(tmp)
       Temp_func_Al1(j,i) = tmp(j,1); 
    end
end

figure(); hold on;
for i = 1:length(x)
    plot(t_Al1,Temp_func_Al1(:,i))
end

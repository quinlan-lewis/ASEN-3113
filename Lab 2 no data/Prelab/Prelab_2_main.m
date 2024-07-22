%% Prelab 2 - 3113: Thermodynamics and Heat Transfer
% Author: Quinlan Lewis
% Last Modified: 3/7/21

clear all; clc; close all;

%% Givens
alphas = [(4.82e-5); (4.05e-6); (3.56e-5)]; %[m^2/s]
alphas = 1550.*alphas; %[in^2/s]
L = 1.375 + .5*7 + 1; %[in]

%% Question 3: Interpolation for T_0
Th_readings = [18.53; 22.47; 26.87; 30.05; 35.87; 38.56; 41.50; 46.26]; %[deg C]
d1 = 1.375; %[in]
Linear_dist = [0; d1; d1+(.5); d1+(.5*2); d1+(.5*3); d1+(.5*4); d1+(.5*5); d1+(.5*6); d1+(.5*7)]; %[in]

%least square fitting for the temperature
coefs = polyfit(Linear_dist(2:end,:),Th_readings,1);
H = coefs(1);
fitted = polyval(coefs,Linear_dist);
figure()
hold on;
plot(Linear_dist(2:end,:), Th_readings)
plot(Linear_dist,fitted)
xlabel('Linear distance [in]'); ylabel('Temperature [Celcius]');
title('Distance vs Temperature')
legend('Given Data','Fitted Line','Location','SouthEast')

%initial temp occurs at 0 linear distance
T_0 = fitted(1);
fprintf('Initial Temperature: %f \nSlope of Curve: %f\n',T_0,H)

%% Question 5: Calculating temperature as a function of time and position
t1= 1; %[sec]
t2 = 1000; %[sec]
x = L - 1; %[in]
u1 = zeros(10:1);
u2 = zeros(10:1);

for j = 1:10
    sum1 = 0;
    sum2 = 0;
    for i = 1:j
        b_n = (-8*H*L/(pi^2))*((-1)^(i+1)/(2*i-1)^2);
        lambda = (2*i-1)*pi/(2*L);
        
        sum1 = sum1 + b_n*sin(lambda*x)*exp(-lambda^2*alphas(1)*t1);
        sum2 = sum2 + b_n*sin(lambda*x)*exp(-lambda^2*alphas(1)*t2);
    end
    u1(j) = T_0 + H*x + sum1;
    u2(j) = T_0 + H*x + sum2;
end

F_01 = alphas(1)*t1/L^2;
F_02 = alphas(1)*t2/L^2;
fprintf('F_0 for 1 sec: %f\nF_0 for 1000 sec: %f\n',F_01,F_02);

figure();
plot(1:10, u1);
xlabel('Fourier Coefficients'); ylabel('Temperature in Celcius');
title('Converging Temperature for t = 1s')
figure();
plot(1:10, u2);
xlabel('Fourier Coefficients'); ylabel('Temperature in Celcius');
title('Converging Temperature for t = 1000s')

%% Question 6: Sensitivity study on Diffusivity
temps = zeros(100:length(alphas));

for j = 1:length(alphas)
   for t =  1:5000
       b_n = (-8*H*L/(pi^2))*((-1)^(1+1)/(2*1-1)^2);
       lambda = (2*1-1)*pi/(2*L);
       
       temps(t,j) = T_0 + H*x + b_n*sin(lambda*x)*exp(-lambda^2*alphas(j)*t);
   end
end

figure(); hold on;
plot(1:5000,temps(:,1));
plot(1:5000,temps(:,2));
plot(1:5000,temps(:,3));
legend(sprintf('Alpha: %f',alphas(1)),sprintf('Alpha: %f',alphas(2)),sprintf('Alpha: %f',alphas(3)),'location','SouthEast');
xlabel('Time [sec]'); ylabel('Temperature [Celcius]');
title('Temperature vs Time for varying Alphas')
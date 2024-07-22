%% ASEN 3113 LAB 1 MAIN SCRIPT

%{

GROUP: 
    QUINLAN LEWIS
    OWEN KAUFMANN
    EVAN MEANY
    SAMUEL OBERG
    FRED TAYLOR

CREATED: 2020-02-05
LAST MODIFIED: 2020-02-12

DESCRIPTION:

%}

clear; close all; clc

%% LOAD DATA:

fileList = dir(fullfile('Data','Team15*'));

N       = length(fileList); % NUMBER OF FILES
data    = cell(N,1); % INIT

for i = 1:N
    data{i} = readtable(fullfile('Data',fileList(i).name));
    data{i}.Properties.VariableNames = {'t_s','P_PSI','ToT','BoT','ToB','BoB','I','Opt_Switch','Str_Num'};
end

%% CALCULATE RPM:
avg_RPM = zeros(N,1); % INIT
for i = 1:N
    avg_RPM(i) = 1/mean(pulseperiod(data{i,1}.Opt_Switch,data{i,1}.t_s))*60;
end

%% AVERAGE TEMPERATURE DIFFERENCE: USE INNER SURFACES TO CALCULATE
avg_dT = zeros(N,1); % INIT
for i = 1:N
    avg_dT(i) = mean(data{i,1}.ToB-data{i,1}.BoT); 
end

% plot(data{1,1}.P_PSI)

%% PULL IN SOLIDWORKS DATA
fileList2 = dir(fullfile('Soldworks_Data','*AngDisp.csv'));
fileList3 = dir(fullfile('Soldworks_Data','*BigLinDisp.csv'));
fileList4 = dir(fullfile('Solidworks_Data','*SmallLinDisp.csv'));

N2       = length(fileList2); % NUMBER OF FILES
data2    = cell(N2,1); % INIT
N3       = length(fileList3); % NUMBER OF FILES
data3    = cell(N3,1); % INIT
N4       = length(fileList4); % NUMBER OF FILES
data4    = cell(N4,1); % INIT

for i = 1:N2
    data2{i} = readtable(fullfile('Soldworks_Data',fileList2(i).name),'HeaderLines',4);
    data2{i}.Properties.VariableNames = {'Frame','Time','Ang_Disp'};
end
for i = 1:N3
    data3{i} = readtable(fullfile('Soldworks_Data',fileList3(i).name),'HeaderLines',4);
    data3{i}.Properties.VariableNames = {'Frame','Time','Disp'};
end
for i = 1:N4
    data4{i} = readtable(fullfile('Solidworks_Data',fileList4(i).name),'HeaderLines',4);
    data4{i}.Properties.VariableNames = {'Frame','Time','Disp'};
end

%Get flywheel data associated with closest dT
degree7wheel = data2{3};
degree10wheel = data2{1};
degree12wheel = data2{2};

%Get small piston data associated with closest dT
smallpis7 = data4{3};
smallpis10 = data4{1};
smallpis12 = data4{2};

%Get big piston data associated with closest dT
Bigpis7 = data3{5};
Bigpis10 = data3{2};
Bigpis12 = data3{3};

%% P-V Plots
Vmin = 0.000172674;
% hpower = 

idx1 = find(data{1,1}.Opt_Switch==1);
idx2 = find(data{2,1}.Opt_Switch==1);
idx3 = find(data{3,1}.Opt_Switch==1);

%% dT = 7
cycles1 = ceil(mean(idx1((15000<idx1) & (idx1<17000))))-1;

cyclef1 = floor(mean(idx1((17000<idx1) & (idx1<18000))))+1;

P1 = data{1,1}.P_PSI(cycles1:cyclef1).*6894.76;
T1 = data{1,1}.t_s(cycles1:cyclef1);
T1 = T1- T1(1);
Tsmpis = smallpis7.Time(1:160);
Dsmpis = smallpis7.Disp(1:160);
minDsmpis = min(-Dsmpis);
Dsmpis = abs(-Dsmpis - minDsmpis).*10^-3;


D1 = interp1(Tsmpis,Dsmpis,T1);
V1 = Vmin + (pi.*D1.*(7.5.*10^(-3)).^2);

figure()
hold on
grid on
set(gca,'fontsize', 15);
title('P-V Diagram for {\Delta}T = 7 {\circ}C','Interpreter','tex')
plot(V1,P1,'LineWidth',2)
ylim([-400 400])
xlabel('Volume (m^{3})','Interpreter','tex')
ylabel('Pressure (Pa)','Interpreter','tex')
hold off

Wnet1 = polyarea(V1,P1);
%% dT = 10
cycles2 = ceil(mean(idx2((11000<idx2) & (idx2<12000))))-1;

cyclef2 = floor(mean(idx2((12000<idx2) & (idx2<13000))))+1;

P2 = data{2,1}.P_PSI(cycles2:cyclef2).*6894.76;
T2 = data{2,1}.t_s(cycles2:cyclef2);
T2 = T2- T2(1);
Tsmpis = smallpis10.Time(1:160);
Dsmpis = smallpis10.Disp(1:160);
minDsmpis = min(-Dsmpis);
Dsmpis = abs(-Dsmpis - minDsmpis).*10^-3;


D2 = interp1(Tsmpis,Dsmpis,T2);
V2 = Vmin + (pi.*D2.*(7.5.*10^(-3)).^2);

figure()
hold on
grid on
set(gca,'fontsize', 15);
title('P-V Diagram for {\Delta}T = 10 {\circ}C','Interpreter','tex')
plot(V2,P2,'LineWidth',2)
ylim([-400 400])
xlabel('Volume (m^{3})','Interpreter','tex')
ylabel('Pressure (Pa)','Interpreter','tex')
hold off

Wnet2 = polyarea(V2,P2);
%% dT = 12
cycles3 = ceil(mean(idx3((16000<idx3) & (idx3<16500))))-1;
cyclef3 = floor(mean(idx3((16500<idx3) & (idx3<17000))))+1;

P3 = data{3,1}.P_PSI(cycles3:cyclef3).*6894.76;
T3 = data{3,1}.t_s(cycles3:cyclef3);
T3 = T3- T3(1);
Tsmpis = smallpis12.Time(1:160);
Dsmpis = smallpis12.Disp(1:160);
minDsmpis = min(-Dsmpis);
Dsmpis = abs(-Dsmpis - minDsmpis).*10^-3;


D3 = interp1(Tsmpis,Dsmpis,T3);
V3 = Vmin + (pi.*D3.*(7.5.*10^(-3)).^2);

figure()
hold on
grid on
set(gca,'fontsize', 15);
title('P-V Diagram for {\Delta}T = 12 {\circ}C','Interpreter','tex')
plot(V3,P3,'LineWidth',2)
ylim([-400 400])
xlabel('Volume (m^{3})','Interpreter','tex')
ylabel('Pressure (Pa)','Interpreter','tex')
hold off

Wnet3 = polyarea(V3,P3);
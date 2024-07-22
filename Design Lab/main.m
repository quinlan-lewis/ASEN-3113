%Shawn Stone
%ASEN 3113: Thermodynamics
%Created: 4/9/21
clc; clear all; close all;

%% Week 2: Question 1

%Area Factor calculations
t = [0:10:24*60*60];
th = deg2rad(23.5); %Earth's inclination
rotation = (2*pi)/(24*60*60); %[rad/s] rotation rate of spacecraft (360 degrees in 24 hours)
nu = rotation*t;
t_eclipse = [40910, 40910+4580]; %start, end in [s]

%full equations (summer, autumn, winter ,spring - per row)
north = zeros(4, length(t));
south = zeros(4, length(t));
%Zenith Face - "front" at alpha = 0, away from Earth
zenith = zeros(4, length(t));
%Velocity Face - "left hand side" at alpha = 90 deg
vel_face = zeros(4, length(t));
%Nadir Face - "back" at alpha = 180 deg, towards Earth
nadir = zeros(4, length(t));
%Radiator Face - "right hand side" at alpha = 270 deg
radiator = zeros(4, length(t));

for(i = 1:4)
    for(j = 1:length(t))
        north(1, j) = sin(th); %[sin(th);   zeros(i, j)];
        south(3, j) = sin(th); %[zeros(2, length(t));   t_vec*sin(th);  zeros(1, length(t))];
        
        if(  (t_eclipse(1) <= t(j)) && (t(j) <= t_eclipse(2)) && (i == 2 || i == 4) )
            zenith(i, j) = 0;
            vel_face(i, j) = 0;
            nadir(i, j) = 0;
            radiator(i, j) = 0;
        else
            zenith(i, j) = cos(th)*cos(nu(j)); %t_vec* [cos(th)*cos(a);   0;   cos(th)*cos(a);  0]; 
            vel_face(i, j) = cos(th)*cos(nu(j)+pi/2); %t_vec* [cos(th)*cos(a+pi/2);   0;  cos(th)*cos(a+pi/2);  0]; 
            nadir(i, j) = cos(th)*cos(nu(j)+pi); %t_vec* [cos(th)*cos(a+pi);  0;  cos(th)*cos(a+pi)];  0; 
            radiator(i, j) = cos(th)*cos(nu(j)+3*pi/2); %t_vec* [cos(th)*cos(a+3*pi/2);  0;  cos(th)*cos(a+3*pi/2);  0]; 
        end
    end
end

% make negative values zero
zenith = max(zenith,0);
radiator = max(radiator,0);
nadir = max(nadir,0);
vel_face = max(vel_face,0);

% Plot each critical day per face
t2 = t/3600; % time in hours

% Summer Area Factors
title_vec = ["Summer Solstice", "Autumn Equinox", "Winter Solstice", "Spring Equinox"];
for(i = 1:4)
    figure('Units','Normalized','Position',[1/4 1/4 1/2 1/2])
    sgtitle('Area Factor On ' +  title_vec(i));
    subplot(2,3,1)
    plot(t2,north(i,:))
    grid on
    xlim([0 24]);
    xlabel('Time (hr)'); ylabel('Area Factor');
    title('North Face');
    subplot(2,3,2)
    plot(t2, south(i,:))
    grid on
    xlim([0 24]);
    title('South Face');
    subplot(2,3,3)
    plot(t2, zenith(i,:))
    grid on
    xlim([0 24]);
    title('Zenith Face');
    subplot(2,3,4)
    plot(t2, radiator(i,:))
    grid on
    xlim([0 24]);
    title('Radiator Face');
    subplot(2,3,5)
    plot(t2, nadir(i,:))
    grid on
    xlim([0 24]);
    title('Nadir Face');
    subplot(2,3,6)
    plot(t2, vel_face(i,:))
    grid on
    xlim([0 24]);
    title('Velocity Face');
end

%% Finding minimum Surface area
rho_s = [1329, 1374, 1421, 1374];
T = -40 + 273; %[temp in K]
A = zeros(4,6);
A_min = zeros(1,4);
for i =1:4
%passing As as a vecotr into the surfarea function
%order for As will be [Nadir, Radiator, north, south, Zenith, VelFace]
As = [min(nadir(i,:)), min(radiator(i,:)), min(north(i,:)), min(south(i,:)), min(zenith(i,:)), min(vel_face(i,:))];
[Nadir, Radiator, Others] = SurfArea(T,As,rho_s(i));

A(i,1) = Nadir;
A(i,2) = Radiator;
A(i,3:end) = Others(1:4);

A_min(1,i) = max(A(i,:));
end
A_min = max(A_min)
A
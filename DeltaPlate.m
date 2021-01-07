%% DELTA PLATE - TMMV01 (2018)
%% Post-Processing of Force Measurements
%
%% 
close all
clc
clear all

% Number of samples
N_samples = 1000;

% AoA range
aoa_i=4;   % Initial AoA
aoa_f=60;  % Final AoA
daoa =4;  % AoA step
aoa=(aoa_i:daoa:aoa_f); % AoA array
j=1;

%% Read raw data 
load('delta.mat')

%% Plot of measured signals
t=.001:.001:1;

scrsz = get(0,'ScreenSize');
figure('Position',[200 200 scrsz(3)/3 scrsz(3)/4]) 
plot(t,Normal(:,1),'r-',  t,Normal(:,15),'b-',  t,Normal0(:,1),'r:',  t,Normal0(:,15),'b:');
legend('N \alpha=4^{\circ}','N \alpha=60^{\circ}','N0 \alpha=4^{\circ}','N0 \alpha=60^{\circ}')
title('Normal Force');
xlabel('Time [s]');
ylabel('Normal Force [kN]');

%% Calculation of the normal force and coefficient
% We have to correct the normal force because we only want the aerodynamic
% part of it:
N=Normal-Normal0;

figure(2);
plot(aoa,mean(N),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2); % Ploting normal aerodynamic force over the angle of attack.
title('Normal Aerodynamic Force vs Angle of Attack');
xlabel('Angle of attack [degrees (\circ)]');
ylabel('Normal aerodynamic force [kN]');

% we need to calculate the normal force coefficient and for that, we need
% to define some constant that we already know and calculate others:
b = 246/1000; % span in meters
c = 265/1000; % chord in meters
S = c*(b/2); % wing surface in squared meters
T = 20; % temperature of the water our wing is immersed in Celsius degrees
rho = 1000 * (1-((T+288.9414)/(508929.2*(T+68.12963)))*(T-3.9863)^2); % density of water at 20 ºC
V = 0.18; % velocity of the stream in meters/second
q_inf = (1/2)*rho*V^2; % dynamic pressure
CN = N/(q_inf*S); % Normal aerodynamic force coefficient

%% Calculation of the theoretical lift 
% We will apply equation (7.61) of Anderson's book in order to calculate
% this, which says that CL is the sum of CL (potential flow) and CL (vortex
% flow). They both depend on the angle of attack, which we have, and a
% variable defined in Figures (7.42) and (7.43) in the book.
% We will define and calculate some parameters we need to know in order to use those
% plot figures:
AR = 2*b/c; % aspect ratio
a = 0; a_c = a/c; 
sweep_angle = 65; % sweep angle
% Therefore, knowing that for our delta wing: AR = 1.857, a = 0 and sweep
% angle = 65º, we can look for Kp and Kv in the figures, obtaining that:
Kp = 2.10; % planform parameter for potential flow
Kv = 3.15; % planform parameter for vortex flow
% Now, we can apply this to the CL formula:
CLp = Kp * sind(aoa) .* cosd(aoa).^2; % CL for potential flow
CLv = Kv * sind(aoa).^2 .* cosd(aoa); % CL for vortex flow
CL = CLp + CLv; % CL total for our delta wing
%% Plot the results
figure(3);
plot(aoa,mean(CN),'LineWidth',2);
hold on;
plot(aoa,CLp,'LineStyle','--','LineWidth',1);
plot(aoa,CLv,'LineStyle','--','LineWidth',1);
plot(aoa,CL,'LineWidth',2);
title('Normal and Lift Coefficients vs Angle of Attack');
xlabel({'Normal Force Coefficient [-]' ; 'Theoretical Lift Force Coefficient [-]'});
ylabel('Angle of attack [degrees (\circ)]');
legend('CN','CL,p','CL,v','CL');

%--------------------------------
% Test case for overland flow code
%--------------------------------
% Franklin Hinckley
% 7 March 2016
%--------------------------------
%
%--------------------------------

%% Clean up
clearvars
close all
clc


%% Define hillslope
dx = 1; % [m]
len = 1000; % [m]
x = 0:dx:len;

% Parabolic hillslope
zH_par = 5e-5.*(len^2-x.^2);

% Linear hillslope
zH_lin = 100 - 0.02*x;

%% Time array
dt = 0.025; % [s]
tRain = 45*60; % time when rainfall stops

%% Call overland flow subroutine for known case
% Make time array
tSim = tRain; % [s]
t = 0:dt:tSim;

% Call subroutine
[Q, H] = overlandFlow(zH_lin,t,tRain,dx,dt);

% Plot results
figure
plot(x,Q)
xlabel('Position [m]')
ylabel('Flux [m^2/s]')
title('Flux')

figure
plot(x,H)
xlabel('Position [m]')
ylabel('Thickness [m]')
title('Flow Thickness')

%% Check results
R = 0.7e-5; % [m/s]
I = 0;
Qref = 0 + (R - I)*x;

figure
hold on
plot(x,Q)
plot(x,Qref)
hold off
xlabel('Position [m]')
ylabel('Flux [m^2/s]')
title('Flux and Reference')

relDiff = (Q - Qref)./Qref;
figure
plot(x,relDiff)
xlabel('Position [m]')
ylabel('Flux [m^2/s]')
title('Flux: Relative Difference')

%% Call overland flow subroutine to see hydrograph 
% Make time array
tSim = 3*tRain; % [s]
t = 0:dt:tSim;

% Call subroutine
[~,~] = overlandFlow(zH_lin,t,tRain,dx,dt);

%% Call overland flow subroutine with parabolic hillslope
% Make time array
tSim = 3*tRain; % [s]
t = 0:dt/10:tSim;

% Call subroutine
[~,~] = overlandFlow(zH_par,t,tRain,dx,dt);


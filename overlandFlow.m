%--------------------------------
% Simulation of overland flow
%--------------------------------
% Franklin Hinckley
% 7 March 2016
%--------------------------------
% Inputs:
%   zH (n x 1 double array): vector for hillslope heights [m]
%   t (m x 1 double array): time points [s]
%   tRain (double): time point where rainfall stops [s]
%   dx (double): position step [m]
%   dt (double): time step [s]
% Outputs:
%   Q (n x 1 double array): Flux array at final time [m^2/s]
%   H (n x 1 double array): Flow thickness array at final time [m]
%--------------------------------

function [Q, H] = overlandFlow(zH,t,tRain,dx,dt)

%% Constants 
% Flow
kappa = 0.408; % van Karman's Constant
z0 = 2e-3/30; % roughness height [m]
g = 9.81; % standard gravity [m/s^2]

%% Initial thickness
H = zeros(size(zH));

%% Initialize arrays
% % Declare arrays
% H_S = zeros(length(x),length(t));
% Q_S = zeros(length(x),length(t));
% 
% % Assign initial values
% H_S(:,1) = H;
% Q_S(:,1) = zeros(length(x),1);

hydro = zeros(length(t),1);

%% Main loop
% Loop
for ii = 1:length(t)
    % Evaluate rainfall
    if t(ii) < tRain
        R = 0.7e-5; % [m/s] roughly 1 in/hr
    else
        R = 0;
    end
    
    % Evaluate infiltration
    I = 0; % disabled for now
    
    % Define energy slope (assume same as surface slope)
    Se = -diff(zH)/dx;
    
    % Get box-centered heights
    HM = (H(1:end-1) + H(2:end))/2;
    
    % Mean flow velocity (Law of the Wall)
    ubar = zeros(size(HM));
    ubar(HM > 0) = ...
        (sqrt(g*HM(HM > 0).*Se(HM > 0))/kappa).*(log(HM(HM > 0)/z0) - 1);
    
    % Flux 
    Q = ubar.*HM;
    
    % Pad flux array (no flux into left edge)
    Q = [0 Q];
    
    % Error trap (stops simulation if there is a numeric crash)
    if any(isnan(Q))
        error('NaN in Q')
    end
    
    % Flux gradient
    dQdx = diff(Q);
        
    % Thickness rate
    dHdt = -dQdx + R - I;
    
    % Update thickness (assume last height is same as one before)
    H = [H(1:end-1) + dHdt*dt H(end) + dHdt(end)*dt];
    
    % Log hydrograph
    hydro(ii) = Q(end);
    
    % Plot
%     figure(1)
%     subplot(2,1,1)
%     title('Flow','Fontsize',14)
%     plot(x,H_S(:,ii))
%     hold on
%     %plot(x,zH)
%     hold off
%     ylabel('Elevation [m]')
%     xlabel('Position [m]')
%     
%     % Plot flux
%     subplot(2,1,2)
%     % Flux
%     plot(x,Q_S(:,ii))
%     %ylim([0 1e6])
%     ylabel('Flux [m^3/s]')
%     xlabel('Position [m]')
    
    % Save thickness and flux
%     H_S(:,ii) = H;
%     Q_S(:,ii) = Q';    
end

%% Plots
figure
plot(t/60,hydro)
xlabel('Time [min]')
ylabel('Flow at Base [m^2/s]')
ylim([0 8e-3])
title('Hydrograph')



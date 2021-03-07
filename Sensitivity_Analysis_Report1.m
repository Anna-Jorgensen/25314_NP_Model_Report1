%% Sensitivity Testing 

% Grid Values
param.depth = 100; % Water column depth (meters)
nopoints = 100;    % Number of grid cells
param.deltaz= param.depth/nopoints; % cell width (meters)
param.z = param.deltaz*0.5:param.deltaz:(param.depth-param.deltaz*0.5); % center of the cell
param.nGrid = length(param.z);  % No. of grid cells

% Other parameter values
param.D = 10;           % Diffusivity (m2/day) 
param.u = 1;            % Settling velocity (m/day)
param.HI= 20*86400;     % half-saturation constant of light-limitead growth
param.HN=0.0425;        % half sat constant of light
param.m = 0.01*24;      % Loss (mortality) 
param.alpha = 1*10^-9;  % nutrient conc. of P
param.eps = 0.5;        % recycling coef. 
param.NB = 5;           % nutirents at the bottom

tRange = [0:1095];      % Time frame

deltaz = param.deltaz;  % they are the same
z = param.deltaz*0.5:param.deltaz:(param.depth-param.deltaz*0.5);

% These are the parameters we will be testing sensitivity for: 
kp = 6*10^-10;           % Absorp coef. of phytoplankton
Kbg = 0.045;             % background turbidity
u = 1;                   % Settling velocity (m/day)
mumax = 0.04*24;         % max production rate
HI = 20*86400;           % half-saturation constant of light-limitead growth

% ---------- Different Light Incident ---------- %

figure
tiledlayout(2,2);
nexttile
hold on;
legends = [];
for Iint=100*86400:100*86400:600*86400;
    
 % Lets run the model:
[t,z,Y] = Sensitivity(tRange,param,deltaz,z,kp,Kbg,Iint,u,HI,mumax);


N= Y(:,1:param.nGrid);     % split the solution into components
P = Y(:,param.nGrid+1:2*param.nGrid); % split the solution into components


Plast = P(end,:);

[maxvalue,index]=max(Plast);

plot(z(index),maxvalue,'o','linewidth', 20);
ylabel('Phytoplankton Max Conc. Value (cells/m^3)');
xlabel('Depth (meters)');
legends = [legends, 'Light Int. = '+string(Iint)];

end
hold off;
legend(legends);legend('boxoff');
legend('Location','SouthEast');
title('(A) Different incident light intensities (µmol photons/m^2day)');
set(gca,'FontSize',12);

% ----------  Different Kbg (Turbidity)  ---------- %

deltaz = param.deltaz;
z = param.deltaz*0.5:param.deltaz:(param.depth-param.deltaz*0.5);

Iint= 200*86400;  % Light incident
kp = 6*10^-10;    % Absorp coef. of phytoplankton
u = 1;            % Settling velocity (m/day)
HI = 20*86400;    % half-saturation constant of light-limitead growth
mumax = 0.04*24;  % max production rate

nexttile
hold on;
legends = [];
for Kbg=0.02:0.01:0.07;
    
    
 % Lets run the model:
[t,z,Y] = Sensitivity(tRange,param,deltaz,z,kp,Kbg,Iint,u,HI,mumax);


N= Y(:,1:param.nGrid);     % split the solution into components
P = Y(:,param.nGrid+1:2*param.nGrid); % split the solution into components


Plast = P(end,:);

[maxvalue,index]=max(Plast);

plot(z(index),maxvalue,'o','linewidth', 20);
ylabel('Phytoplankton Max Conc. Value (cells/m^3)');
xlabel('Depth (meters)');
legends = [legends, 'Background Turbidity = '+string(Kbg)];

end
hold off;
legend(legends);legend('boxoff');
legend('Location','SouthEast');
title('(B) Different Background Turbidities (m^-1)');
set(gca,'FontSize',12);

% ---------- Different sinking velocity ---------- %


deltaz = param.deltaz;
z = param.deltaz*0.5:param.deltaz:(param.depth-param.deltaz*0.5);

Iint= 200*86400; % Light incident
kp = 6*10^-10;   % Absorp coef. of phytoplankton
Kbg = 0.045;     % Background turbidity
HI = 20*86400;   % half-saturation constant of light-limitead growth
mumax = 0.04*24; % max production rate

nexttile
hold on;
legends = [];
for u=.5:.25:2;
    
    
 % Lets run the model:
[t,z,Y] = Sensitivity(tRange,param,deltaz,z,kp,Kbg,Iint,u,HI,mumax);


N= Y(:,1:param.nGrid);     % split the solution into components
P = Y(:,param.nGrid+1:2*param.nGrid); % split the solution into components


Plast = P(end,:);

[maxvalue,index]=max(Plast);

plot(z(index),maxvalue,'o','linewidth', 20);
ylabel('Phytoplankton Max Conc. Value (cells/m^3)');
xlabel('Depth (meters)');
legends = [legends, 'Sinking Velocity = '+string(u)];

end
hold off;
legend(legends);legend('boxoff');
legend('Location','SouthWest');
title('(C) Different Sinking Velocities (m/day)');
set(gca,'FontSize',12);

% ---------- Different H_I ---------- %

deltaz = param.deltaz;
z = param.deltaz*0.5:param.deltaz:(param.depth-param.deltaz*0.5);

Iint= 200*86400; % light incident
kp = 6*10^-10;   % Absorp coef. of phytoplankton
Kbg = 0.045;     % Background turbidity
u = 1;           % Settling velocity (m/day)
mumax = 0.04*24; % max production rate

nexttile
hold on;
legends = [];
for HI=10*86400:20*86400:80*86400;
    
    
 % Lets run the model:
[t,z,Y] = Sensitivity(tRange,param,deltaz,z,kp,Kbg,Iint,u,HI,mumax);


N= Y(:,1:param.nGrid);     % split the solution into components
P = Y(:,param.nGrid+1:2*param.nGrid); % split the solution into components


Plast = P(end,:);

[maxvalue,index]=max(Plast);

plot(z(index),maxvalue,'o','linewidth', 20);
ylabel('Phytoplankton Max Conc. Value (cells/m^3)');
xlabel('Depth (meters)');
legends = [legends, 'Half Sat. constant (Light) = '+string(HI)];

end
hold off;
legend(legends);legend('boxoff');
legend('Location','NorthWest');
title('(D) Different Half Sat. constant of light limited growth (µmol photons/m^2day)');
set(gca,'FontSize',12);

%% Sensitivity for mumax - not plotting this. only going to show 4 in the report but no one did this one in class. 

% Grid Values
param.depth = 100;
nopoints = 100;
param.deltaz= param.depth/nopoints; % cell width (m)
param.z = param.deltaz*0.5:param.deltaz:(param.depth-param.deltaz*0.5);
param.nGrid = length(param.z);  % No. of grid cells

% Other Values
param.D = 10;  %; % Diffusivity (m2/day) %Same as ex. 2
param.z = param.deltaz*0.5:param.deltaz:(param.depth-param.deltaz*0.5);
param.nGrid = length(param.z);  % No. of grid cells
param.HN=0.0425;
param.m = 0.01*24; % Loss 
param.alpha = 1*10^-9;
param.eps = 0.5;
param.NB = 5; % nutirents at the bottom

tRange = [0:1095]; % Time frame


deltaz = param.deltaz;
z = param.deltaz*0.5:param.deltaz:(param.depth-param.deltaz*0.5);

Iint= 200*86400; % light incident
kp = 6*10^-10;   % Absorp coef. of phytoplankton
Kbg = 0.045;     % Background turbidity
u = 1;           % Settling velocity (m/day)
HI = 20*86400;   % half-saturation constant of light-limitead growth

hold on;
legends = [];
for mumax =0.04*24:0.02*24:0.12*24;
    
    
 % Lets run the model:
[t,z,Y] = Sensitivity(tRange,param,deltaz,z,kp,Kbg,Iint,u,HI,mumax);


N= Y(:,1:param.nGrid);     % split the solution into components
P = Y(:,param.nGrid+1:2*param.nGrid); % split the solution into components


Plast = P(end,:);

[maxvalue,index]=max(Plast);

plot(z(index),maxvalue,'o','linewidth', 20);
ylabel('Phytoplankton Max Conc. Value (cells/m^3)');
xlabel('Depth (meters)');
legends = [legends, 'µmax = '+string(mumax)];

end
hold off;
legend(legends);legend('boxoff');
legend('Location','NorthWest');
title('Different µmax values (1/day)');
set(gca,'FontSize',12);

%% Function for Testing Sensitivity 
% ======================================================
% Run a advection-diffusion equation
% ======================================================

function [t,z,Y] = Sensitivity(tRange, param,deltaz,z,kp,Kbg,Iint,u,HI,mumax) % dz is the grid spacing
%
P0=linspace(0, 0, param.depth/param.deltaz);  % A vector the size of the grid cells
P0(1:length(P0))=1e07;

N0=linspace(0, 0, param.depth/param.deltaz);  % A vector the size of the grid cells
N0(1:length(P0))=13;

% Run model
%
[t, Y] = ode45(@PNmodelDeriv, tRange, [N0 P0], [], param);

z = param.z;

    % ---------------------------------------------------------
    % Derivative function
    % ---------------------------------------------------------
    function dYdt = PNmodelDeriv(t,Y,param)
        %
        N = Y(1:param.nGrid);
        P = Y(param.nGrid+1:2*param.nGrid);
        
         %  ------------ Nutrients Fluxes ------------ %
        
        % Only have mixing (diffusion), no sinking (advection)
        
        ix = 2:(param.nGrid);

        JNdiff(ix) = -param.D*(N(ix)-N(ix-1))/param.deltaz;
        
        JNdiff(1) = 0; % Closed surface boundary...
        
        JNdiff(param.nGrid+1) = -param.D*(param.NB-N(end)/param.deltaz);  % the bottom we have a flux
        
        % Total Flux for Nutrients (just diffusion):
        
        JN =  JNdiff;
        
        %  ------------ Plankton Fluxes ------------ %
        
        % Advective fluxes for P
                JPadv(ix) = u*P(ix-1);
                JPadv(1) = 0;             % Closed surface boundary
                JPadv(param.nGrid+1) = 0; % Closed surface bottom
        
        % Diffusive fluxes for P:
               JPdiff(ix) = -param.D*(P(ix)-P(ix-1))/param.deltaz;
               JPdiff(1) = 0;             % No flux at the surface...
               JPdiff(param.nGrid+1) = 0; % ...or the bottom
        
        % Total Flux for Plankton (adv + diff):
        
        JP = JPadv + JPdiff;
        
        % Growth of P
        
        I=calclight(z,P,Iint,deltaz,kp,Kbg);  % calc I --> light
        
        N=N'; % change to column vector
        
        mu = mumax*min((N./(N+param.HN)),(I./(I+HI)));    % now the functional response
            
        % Calc. Derivatives
        
        % ------------------- Nutrients ------------------- %
        
          %          -       uptake   +       % recycling               % - % mixing                                                      
        dNdt =   - param.alpha*mu.*P' + param.eps*param.alpha*param.m.*P' - ((JN(2:(param.nGrid+1))-JN(1:param.nGrid))/param.deltaz);
       % ------------------- Phytoplankton ------------------- %
        
                        % - sinking and mixing                           + growth - loss/death
        dPdt = (-(JP(2:(param.nGrid+1))-JP(1:param.nGrid))/param.deltaz) + mu.*P' - param.m.*P';

        % Make dPdt a column vector:
        dYdt = [dNdt dPdt]';
        
    end

end


%%
function I=calclight(z,P,Iint,deltaz,kp,Kbg)

n = length(z); % No. of grid cells

int=0;

for i=2:n;
  
  int(i)=int(i-1)+kp*P(i)*deltaz; % sum of kpPDeltaZ
  
end
    
 I=Iint*exp(-Kbg*z).*exp(-int);  % now calc I (light) at the diff grids

 I(1) = Iint; % First cell should be the initial light value
 
end
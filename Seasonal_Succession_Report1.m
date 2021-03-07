%% Seasonal Sucession Plots

% Parameters:

% Grid Values
param.depth = 100;                  % Water column depth (meters)
nopoints = 100;                     % Number of grid cells
param.deltaz= param.depth/nopoints; % cell width (m)
param.z = param.deltaz*0.5:param.deltaz:(param.depth-param.deltaz*0.5);
param.nGrid = length(param.z);      % No. of grid cells

% Other  Parameter Values
param.D = 10;                     % Diffusivity (m2/day) %Same as ex. 2
param.u = 1;                     % Sinking velocity (m/day)
param.HI= 20*86400;              % half-saturation constant of light-limitead growth
param.HN=0.0425;                 % half sat constant of light
param.m = 0.01*24;               % Loss Mortality 
param.alpha = 1*10^-9;           % nutrient conc. of P
param.eps = 0.5;                 % recycling coef. 
param.NB = 5;                    % nutirents at the bottom
param.mumax= 0.9600;             % max production rate

% Light function values 
param.Iint= 300*86400;  % Initial light intensity
param.kp= 6*10^-10;     % specific light attenuation of phytoplankton
param.Kbg = 0.045;      % background turbidity

tRange = [0:1825];      % Time frame

% Lets run the model:
[t,z,Y] = PNmodel_season(tRange, param);

N= Y(:,1:param.nGrid);                % split the solution into Nutrient
P = Y(:,param.nGrid+1:2*param.nGrid); % split the solution into Phytoplankton

%% Plot - Diffusion and Light Seasonality

clf
t = tiledlayout(1,2,'TileSpacing','Compact');
nexttile
startDate = datenum('01-01-2020');
endDate = datenum('12-31-2020');
xData = linspace(startDate,endDate,12);
D = @(t) param.D*(1+cos(2*pi*(t/365)))
fplot(D,[0 730],'b','LineWidth',2.5)
datetick('x','mmm')
xlabel('Time (month)')
ylabel('Diffusion (m^2/day)')
set(gca,'FontSize',15)
title('(A) Diffusion')

% We want D to be the highest in Jan and lowest in Jul

% Plot - Light Seasonality
nexttile
startDate = datenum('01-01-2020');
endDate = datenum('12-31-2020');
xData = linspace(startDate,endDate,12);
Iint = @(t) param.Iint*(1-.8*cos(2*pi*(t/365))) 
fplot(Iint,[0 730],'m','LineWidth',2.5)
datetick('x','mmm')
xlabel('Time (month)')
ylabel('Light Intensity (µmol photons/m^2day)')
title('(B) Light Intensity')
set(gca,'FontSize',15)

% We want Light to be the highest in Jul and lowest in Jan

%% Nice Surface Plots 

figure
tiledlayout(2,2);
nexttile
surface(t,-z,P')
colorbar
xlabel('Time (days)')
ylabel('Depth (meters)')
title('(A) Phytoplankton surface plot over time')
set(gca,'FontSize',15)
shading interp
axis tight

nexttile
surface(t,-z,N')
colorbar
xlabel('Time (days)')
ylabel('Depth (meters)')
title('(B) Nutrient surface plot over time')
set(gca,'FontSize',15)
shading interp
axis tight

nexttile
plot(t,P(:,z==.5),'DisplayName',' Surface (1m depth)','LineWidth',3)
hold on
plot(t,P(:,z==49.5),'DisplayName','Middle (50m depth)','LineWidth',3)
hold on
plot(t,P(:,z==99.5),'DisplayName','Bottom (100m depth)','LineWidth',3)
xlabel('Time (days)')
ylim([-0 2.5e09])
set(gca,'FontSize',15)
ylabel('Phytoplankton Concentration (cells/m^3)')
title('(C) Phytoplankton concentration over time')
xlim([0 1825])
hold on
hold off
legend('Location','NorthEast')

nexttile
plot(t,N(:,z==.5),'DisplayName',' Surface (1m depth)','LineWidth',3)
hold on
plot(t,N(:,z==49.5),'DisplayName','Middle (50m depth)','LineWidth',3)
hold on
plot(t,N(:,z==99.5),'DisplayName','Bottom (100m depth)','LineWidth',3)
xlabel('Time (days)')
ylabel('Nutrient Concentration (mmol nutrients/m^2day)')
set(gca,'FontSize',12)
title('(D) Nutrient concentration over time')
xlim([0 1825])
hold off
legend('Location','NorthEast')
set(gca,'FontSize',15)

%% What is Limiting in the last time step - winter?

figure
x1 =  P(end,:);
y1 = param.z;
plot(x1,y1,'g','LineWidth',1.5,'DisplayName','Phytoplankton')
ylabel('Depth (meters)')
xlabel('Phytoplankton Concentration (cells/m^3)')

ax1 = gca; % current axes
ax1.XColor = 'b';
ax1.YColor = 'b';
axis ij

ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');

set(gca,'ytick',[]);
set(gca,'XLim',[0 1])
set(gca,'XTick',(0:0.5:1))

x2 = 0:.5:1;

xline(param.m/param.mumax,'k--','DisplayName','m mumax^-1','LineWidth',1.5)
hold on

xlim([0 1])

NLim = N./(N+param.HN)
NLim2 = NLim(end,:);

P0 = P(end,:);
I1 = Paramcalclight(P0,param)*(1-.8*cos(2*pi*(365/365))) % calc I --> light %day 365 is last day of dec. and ~185 day is jul
ILim1 = I1./(I1+param.HI)
ILim1'

plot(ILim1,y1,'m-.','LineWidth',1.5,'DisplayName',' Limiting Light Term')
axis ij
hold on
plot(NLim2,y1,'r-.','LineWidth',1.5,'DisplayName','Nutrient Limiting Term')
legend('Location','SouthEast')
hold off
title('(A) Limitation by light and nutrients: winter season')

%% What is Limiting in the last time step - summer?

% summer day is in calender ~ day 185, 365-185 = 180, 1825-180...1645
% so the last time step the summer is day 1645

figure
x1 = P(1645,:);
y1 = param.z;
plot(x1,y1,'g','LineWidth',1.5,'DisplayName','Phytoplankton');
ylabel('Depth (meters)');
xlabel('Phytoplankton Concentration (cells/m^3)');

ax1 = gca; % current axes
ax1.XColor = 'b';
ax1.YColor = 'b';
axis ij

ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');

set(gca,'ytick',[]);
set(gca,'XLim',[0 1]);
set(gca,'XTick',(0:0.5:1));
x2 = 0:.5:1;


xline(param.m/param.mumax,'k--','DisplayName','mumax^-1','LineWidth',1.5);
hold on;

xlim([0 1]);

NLim = N./(N+param.HN);
NLim3 = NLim(1645,:);

P0 = P(1645,:);
I2 = Paramcalclight(P0,param)*(1-.8*cos(2*pi*(185/365))); % calc I --> light %day 365 is last day of dec. and ~185 day is jul
ILim2 = I2./(I2+param.HI);
ILim2';

plot(ILim2,y1,'m-.','LineWidth',1.5,'DisplayName',' Limiting Light Term');
axis ij;
hold on;
plot(NLim3,y1,'r-.','LineWidth',1.5,'DisplayName','Nutrient Limiting Term');
legend('Location','SouthEast');
title('(B) Limitation by light and nutrients: summer season');


%% Seasonality Function
% ======================================================
% Run a advection-diffusion equation
% ======================================================

function [t,z,Y] = PNmodel_season(tRange, param) % dz is the grid spacing
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

        JNdiff(ix) = -param.D*(1+cos(2*pi*(t/365)))*(N(ix)-N(ix-1))/param.deltaz;
        
        JNdiff(1) = 0; % No flux at the surface...
        
        JNdiff(param.nGrid+1) = -param.D*(1+cos(2*pi*(t/365)))*(param.NB-N(end)/param.deltaz);  % the bottom we have a flux
        
        % Total Flux for Nutrients (just diffusion):
        
        JN =  JNdiff;
        
        
        %  ------------ Plankton Fluxes ------------ %
        
        % Advective fluxes for P
                JPadv(ix) = param.u*P(ix-1);
                JPadv(1) = 0;  % No input from the surface
                JPadv(param.nGrid+1) = 0; % Closed bottom
        
        % Diffusive fluxes for P:
               JPdiff(ix) = -param.D*(1+cos(2*pi*(t/365)))*(P(ix)-P(ix-1))/param.deltaz;
               JPdiff(1) = 0;              % No flux at the surface...
               JPdiff(param.nGrid+1) = 0;  % ...or the bottom
        
        % Total Flux for Plankton (adv + diff):
        
        JP = JPadv + JPdiff;
        
        % Growth of P
        
        I = Paramcalclight(P,param)*(1-.8*cos(2*pi*(t/365)));  % calc I --> light with seasonality
        
        N=N'; % change to column vector
        
        mu = param.mumax*min((N./(N+param.HN)),(I./(I+param.HI)));    % now the functional response
            
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

%% Ligth Function

function I = Paramcalclight(P,param)

n = length(param.z); % grid cells

int=0;

for i=2:n;
  
  int(i)=int(i-1)+param.kp*P(i)*param.deltaz; % sum of kpPDeltaZ
  
end
    
 I=param.Iint*exp(-param.Kbg*param.z).*exp(-int);  % now calc I at the diff grids

 I(1) = param.Iint; % First cell should be the initial light value
 
end

%% Function for plotting limiting light

function I = ParamcalclightP0(P0,param)

n = length(param.z); % grid cells

int=0;

for i=2:n;
  
  int(i)=int(i-1)+param.kp*P0(i)*param.deltaz; % sum of kpPDeltaZ
  
end
    
 I=param.Iint*exp(-param.Kbg*param.z).*exp(-int);  % now calc I at the diff grids

 I(1) = param.Iint % First cell should be the initial light value
 
end

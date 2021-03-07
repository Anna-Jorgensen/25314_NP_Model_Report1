%%  Basic Model Output _ No Seasonality

% Parameters:

% Grid Values
param.depth = 100;
nopoints = 100;
param.deltaz= param.depth/nopoints; % cell width (m)
param.z = param.deltaz*0.5:param.deltaz:(param.depth-param.deltaz*0.5);
param.nGrid = length(param.z);  % No. of grid cells

% Other Values
param.D = 10;           % Diffusivity (m2/day) 
param.u = 1;            % Settling velocity (m/day)
param.HI= 20*86400;     % half-saturation constant of light-limitead growth
param.HN=0.0425;        % half sat constant of light
param.m = 0.01*24;      % Loss 
param.alpha = 1*10^-9;  % nutrient conc. of P
param.eps = 0.5;        % recycling coef. 
param.NB = 5;           % nutirents at the bottom
param.mumax= 0.96;    % max production rate

% Light function values 
param.Iint= 600*86400;  % Initial light intensity
param.kp= 6*10^-10;     % specific light attenuation of phytoplankton
param.Kbg = 0.045;      % background turbidity

tRange = [0:1460];      % Time frame

% Lets run the model:
[t,z,Y] = NO_seas_PNmodel(tRange, param);

N= Y(:,1:param.nGrid);                % split the solution into Nutrient
P = Y(:,param.nGrid+1:2*param.nGrid); % split the solution into Phytoplankton

%% Lets plot the figure -> surface plot

figure;
tl = tiledlayout(2,2);
nexttile;
surface(t,z,P');
colorbar;
xlabel('Time (days)');
ylabel('Depth (meters)');
title('(A) Phytoplankton surface plot over time');
shading interp;
axis ij;
axis tight;

nexttile;
surface(t,z,N');
colorbar;
xlabel('Time (days)');
ylabel('Depth (meters)');
title('(B) Nutrient surface plot over time');
shading interp;
axis ij;
axis tight;

nexttile;
plot(t,P(:,z==.5),'DisplayName','Surface (1m depth)','LineWidth',3);
hold on;
plot(t,P(:,z==49.5),'DisplayName','Middle (50m depth)','LineWidth',3);
hold on;
plot(t,P(:,z==99.5),'DisplayName','Bottom (100m depth)','LineWidth',3);
title('(C) Phytoplankton concentration over time');
legend('Location','NorthEast');
ylabel('Phytoplankton Concentration (cells/m^3)');
xlabel('Time (days)');
%xline(2555,'k--','DisplayName','Converged (end of year 7)','LineWidth',1.5)
xlim([0 1460])

nexttile;
plot(t,N(:,z==.5),'DisplayName','Surface (1m depth)','LineWidth',3);
hold on;
plot(t,N(:,z==49.5),'DisplayName','Middle (50m depth)','LineWidth',3);
hold on;
plot(t,N(:,z==99.5),'DisplayName','Bottom (100m depth)','LineWidth',3);
xlabel('Time (days)')
ylabel('Nutrient Concentration (mmol nutrients/m^2day)');
title('(D) Nutrient concentration over time');
legend('Location','NorthEast');
%xline(2555,'k--','DisplayName','Converged (end of year 7)','LineWidth',1.5)
xlim([0 1460])

title(tl,'Basic Model Output: No seasonal effects');

%% Delta Indep?

% Grid Values
param.depth = 100;      % Depth (meters)

% Other Parameter Values
param.D = 10;           % Diffusivity (m2/day) 
param.u = 1;            % Settling velocity (m/day)
param.HI= 20*86400;     % half-saturation constant of light-limitead growth
param.HN=0.0425;        % half sat constant of light
param.m = 0.01*24;      % Loss (mortality) 
param.alpha = 1*10^-9;  % nutrient conc. of P
param.eps = 0.5;        % recycling coef.
param.NB = 5;           % nutirents at the bottom
param.mumax= 0.9600;    % max production rate

% Light function values 
param.Iint= 600*86400; % Initial light intensity
param.kp= 6*10^-10;    % specific light attenuation of phytoplankton
param.Kbg = 0.045;     % background turbidity

tRange = [0:3000];     % Time frame

%%%%%% delta is now 1 %%%%%
nopoints = 100;
param.deltaz= param.depth/nopoints; % cell width (m)
param.z = param.deltaz*0.5:param.deltaz:(param.depth-param.deltaz*0.5);
param.nGrid = length(param.z);  % No. of grid cells

% Lets run the model:
[t,z,Y] = NO_seas_PNmodel(tRange, param,nopoints);

P1 = Y(:,param.nGrid+1:2*param.nGrid); 
P1_end=P1(end,:);
plot(P1_end,-z,'m:','LineWidth',1.5,'DisplayName','Delta = 1');
hold on;

%%%%%% delta is now 2 %%%%%
nopoints = 50;
param.deltaz= param.depth/nopoints; 
param.z = param.deltaz*0.5:param.deltaz:(param.depth-param.deltaz*0.5);
param.nGrid = length(param.z);  

% Lets run the model:
[t,z,Y] = NO_seas_PNmodel(tRange, param,nopoints);

P2 = Y(:,param.nGrid+1:2*param.nGrid); 
% Lets plot the figure
P2_end=P2(end,:);
plot(P2_end,-z,'g--','LineWidth',1.5,'DisplayName','Delta = 2');
hold on;

%%%%% delta is now .5 %%%%%
nopoints = 200;
param.deltaz= param.depth/nopoints; 
param.z = param.deltaz*0.5:param.deltaz:(param.depth-param.deltaz*0.5);
param.nGrid = length(param.z);  

% Lets run the model:
[t,z,Y] = NO_seas_PNmodel(tRange, param,nopoints);

P3 = Y(:,param.nGrid+1:2*param.nGrid); % split the solution into components
P3_end=P3(end,:);
plot(P3_end,-z,'r*','LineWidth',1.5,'DisplayName','Delta = .5');
hold on 

%%%%% delta is now .75 %%%%%
nopoints = 75;
param.deltaz= param.depth/nopoints; 
param.z = param.deltaz*0.5:param.deltaz:(param.depth-param.deltaz*0.5);
param.nGrid = length(param.z);  

% Lets run the model:
[t,z,Y] = NO_seas_PNmodel(tRange, param,nopoints);

P4 = Y(:,param.nGrid+1:2*param.nGrid); 
P4_end=P4(end,:);

plot(P4_end,-z,'c-.','LineWidth',1.5,'DisplayName','Delta = .75');
hold off;
xlabel('Phytoplankton Concentration (cells/m^3)');
ylabel('Depth (meters)');
legend('Location','NorthEast');
title('Phytoplakton depth profile at the last time step when converged');

%% What is Limiting in the last time step ?

% Grid Values
param.depth = 100;
nopoints = 100;
param.deltaz= param.depth/nopoints; % cell width (m)
param.z = param.deltaz*0.5:param.deltaz:(param.depth-param.deltaz*0.5);
param.nGrid = length(param.z);  % No. of grid cells

% Other Values
param.D = 10;           % Diffusivity (m2/day) 
param.u = 1;            % Settling velocity (m/day)
param.HI= 20*86400;     % half-saturation constant of light-limitead growth
param.HN=0.0425;        % half sat constant of light
param.m = 0.01*24;      % Loss 
param.alpha = 1*10^-9;  % nutrient conc. of P
param.eps = 0.5;        % recycling coef. 
param.NB = 5;           % nutirents at the bottom
param.mumax= 0.9600;    % max production rate

% Light function values 
param.Iint= 600*86400;  % Initial light intensity
param.kp= 6*10^-10;     % specific light attenuation of phytoplankton
param.Kbg = 0.045;      % background turbidity

tRange = [0:1460];      % Time frame

% Lets run the model:
[t,z,Y] = NO_seas_PNmodel(tRange, param,nopoints);

% Nice Plots 

N= Y(:,1:param.nGrid);      % split the solution into components
P = Y(:,param.nGrid+1:2*param.nGrid); % split the solution into components

figure

x1 = P(end,:);
y1 = param.z;
plot(x1,y1,'g','LineWidth',1.5,'DisplayName','Phytoplankton');
ylabel('Depth (meters)');
xlabel('Phytoplankton Concentration (cells/m^3)');

ax1 = gca; % current axes
ax1.XColor = 'b';
ax1.YColor = 'b';
axis ij;

ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');

set(gca,'ytick',[]);
set(gca,'XLim',[0 1]);
set(gca,'XTick',(0:0.5:1));

x2 = 0:.5:1;


xline(param.m/param.mumax,'k--','DisplayName','m mumax^-1','LineWidth',1.5);
hold on;

xlim([0 1]);

NLim = N./(N+param.HN);
NLim2 = NLim(end,:);


P0 = P(end,:);
I1 = Paramcalclight(P0,param); 
ILim1 = I1./(I1+param.HI);
ILim1'

plot(ILim1,y1,'m-.','LineWidth',1.5,'DisplayName',' Limiting Light Term')
axis ij
hold on
plot(NLim2,y1,'r-.','LineWidth',1.5,'DisplayName','Nutrient Limiting Term')
legend('Location','SouthEast')
hold off
title('Limitation by light and nutrients: no seasonality')

%% No Seasonality Function
% ======================================================
% Run a advection-diffusion equation
% ======================================================

function [t,z,Y] = NO_seas_PNmodel(tRange, param,nopoints) % dz is the grid spacing
%
P0=linspace(0,0, param.depth/param.deltaz);  % A vector the size of the grid cells
P0(1:length(P0))=1e07;

N0=linspace(0,0, param.depth/param.deltaz);  % A vector the size of the grid cells
N0(1:length(N0))=13;

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
        
        JNdiff(1) = 0; % No flux at the surface...
        
        JNdiff(param.nGrid+1) = -param.D*(param.NB-N(end)/param.deltaz);  % the bottom we have a flux dep. on the NB
        
        % Total Flux for Nutrients (just diffusion):
        
        JN =  JNdiff;
        
        %  ------------ Plankton Fluxes ------------ %
        
        % Advective fluxes for P
                JPadv(ix) = param.u*P(ix-1);
                JPadv(1) = 0;             % No input from the surface
                JPadv(param.nGrid+1) = 0; % Closed bottom
        
        % Diffusive fluxes for P:
               JPdiff(ix) = -param.D*(P(ix)-P(ix-1))/param.deltaz;
               JPdiff(1) = 0;              % No flux at the surface...
               JPdiff(param.nGrid+1) = 0;  % ...or the bottom
        
        % Total Flux for Plankton (adv + diff):
        
        JP = JPadv + JPdiff;
        
        % Growth of P
        
        I = Paramcalclight(P,param);  % calc I --> light
        
        N=N'; % change to column vector
        
        mu = param.mumax*min((N./(N+param.HN)),(I./(I+param.HI)));    % now the functional response
            
        % Calc. Derivative
        
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
function I = Paramcalclight(P,param)

n = length(param.z); % grid cells

int=0;

for i=2:n;
  
  int(i)=int(i-1)+param.kp*P(i)*param.deltaz; % sum of kpPDeltaZ
  
end
    
 I=param.Iint*exp(-param.Kbg*param.z).*exp(-int);  % now calc I at the diff grids

 I(1) = param.Iint; % First cell should be the initial light value
 
end
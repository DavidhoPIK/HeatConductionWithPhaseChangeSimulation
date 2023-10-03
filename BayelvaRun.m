% Bayelva permafrost research site

%% Setup Simulation
load("Ny_Alesund_ERA5_downscaled_3h_original.mat")

Tair=FORCING.data.Tair; % air temperature as dirichlet condition
Tair_handl= @(x) Tair(floor(x/(3*60*60)+1))
k=20; % num of gidpoints frst layer; 100 used for masterthesis
u=20; % " second layer
s=20; % " third layer
grid= create_grid(k,u,s)

%phiscal constants that are from cryogrid constants
c_w=4200000; %heat cap water
c_i=1900000; %heat cap ice
c_o= 2500000;% " organic components
c_m=2000000; % mineral comp
L_w= 334000000; %latent heat water
k_w=0.57; % thermal conductivity water
k_i=2.2; % " ice
k_o=0.25; % " organic comp
k_m=3; % " minerl

%set up values for each layer
%layer 1,2 have 0.4 water, 0.6 minerals, layer three has 0.03 water, 0.97 minerals
c_f_layer= [0.4*c_i+0.6*c_m, 0.4*c_i+0.6*c_m, 0.03*c_i+0.97*c_m]; % stores values for each of the three layes
c_u_layer= [0.4*c_w+0.6*c_m, 0.4*c_w+0.6*c_m, 0.03*c_w+0.97*c_m];
k_f_layer= [k_w^(0.4)+k_m^(0.6), k_w^(0.4)+k_m^(0.6), k_w^(0.03)+k_m^(0.97)]; % geometric mean
k_u_layer= [k_i^(0.4)+k_m^(0.6), k_i^(0.4)+k_m^(0.6), k_i^(0.03)+k_m^(0.97)];
L_layer=[0.4*L_w, 0.4*L_w, 0.03*L_w];
ph= create_physical_data(k,u,s, c_f_layer,c_u_layer,k_f_layer,k_u_layer,L_layer); % puts the physical data to the grid

%% Run simulation
sim = StefanSim2(grid', ... %grid
                 60*60*24*365*4, ... %simulation length
                 ph, ... % physical data
                 Tair_handl, ... % dirichlet condition
                 @(x)0, ... % nuemann condition 
                 repmat(-10000000,length(grid)-1,1), ... % initial condition
                 50, ... % factor to define scale the timestep
                 0, ... % do not plot diagnose info while while running
                 1) % use backward euler 

sim.initialize

sim.simulate % core simulation happening here

%% Plot Results 1, Lineplot
starttime=FORCING.data.t_span(1);
plot((1:floor(length(Tair)/18))/(8)+starttime, Tair(1:floor(length(Tair)/18)), color="#909090" )
hold on
p1= sim.plotThist(0.2,starttime)
set(p1,'Color','#f5cd5f', LineWidth=1)
p2=sim.plotThist(0.6,starttime)
set(p2,'Color','#AF531E', LineWidth=1)
%p2=sim.plotThist(1,"days") 
%set(p2,'Color','#FF531E', LineWidth=0.9)
p3=sim.plotThist(1.5,starttime)
set(p3,'Color','#494C66', LineWidth=1)
p4=sim.plotThist(6,starttime)
set(p4,'Color','#32A691', LineWidth=1)
hold off
%}
line(xlim, [0,0], 'Color', '#959595', 'LineStyle', '--', LineWidth=0.5 );

legend("Tair","0.2m","0.6m","1.5m","6m")

xlim([starttime-1,starttime+2*365]);
datetick('x','mm/yy','keeplimits');
ylim([-40,15])


%% Plot Results 2, 2d tile plot
T=TH((sim.S),ph.c_fro(1:end),ph.c_nor(1:end),ph.L(1:end)); % get temperature

% Assume timesteps are equally spaced, we create vectors for them.
time_indays = ((1:size(T, 2))*sim.dt)/(60*60*24);

% Add an extra column and row filled with NaN values
T_pcolor = [T, NaN(size(T, 1), 1)];
T_pcolor = [T_pcolor; NaN(1, size(T_pcolor, 2))];

% Create the figure and the axes
figure;
axes;

% Create the pseudocolor plot with colorbar
pcolor([time_indays, time_indays(end)+1]+starttime, [grid(2:end), grid(end)+1], T_pcolor);

% Adding vertical lines and labels for each year
num_years = 4;
for i = 1:num_years
    % Add a vertical line at the start of each year
    line([i*365+starttime, i*365+starttime], ylim, 'Color', '#E5E5E5', 'LineStyle', '-', LineWidth=1.4 );
    
    % Add a label for each year
    text(i*365, 0, num2str(i), 'VerticalAlignment', 'top');
end
% Configure the y-axis to be in reverse order
set(gca, 'YDir','reverse');

% Adding colorbar
% Define the transition points
transitions = [0, 0.6, 1];  % normalized transition points

% Define the colors
colors = [0, 33, 80;
           199, 199, 168;
            255, 233, 100 % blue
            % white
          ]/255; % red


% Create the colormap
custom_colormap = interp1(transitions, colors, 0:0.01:1);
custom_colormap = imadjust(custom_colormap, [0.1; 0.8], []);
colormap(custom_colormap);
c=colorbar;
title(c,'°C')

%brighten(0.2)
%contrast(8)

% Set y-axis limits to display only the first few meters
ylim([0 10]);

% Label axes
%xlabel('Time in days');
ylabel('Soil depth in m');

% Title
%title('Evolution of soil temperatures');

% Remove grid lines
shading flat;

xlim([starttime-1,starttime+4*365]);
datetick('x','mm/yy','keeplimits');

%% Helper functions
function grid = create_grid(k, s, u)
    % create a vector of equally spaced points between 0 and 0.1
    grid1 = linspace(0, 1, k+1); 

    % create a vector of equally spaced points between 0.1 and 5
    % don't repeat the point 0.1, hence start from grid2(2)
    grid2 = linspace(1, 15, s+1);
    grid2 = grid2(2:end);

    % create a vector of equally spaced points between 5 and 100
    % don't repeat the point 5, hence start from grid3(2)
    grid3 = linspace(15, 100, u+1);
    grid3 = grid3(2:end);

    % concatenate all vectors to get the final grid
    grid = [grid1, grid2, grid3];
end

function pd = create_physical_data(k,u,s,c_fro_layer,c_un_layer,k_fro_layer,k_un_layer,L_layer)
    % grid is the vector of grid points
    % each layer parameter is a 1x3 vector, representing the value for each layer
    c_fro = [repmat(c_fro_layer(1),k,1); repmat(c_fro_layer(2),u,1); repmat(c_fro_layer(3),s,1)];
    c_un = [repmat(c_un_layer(1),k,1); repmat(c_un_layer(2),u,1); repmat(c_un_layer(3),s,1)];
    L = [repmat(L_layer(1),k,1); repmat(L_layer(2),u,1); repmat(L_layer(3),s,1)];
    k_fro = [repmat(k_fro_layer(1),k,1); repmat(k_fro_layer(2),u,1); repmat(k_fro_layer(3),s,1)];
    k_un = [repmat(k_un_layer(1),k,1); repmat(k_un_layer(2),u,1); repmat(k_un_layer(3),s,1)];
    pd = physicalData(c_fro, c_un, k_fro, k_un, L);
end

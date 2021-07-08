%% GRUBER: terrain based approximation of gravitationnally driven snow transport
%%% Based on Gruber 2007.                                                    
%%% Incorporates snow drift following winstral 2002 (INIT_WINSTRAL function
%%% by Rickard Pettersson). 

%% Yoram Terleth 2021 
%%% example weather data from SMHI tarfala & nikkaluokta wx stations
%%% wind direction from ECMWF ERA 5 re-analysis
%%% DEM from Mercer 2016 and Stockholm University Bolin Centre 
%%% for full references and model explanation please find Terleth 2021,
%%% MSc. Thesis
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INIT parameters 

clearvars

% Time parameters
time.ts = '1-Sep-2007 00:00';                                              % Date and time start run
time.te = '15-Sep-2007 00:00';                                             % Date and time end run
time.dt = .5;                                                             % Timestep (days)
time.tn = round((datenum(time.te)- ...                                     % Nr. of time-steps
    datenum(time.ts))/time.dt)+1;
time.dT_UTC = 1;                                                            % time difference relative to UTC (hours)          

% DEM to be used 
H.D_elev_m = 'tarfala.mat' ; 

% headwall dem parameters
H.cs = 10 ;                                                                   % cell size of DEM
H.x_bounds = [396150 399700];                                                 % UTM easting limits of zone to compute grav. transport for
H.y_bounds = [7533150 7535550];                                               % UTM northing bounds of zone to compute grav. transport for

% Wind drift parameters
WD.thresholdWS = 5 ;                                                         % Threshold wind speed for snow drift 
WD.Sdistance = 20;                                                          % the distance to which terrain is sheltered: 40m or 100m

WD.windDirs = 270.1 ; % e.g.[0.1:22.5:360]  ;                                     % wind directions to consider

% initialise new wind sheltering file?
NewWS = 'no' ; % 'yes' or 'no' 

% Avalanche parameters beta max and d max 
H.Blim =30 ;       % [degrees]                                                              % see willibald et al. 2020 for angle of repose
H.Dlim = .05 ;       % [ m w e ]                                                              % [ m we ] % to be est.: look at gruber: should be no prob to use a large
                                                                                           %  value if you use the max runnout inclusion
H.runout_a = 5 ;   % [degrees]                                                              % runnout angle: this is the lowest angle at which deposits occur, "headwall" only considers terrain above this angle 

% constants
C.precipZ = 0.3 ; % linear precip. / elevation grad.
C.rainsnowT = 272.16 ; % rain snow transition temp.

% climate data: 
% should contain precip, wind speed, air temp at temp. resolution below timestep  
load(['wxdata.mat']);

%% GRID 
[H,DEM] = INIT_GRUBER_headwall(H);

%% WIND 

% re-initialise terrain derived shetering
if strcmp(NewWS,'yes')
[WD]= INIT_WINSTRAL(DEM,WD,H); 
end 

% read in sheltering index
[WD] = INIT_GRUBER_wind_drift(WD); 

%% DYNAMIC PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:time.tn
    
    %% PRINT TIME
    [time] = GRUBER_print_time(t,time);
    
    %% CLIMATE INPUT 
    [H] = GRUBER_headwall_climate(H,time,C,WD,DEM,wx); 
    
    %% GRAVITATIONAL TRANSPORT 
    [H] = GRUBER_headwall(H,C); 
    
    %% Variable Storage
    Mass(:,:,t) = H.M ;
    Acc(:,:,t) = H.Acc ;
    Dep(:,:,t) = H.D ;
    
 
end 

%% sum of steps 
totalM = sum(Mass,3); 
totalD = sum(Dep,3); 
totalAcc = sum(Acc,3); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VISUALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
var_to_plot = totalAcc; % can be modified to ex Fnb(:,:,3) to visualise inter. steps
z1=[min(DEM.z(:)):20:max(DEM.z(:))]; % contour width
var_to_plot(var_to_plot==0) = NaN ; 
box on, 
contour(DEM.x,DEM.y,DEM.z,z1,'LineWidth',0.5,'color',[0,0,0]+0.5),hold on 
S = pcolor(DEM.x,DEM.y,var_to_plot); shading flat;
cRange = S.CData ; 
alpha 0.5
cmp = (abs(cbrewer('seq','YlGnBu',128,'spline')));
cmp(cmp>1)=1;
colormap(cmp)
c = colorbar; 
caxis([min(cRange(:)),max(cRange(:))]);
ylabel(c,'transferred mass fraction');
axis equal
xlim(H.x_bounds),ylim(H.y_bounds)
xlabel('UTM easting (m)'),ylabel('UTM northing (m)')
set(gca, 'Box', 'on','LineWidth',1);
ax = gca ; 
ax.XColor = 'k';
ax.YColor = 'k'; 

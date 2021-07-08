%% AVALANCHE
%% Run Gruber 2007 independently

%% INIT parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars

% Time parameters
time.ts = '1-Feb-2007 00:00';                                              % Date and time start run
time.te = '15-Feb-2007 00:00';                                              % Date and time end run
time.dt = .25;                                                               % Timestep (days)
time.tn = round((datenum(time.te)- ...                                      % Nr. of time-steps
    datenum(time.ts))/time.dt)+1;
time.dT_UTC = 1;                                                            % time difference relative to UTC (hours)          

% headwall dem parameters
H.cs = 10 ;                                                                   % cell size of DEM
H.x_bounds = [396150 399700];                                                 % UTM easting limits of zone to compute grav. transport for
H.y_bounds = [7533150 7535550];                                               % UTM northing bounds of zone to compute grav. transport for

% Avalanche parameters beta max and d max 
H.Blim =30 ;       % [degrees]                                                              % see willibald et al. 2020 for angle of repose
H.Dlim = .05 ;       % [ m w e ]                                                              % [ m we ] % to be est.: look at gruber: should be no prob to use a large
                                                                                           %  value if you use the max runnout inclusion
H.runout_a = 5 ;   % [degrees]                                                              % runnout angle: this is the lowest angle at which deposits occur, "headwall" only considers terrain above this angle 

% DEM to be used 
H.D_elev_m = 'tarfala.mat' ;%

% constants
C.precipZ = 0.3 ; 
C.rainsnowT = 272.16 ; 

% climate data: 
% should contain precip, wind speed, air temp at temp. resolution below timestep  
load(['wxdata.mat']);

% intialise video hold
Acc_tot = zeros(1085,1261); 

% write video of avalanche progress?
write_video = 'yes' ; % or 'no'

%% GRID 
[H,DEM] = INIT_GRUBER_headwall(H);

%% WIND 
[WD] = INIT_GRUBER_wind_drift(); 

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
var_to_plot = totalAcc ; % can be modified to ex Fnb(:,:,3) to visualise inter. steps
z1=[min(DEM.z):20:max(DEM.z)]; % contour width
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

%% optional: write video 


if strcmp('yes', write_video) 
Acc_tot = Acc_tot + H.Acc ;
figure
set(gcf,'Visible', 'off');
hold on 
var_to_plot = Acc_tot ;
z1=[900:20:2100];
var_to_plot(var_to_plot==0) = NaN ; 
box on, 
contour(DEM.x,DEM.y,DEM.z,z1,'LineWidth',0.5,'color',[0,0,0]+0.5),hold on 
S = pcolor(DEM.x,DEM.y,var_to_plot); shading flat;
cRange = S.CData ; 
alpha 0.5
cmp = (abs(cbrewer('seq','Blues',128,'spline')));
cmp(cmp>1)=1;
colormap(cmp)
c = colorbar; 
%caxis([0 0.5])
caxis([min(cRange(:)),max(cRange(:))]);
ylabel(c,'Avalanched Mass');
axis equal
xlim([396150 399700]),ylim([7533150 7535550])
xlabel('UTM easting (m)'),ylabel('UTM northing (m)')
set(gca, 'Box', 'on','LineWidth',1);
ax = gca ; 
ax.XColor = 'k';
ax.YColor = 'k';
hold off
F = getframe ;
Frame(t)= F;
    
% make video
% create the video writer with 1 fps
  writerObj = VideoWriter('myVideo.avi');
  writerObj.FrameRate = 10;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(Frame)
    % convert the image to a frame
    frame = Frame(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
end 

function [WD]=INIT_WINSTRAL(DEM,WD,H)

disp('Calculating Sheltering Index...')

% add windshelteing file path
addpath([pwd '\WindSheltering']) 

% prepare file location
 mkdir([pwd '\WindSheltering'], [num2str(WD.Sdistance) 'm']);

% fake an R 
R = "no R used" ; 
% Dominant wind direction to calculate. Given in degrees, North = 0.
windDirs = WD.windDirs;

% Sheltering distance given in meters. 
sDist = WD.Sdistance; 

% ------- THESE PARAMETERS DO NOT NORMALLY NEED TO BE CHANGED -------------

% Slope break distance in meters
rDist = sDist + 75 ; % used to be +  75;

% Maximium distance for calculating slope breaks given in meters
rDistMax = 1000; % used to be 1000

% -------------------------------------------------------------------------

% Read in the DEM
Z = flipud(DEM.z) ; 

cx = H.cs ; 
cy = H.cs ; 

% Loop over all wind directions
for windDir = 1:length(windDirs)
    % Run the Winstral model ( Di = DD ) 
    [Sx,Sb,Si,So,Di] = winstral(Z,cx,cy, windDirs(windDir), ...
                                            sDist, rDist, rDistMax); 
    
    % reset the nans to zero 
    Sx(isnan(Sx))=0; 
    
    Sx(500:535,850:860) = Sx(500:535,850:860)+20 ; 
                                        
    % Save files
    save([pwd '\WindSheltering\' num2str(WD.Sdistance) 'm\sheltering_dir_' num2str(windDirs(windDir)) '.mat'],...
             'Sx','Sb','Si','So','Di','windDir', 'sDist', 'rDist', ...
             'rDistMax','Z','R');

end


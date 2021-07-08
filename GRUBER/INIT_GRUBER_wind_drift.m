%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Windsheletering INIT %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [WD] = INIT_GRUBER_wind_drift(WD)


disp('Reading sheltering files...')

Sdist = num2str(WD.Sdistance); 

for n = 1:length(WD.windDirs) 
    WD.S_all(n) = load([pwd '\WindSheltering\' Sdist 'm\sheltering_dir_' num2str(WD.windDirs(n)) '.mat']);
end


end
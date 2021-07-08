%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% precip slope break: By Rickard Petterson (modified to fit in EBFM)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% utilises winstral et al.'s sheltering index to recalculate snow on the
% ground 

function [A_tot] = GRUBER_precipslopebreak(CurP, DEM, wx,time,WD,C)

% AWS location: 
UTMx = 399733.53 ; % coordinates of the AWS location where windspeed is measured
UTMy = 7535024.61 ; 

% Wind directions that have been calculated for the sheltering index
%windDirs = [0.1:10:350.1];

% ------------------------------------------------------------------------
% Read in the AWS values
timeW = time.TCUR_DT ;  

x = ones(length(timeW),1).* UTMx ; % UTM easting = 399733.53
y = ones(length(timeW),1).* UTMy ; % UTM northing = 7535024.61      % UTM zone = 34 W 
p = CurP ;                          % current Precip
 
wind_dir = wx.wind_hourly ; 
%v = wind_dir.Vindriktning(wind_dir.time == time.TCUR_DT); % wind direction at moment
v = nanmedian(wind_dir.Vindriktning((find(wind_dir.time == time.TCUR_DT)-...
        (time.dt*24-1)):find(wind_dir.time == time.TCUR_DT))); % we use the median as main direction 
v(isnan(v))=nanmean(wind_dir.Vindriktning) ;              % of timestep (no averaging)

%%

% Find all unique times to process
timestep = unique(timeW);

% Load the specific sheltering index file
% for n = 1:length(windDirs) 
%     % S_all(n) = load(['sheltering_dir_' num2str(windDirs(n)) '_fix.mat']);
%     S_all(n) = load([io.homedir '\WindSheltering\sheltering_dir_' num2str(windDirs(n)) '.mat']);
% end

% Initialize output
A_tot = [];

%tic
for n = 1:length(timestep)
    %disp(['Date: ' datestr(timestep(n) + datenum(1899,12,30))]);
    
    % Find all entries of AWS stations for this timestep
    i = find(timeW == timestep(n));
    
    % Determine dominant wind direction
    vm = mean(v(i), 'omitnan');
    
    % Determine the dominant wind direction
    [~,j] = min(abs(vm - WD.windDirs));
   
    % Get the specific sheltering index file
    S = WD.S_all(j);
    
     % Rescale the sheltering index
    if verLessThan('matlab','9.3')
        Sx = rescaleSX(S.Sx,abs(min(S.Sx(:))/max(S.Sx(:))).*-1,1);
    else
        Sx = rescale(S.Sx,abs(min(S.Sx(:))/max(S.Sx(:))).*-1,1);
    end
    Sx(isnan(Sx)) = 0;
    S.Di(isnan(S.Di)) = 0;
    
    % Adjust for slope break
    Sx = Sx .* S.Di; 
    
    % Set right side up 
    Sx =flipud(Sx);
       
    % use the dem to spread out over elevation     
    P = p(i) .* (1 + C.precipZ .* (DEM.z./100));

    % Remove any negative precipitation
    P(P < 0) = 0;
    %Remove NaNs
    P(isnan(P)) = 0 ; 
    
    % Rescale the variability so it can be added to the sheltering index
    if max(P(:)) > 0
        if verLessThan('matlab','9.3')
            Ps = rescaleSX(P,min(P(:))/max(P(:)),1);
        else
            Ps = rescale(P,min(P(:))/max(P(:)),1);
        end
    else
        Ps = zeros(size(P));
    end
    
    %re-remove NaNs
    Ps(isnan(Ps))=0 ;
    
    % Add the rescaled spatial precipitation variability
    % forming an accumulation factor
    Af = Sx + Ps;
    
    % Negative accumulation factor means erosion, 
    % so set it to zero
    Af(Af < 0) = 0;
    
    % Normalize the accumulation factor to the total 
    % so the sum will equal 1.
    Af = Af./sum(Af(:));

    % Calculate the total precipitation volume
    Pv = sum(P(:));

    % Multiply with precipitation volume to get distributed precipitation
    A = Af .* Pv;
       
    % Accumulate the precipitation
    if (isempty(A_tot))
        A_tot = zeros(size(A));
    end
    
   
    % Accumulate precipitation
    A_tot = A_tot + A;
      
end
%toc 

% Write output files
% save('Accumulation.mat','-v7.3');   

end 

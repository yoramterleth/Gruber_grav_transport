%% short function to produce precip and temp for testing of avalanches 
function [H] = GRUBER_headwall_climate(H,time,C,WD,DEM,wx) 


%% air temp (K)
                                                              
% INPUT = temp at Nikka station (K)                                    
T_lapse_rate = -0.007 ;                                  % Temperature lapse rate (K m-1) (empirically derived)
NIKKA_2d_z = (DEM.z - wx.info.nikka.alt) ;                      % wx station - dem cells elev diff


% use the mean of the previous timestep
CurT = nanmean(wx.NIKKA_temp_hourly.Lufttemperatur...                   % identify the temperature at that moment: might need a
    (find(wx.NIKKA_temp_hourly.time == time.TCUR_DT)-(time.dt*24-1):find(wx.NIKKA_temp_hourly.time == time.TCUR_DT)) + 273.15);


% nan values due to lack of data get replaced with average temp
if isnan(CurT)
    disp('replaced NaN in temperature on this timestep!')
    CurT = nanmean(wx.NIKKA_temp_hourly.Lufttemperatur + 273.15) ;
end

H.T = CurT + T_lapse_rate .* NIKKA_2d_z ;                                 % apply over elevation cells

%% wind speed (m/s)

% wind speed average over prvious timestep
CurWS = nanmean(wx.wind_hourly.Vindhastighet((find(wx.wind_hourly.time == time.TCUR_DT)-...
        (time.dt*24-1)):find(wx.wind_hourly.time == time.TCUR_DT)));

% nan values due to lack of data get replaced with random WS between 0-10
% m/s
if isnan(CurWS)
    disp('replaced NaN in wind speed on this timstep!')
    max_WS = 10;   
    CurWS = max_WS*rand;
end

%% precip

% tranform mm to m we 
precip_mwe = wx.NIKKA_precip_hourly.Nederbrdsmngd./1000; % vs/1000 

% use sum of previous hours in the timestep -1 (to avoid
    % overlap)
CurP = nansum(precip_mwe((find(wx.NIKKA_precip_hourly.time == time.TCUR_DT)-(time.dt*24-1)): find(wx.NIKKA_precip_hourly.time == time.TCUR_DT)));

if isnan(CurP) % rplace with 0 precip if missing data
    disp('replaced NaN in precip on this timestep!')
    CurP = 0; 
end 
clear precip_mwe

% wind sheltering: (Winstral 2002)  
if  max(H.T(:))< C.rainsnowT && CurWS >= WD.thresholdWS                   % turned on ? % all precip falls as snow? Wind speed above 0 ?
   [CurA]=GRUBER_precipslopebreak(CurP,DEM,wx,time,WD,C);
   H.P = CurA; 
else 
    H.P = CurP .* (1 + C.precipZ .* NIKKA_2d_z/100); 
end 

end 

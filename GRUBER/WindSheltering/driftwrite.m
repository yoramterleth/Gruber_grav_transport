%% write drift zones to a tiff

geotiffwrite(['DriftZones.tif'],Di,S.R, 'CoordRefSysCode', 'EPSG:3006');
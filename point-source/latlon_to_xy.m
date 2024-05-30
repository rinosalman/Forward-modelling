function [x, y] = latlon_to_xy (lat, lon, lat0, lon0)

%% Convert lat-lon-height to local xy via a given origin
%% Based on code from PCAIM (Kositsky/Perfettini)
%% Edited by Emma, 19 Mar 2011
%%
%% Inputs list of lat/lon coords and origin coords

Nstn = length(lat);

%% Convert to decimal seconds
lat = 3600 * lat;
lon = 3600 * lon;
lat0 = 3600 * lat0;

diffLon = 3600 * lon0*ones(size(lon)) - lon;

%% Make lat or lon values that are exactly zero a very small number, otherwise it gives NaN
lat(find(lat==0)) = 0.000000001;
lon(find(lon==0)) = 0.000000001;

%% Initialize output
xy = zeros(Nstn,2);

%% Loop through and do projection
for k = 1:Nstn
     xy(k,:) = polyconic(lat(k), diffLon(k), lat0);
end

%% Convert units from m into kilometer and flip x-axis
x = -xy(:,1) / 1000.0;
y =  xy(:,2) / 1000.0;





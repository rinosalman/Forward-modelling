%% Generate synthetic surface deformation due to Java megathrust locking
%% July, 2024
%% Rino Salman, EOS-RS Lab, NTU
close all;clear all;clc

%% targeted sites
fid1 = fopen('gps_subsidence_susilo+2023','rt');
readdata = textscan(fid1,'%f%f%f%s','HeaderLines',1);
[lonSta,latSta,veloObsUp,site] = readdata{1:4};
fclose(fid1);


%% fault parameters
fid2 = fopen('fault_parameters.txt','rt');
readdata = textscan(fid2,'%f%f%f%f%f%f%f%f%f','HeaderLines',1);
[lonTopLeft,latTopLeft,depthTop,length,width,str,dip,slip,rake] = readdata{1:9};
fclose(fid2);

%% C O M P U T A T I O N %%
% origin
lon0 = mean(lonSta);
lat0 = mean(latSta);

% convert geographic coordinate to cartesian
[xs, ys] = latlon_to_xy (latSta, lonSta(:), lat0, lon0);
loc = [xs, ys]; 

% convert geographic coordinate of fault plane to cartesian (unit m)
[xF1, yF1] = latlon_to_xy (latTopLeft, lonTopLeft, lat0, lon0);

% synthetic surface displacement (unit m)
mInitial1 = [xF1,yF1,depthTop,str,dip,rake,slip,length,width];
syn_enu = zeros(numel(lonSta),3);
for j=1:size(mInitial1,1)
    [enu] = los_greens_function(mInitial1(j,:),loc);
    temp =  [enu(:,1), enu(:,2), enu(:,3)];
    syn_enu = syn_enu+temp;
end

% make the unit mm
syn_enu = syn_enu.*1000;

% combine with coordinate
lonlat_synenu = [lonSta,latSta,syn_enu,veloObsUp];
lonlat_synenu_site = [num2cell(lonlat_synenu),site];

% plot to check
figure
quiver(lonlat_synenu(:,1),lonlat_synenu(:,2),lonlat_synenu(:,3),lonlat_synenu(:,4))
title('Synthetic horizontal displacement')
figure
quiver(lonlat_synenu(:,1),lonlat_synenu(:,2),zeros(numel(lonlat_synenu(:,3)),1),lonlat_synenu(:,5));hold on
title('Synthetic vertical displacement')


% save the synthetic displacement
mytitle={'#lon','lat','syntheticVeloEastMMPERYEAR','syntheticVeloNorthMMPERYEAR','syntheticVeloUpMMPERYEAR','obsVeloUpMMPERYEAR','siteName'};
comb=[mytitle;lonlat_synenu_site];
writecell(comb,'synthetic_enu_obs_ver.txt','Delimiter','space')
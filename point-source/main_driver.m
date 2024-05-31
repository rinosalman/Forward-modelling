%% Generate synthetic surface deformation of an earthquake
%% August, 2022
%% add lines for plotting, June 2024
%% Rino Salman, EOS, NTU
close all;clear all


%% Provide required information here
% range of a region to generate synthetic locations
S = 0;
N = 1;
W = 95;
E = 99;
incEW = 0.0625;
incNS = 0.0225;

% fault plane solutions (either from GCMT or USGS)
% strike, dip, rake (degrees)
plane1 = [ 0,90,0];
plane2 = [90,90,0];

% hypocenter location
hypoLon = 97;
hypoLat = 0.5;
hypoDep = 1000; %meter

% magnitude
mw = 6.5;

% earthquake type
eqType = 'strike-slip'; %options: 'normal', 'thrust', 'strike-slip'
laws = 'Wells_and_Coppersmith'; %options: 'Wells_and_Coppersmith', 'Blaser_et_all'



%% C O M P U T A T I O N %%
% generate the synthetic locations
[extraLon,extraLat] = meshgrid(W:incEW:E,S:incNS:N);

% origin
lon0 = mean(extraLon(:));
lat0 = mean(extraLat(:));

% convert geographic coordinate to cartesian
[xs, ys] = latlon_to_xy (extraLat(:), extraLon(:), lat0, lon0);
loc = [xs*1000, ys*1000]; % to meter

% unit vectors in east, north, and up components
% ALOS-2
%uvENUascA2 = [-0.510910046564734, -0.0999752952345185, 0.853800766695691];
%uvENUdscA2 = [ 0.621778720611375, -0.11081923838962, 0.775312601643633];
% Sentinel-1
uvENUascS1 = [-0.591690104047, -0.127004014368, 0.796098398603];
uvENUdscS1 = [ 0.63818047602, -0.137820596558, 0.75744938472];

% unit vectors
Ndata=numel(extraLon(:));
uvEasc = ones(Ndata,1).*uvENUascS1(1);
uvNasc = ones(Ndata,1).*uvENUascS1(2);
uvUasc = ones(Ndata,1).*uvENUascS1(3);
uvEdsc = ones(Ndata,1).*uvENUdscS1(1);
uvNdsc = ones(Ndata,1).*uvENUdscS1(2);
uvUdsc = ones(Ndata,1).*uvENUdscS1(3);

% get length, width, and expected slip
[L, Wi, slip] = rupture_size_slip(eqType, mw, laws);

% generate fault planes coordinates
% plane 1
% ulc = upper left corner
% trc = top right corner
ulcLonF1 = hypoLon + km2deg(0.5*L/1000)*sind(plane1(1)+180);
ulcLatF1 = hypoLat + km2deg(0.5*L/1000)*cosd(plane1(1)+180);
trcLonF1 = hypoLon + km2deg(0.5*L/1000)*sind(plane1(1));
trcLatF1 = hypoLat + km2deg(0.5*L/1000)*cosd(plane1(1));
% plane 2
ulcLonF2 = hypoLon + km2deg(0.5*L/1000)*sind(plane2(1)+180);
ulcLatF2 = hypoLat + km2deg(0.5*L/1000)*cosd(plane2(1)+180);
trcLonF2 = hypoLon + km2deg(0.5*L/1000)*sind(plane2(1));
trcLatF2 = hypoLat + km2deg(0.5*L/1000)*cosd(plane2(1));

% convert geographic coordinate of fault planes to cartesian (unit m)
% plane 1
[xF1, yF1] = latlon_to_xy (ulcLatF1, ulcLonF1, lat0, lon0);
xF1 = xF1 *1000;
yF1 = yF1 *1000;
% plane 2
[xF2, yF2] = latlon_to_xy (ulcLatF2, ulcLonF2, lat0, lon0);
xF2 = xF2 *1000;
yF2 = yF2 *1000;

% synthetic surface displacement due to plane 1 (unit m)
mInitial1 = [xF1,yF1,hypoDep,plane1(1),plane1(2),plane1(3),slip,L,Wi];
[synENU1] = los_greens_function(mInitial1,loc);
temp = [extraLon(:), extraLat(:), synENU1(:,1), synENU1(:,2), synENU1(:,3)];
writematrix(temp,'synthetic_enu_due_to_plane1.txt','Delimiter','space')

% convert ENU to LOS
% ascending
losAsc1=[];
for i=1:length(uvEasc)
    temp = uvEasc(i) * synENU1(i,1) + uvNasc(i) * synENU1(i,2) + uvUasc(i) * synENU1(i,3);
    losAsc1 = [losAsc1;temp];
end
%save
temp1 = [extraLon(:), extraLat(:), uvEasc, uvNasc, uvUasc, losAsc1];
writematrix(temp1,'synthetic_los_S1_due_to_plane1_ascending.txt','Delimiter','space')

% descending
losDsc1=[];
for i=1:length(uvEdsc)
    temp = uvEdsc(i) * synENU1(i,1) + uvNdsc(i) * synENU1(i,2) + uvUdsc(i) * synENU1(i,3);
    losDsc1 = [losDsc1;temp];
end
%save
temp2 = [extraLon(:), extraLat(:), uvEdsc, uvNdsc, uvUdsc, losDsc1];
writematrix(temp2,'synthetic_los_S1_due_to_plane1_descending.txt','Delimiter','space')

% synthetic surface displacement due to plane 2 (unit m)
mInitial2 = [xF2,yF2,hypoDep,plane2(1),plane2(2),plane2(3),slip,L,Wi];
[synENU2] = los_greens_function(mInitial2,loc);
temp = [extraLon(:), extraLat(:), synENU2(:,1), synENU2(:,2), synENU2(:,3)];
writematrix(temp,'synthetic_enu_due_to_plane2.txt','Delimiter','space')

% convert ENU to LOS
% ascending
losAsc2=[];
for i=1:length(uvEasc)
    temp = uvEasc(i) * synENU2(i,1) + uvNasc(i) * synENU2(i,2) + uvUasc(i) * synENU2(i,3);
    losAsc2 = [losAsc2;temp];
end
%save
temp3 = [extraLon(:), extraLat(:), uvEasc, uvNasc, uvUasc, losAsc2];
writematrix(temp3,'synthetic_los_S1_due_to_plane2_ascending.txt','Delimiter','space')

% descending
losDsc2=[];
for i=1:length(uvEdsc)
    temp = uvEdsc(i) * synENU2(i,1) + uvNdsc(i) * synENU2(i,2) + uvUdsc(i) * synENU2(i,3);
    losDsc2 = [losDsc2;temp];
end
%save
temp4 = [extraLon(:), extraLat(:), uvEdsc, uvNdsc, uvUdsc, losDsc2];
writematrix(temp4,'synthetic_los_S1_due_to_plane2_descending.txt','Delimiter','space')


%% plot
% unwrapped LOS
figure
sz=3;
tocm=100;
t=tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
nexttile
scatter(temp1(:,1),temp1(:,2),sz,temp1(:,6)*tocm,'filled')
colormap(redblue(128))
box on
nexttile
scatter(temp2(:,1),temp2(:,2),sz,temp2(:,6)*tocm,'filled')
colormap(redblue(128))
box on
nexttile
scatter(temp3(:,1),temp3(:,2),sz,temp3(:,6)*tocm,'filled')
colormap(redblue(128))
box on
nexttile
scatter(temp4(:,1),temp4(:,2),sz,temp4(:,6)*tocm,'filled')
colormap(redblue(128))
box on
h = axes(t,'visible','off');
c = colorbar(h,'Position',[0.965 0.330 0.022 0.35]);

% wrapped LOS
figure
%wvlA2=24; % wavelength of ALOS-2, in cm
wvlS1=5.6; % wavelength of sentinel-1, in cm
temp1(:,6)=mod(temp1(:,6)*tocm,wvlS1/2);
temp2(:,6)=mod(temp2(:,6)*tocm,wvlS1/2);
temp3(:,6)=mod(temp3(:,6)*tocm,wvlS1/2);
temp4(:,6)=mod(temp4(:,6)*tocm,wvlS1/2);
t=tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
nexttile
scatter(temp1(:,1),temp1(:,2),sz,temp1(:,6),'filled')
colormap(redblue(128))
box on
nexttile
scatter(temp2(:,1),temp2(:,2),sz,temp2(:,6),'filled')
colormap(redblue(128))
box on
nexttile
scatter(temp3(:,1),temp3(:,2),sz,temp3(:,6),'filled')
colormap(redblue(128))
box on
nexttile
scatter(temp4(:,1),temp4(:,2),sz,temp4(:,6),'filled')
colormap(redblue(128))
box on

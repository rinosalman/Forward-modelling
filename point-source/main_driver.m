%% Generate synthetic surface deformation due to an earthquake provided by user
%% August, 2022
%% Rino Salman, EOS, NTU
close all;clear all


%% Provide required information here
% region range to generate synthetic locations
S = 0;
N = 1;
W = 95;
E = 99;
incEW = 0.25;
incNS = 0.25;

% fault plane solutions (either from GCMT or USGS)
% strike, dip, rake (degrees)
plane1 = [ 0,90,0];
plane2 = [90,90,0];

% hypocenter location (either from GCMT or USGS)
hypoLon = 97;
hypoLat = 0.5;
hypoDep = 0; %meter

% magnitude (either from GCMT or USGS)
mw = 6.5;

% earthquake type
eqType = 'strike-slip'; %options: 'normal', 'thrust', 'strike-slip'
laws = 'Wells_and_Coppersmith'; %options: 'Wells_and_Coppersmith', 'Blaser_et_all'

% unit vectors in east, north, and up components
% only if you intend to generate LOS deformation and you know the unit vectors for your areas, currently I don't know how to generalize the unit vectors that would work for any areas
uvENUasc = [-0.591690104047, -0.127004014368, 0.796098398603]; 
uvENUdsc = [ 0.63818047602, -0.137820596558, 0.75744938472]; 




%% C O M P U T A T I O N %%
% generate the synthetic locations
[extraLon,extraLat] = meshgrid(W:incEW:E,S:incNS:N);

% origin
lon0 = mean(extraLon(:));
lat0 = mean(extraLat(:));

% convert geographic coordinate to cartesian
[xs, ys] = latlon_to_xy (extraLat(:), extraLon(:), lat0, lon0);
loc = [xs*1000, ys*1000]; % to meter

% unit vectors
Ndata=numel(extraLon(:));
uvEasc = ones(Ndata,1).*uvENUasc(1);
uvNasc = ones(Ndata,1).*uvENUasc(2);
uvUasc = ones(Ndata,1).*uvENUasc(3);
uvEdsc = ones(Ndata,1).*uvENUdsc(1);
uvNdsc = ones(Ndata,1).*uvENUdsc(2);
uvUdsc = ones(Ndata,1).*uvENUdsc(3);

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
temp = [extraLon(:), extraLat(:), uvEasc, uvNasc, uvUasc, losAsc1];
writematrix(temp,'synthetic_los_due_to_plane1_ascending.txt','Delimiter','space')

% descending
losDsc1=[];
for i=1:length(uvEdsc)
    temp = uvEdsc(i) * synENU1(i,1) + uvNdsc(i) * synENU1(i,2) + uvUdsc(i) * synENU1(i,3);
    losDsc1 = [losDsc1;temp];
end
%save
temp = [extraLon(:), extraLat(:), uvEdsc, uvNdsc, uvUdsc, losDsc1];
writematrix(temp,'synthetic_los_due_to_plane1_descending.txt','Delimiter','space')

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
temp = [extraLon(:), extraLat(:), uvEasc, uvNasc, uvUasc, losAsc2];
writematrix(temp,'synthetic_los_due_to_plane2_ascending.txt','Delimiter','space')

% descending
losDsc2=[];
for i=1:length(uvEdsc)
    temp = uvEdsc(i) * synENU2(i,1) + uvNdsc(i) * synENU2(i,2) + uvUdsc(i) * synENU2(i,3);
    losDsc2 = [losDsc2;temp];
end
%save
temp = [extraLon(:), extraLat(:), uvEdsc, uvNdsc, uvUdsc, losDsc2];
writematrix(temp,'synthetic_los_due_to_plane2_descending.txt','Delimiter','space')


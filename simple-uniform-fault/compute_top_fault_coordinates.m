%% compute top fault mid coordinates
format long g

% load coordinates drawn from Google Earth
data=load('fault_top_coordinates.gmt');
N=size(data,1);
Nm1=N-1;
temp=[];

% defined common parameters (units: degree and meter)
DepTop=6000;
faultWidth=200000;
faultDip=16;
slip=0.07;
rake=-90;

% compute
temp1=[];
for i=1:Nm1
    [ARCLEN,az]=distance(data(i,2),data(i,1),data(i+1,2),data(i+1,1));
    distinkm=deg2km(ARCLEN)*1000;%to meter
    lonTopLeft=data(i,1);
    latTopLeft=data(i,2);
    faultLength=distinkm;
    temp2=[lonTopLeft,latTopLeft,DepTop,faultLength,faultWidth,az,faultDip,slip,rake];
    temp1=[temp1;temp2];
end

% save
mytitle={'#faultlonTopLeft','faultlatTopLef','faultDepthTop','faultLength','faultWidth','faultStrike','faultDip','slip','rake'};
comb=[mytitle;num2cell(temp1)];
writecell(comb,'fault_parameters.txt','Delimiter','space')
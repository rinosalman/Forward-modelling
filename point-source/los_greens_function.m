function [enuDisp] = los_greens_function(model_parameters,loc)

%% extract model parameters
xfault = model_parameters(1);
yfault = model_parameters(2);
depth  = model_parameters(3);
strike = model_parameters(4);
dip    = model_parameters(5);
rake   = model_parameters(6);
slip   = model_parameters(7);
fault_length = model_parameters(8);
fault_width = model_parameters(9);

%% convert from degree to radian
strike = deg2rad(strike);
dip = deg2rad(dip);
rake = deg2rad(rake);

%% Split strike slip and dip slip parts.
strike_slip = slip*cos(rake);
dip_slip    = slip*sin(rake);

%% fault dimension (in meters) 
fault_depth  = depth;

%% distance from top left fault to observations points
xd = loc(:,1) - xfault;
yd = loc(:,2) - yfault;

%% run okada
nu = 0.25;
[enudisp_strike_slip(:,1), enudisp_strike_slip(:,2), enudisp_strike_slip(:,3)] = computeDisplacementOkada85(strike_slip,xd(:),yd(:),nu,dip,fault_depth,fault_length,fault_width,1,strike);
[enudisp_dip_slip(:,1), enudisp_dip_slip(:,2), enudisp_dip_slip(:,3)]          = computeDisplacementOkada85(dip_slip   ,xd(:),yd(:),nu,dip,fault_depth,fault_length,fault_width,2,strike);

%% compute total motion 
enuDisp = enudisp_strike_slip + enudisp_dip_slip; 

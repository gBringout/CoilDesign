%clear all;
%close all;

% Coding rules :
% Function start with a Capital letter
% variable start with a lower-case letter
%    but : only the MATRIX used for matrix calculation are completly named in CAPITAL


%% Import meshing
% If reduction = 1, the matrix system will be reduced, in order to set all
% the border node to the same value, in order to respect the divergence free criteria
% i.e. the node vector will be reorgenized in order to have all the border
% node on top, in order to facilitate the reduction of all the system

[coil.listNode,coil.listTriangle,coil.tri] = importMeshWavefront(fullfile('.','data','30x13_R145_D314.obj'));
coil.center = [0 0 0];
coil.reduction = 1;
coil.rateIncreasingWire = 1;

%% Target point definition

degreeMax = 15;
orderMax = 15;
rhoReference = 0.05/2; % radius of reference
rk = createTargetPointGaussLegendreAndRectangle7(rhoReference,degreeMax,orderMax);
%% Then we habe to calculate the field in a given direction

% Initialize the target ampltiude
bc(1).coefficient = zeros(degreeMax+1,orderMax+1);
bs(1).coefficient = zeros(degreeMax+1,orderMax+1);
bc(2).coefficient = zeros(degreeMax+1,orderMax+1);
bs(2).coefficient = zeros(degreeMax+1,orderMax+1);
bc(3).coefficient = zeros(degreeMax+1,orderMax+1);
bs(3).coefficient = zeros(degreeMax+1,orderMax+1);

% define the target ampltiude
targetCoil = 'DriveZ';
bc(3).coefficient(1,1) = 0.003; % MPI - Drive Z - Homogeneous Field in Z Direction

B  = RebuildField7bis(bc,bs,rhoReference, rk,'sch');
coil.btarget = [B(1,:) B(2,:) B(3,:)];

%DisplayFieldOnSphere( B,rk,'TargetField' )

clear('B');
 
 %% In order to calculate the resistance matrix, we have to providfe some data :

% For the coil
coil.wireThickness = 0.00395; % (meter) Thickness of the conductor
coil.wireWidth = 0.00395; % (meter) Thickness of the conductor
coil.wireSurface = coil.wireThickness*coil.wireWidth; % in meter %5mmx5mm is equivalent to the number used in Timo's coil or the 7.5*7.5 litz wire
coil.fillFactor = 1;
coil.rhoCopper = 1.68*10^-8; % (Ohm*m) resistivity of the copper
coil.wireResistivity = coil.rhoCopper/coil.fillFactor;  % (Ohm*m) resistivity of the wire

%% Part to calculate
%optimizationType = 'QP';
%optimizationType = 'standardTikhonov';
optimizationType = 'generalizedTikhonov';
reg=4.0*10^-8;
calculateR = 0;
calculateL = 0;
calculateLwp = 1;

coil.startingWireNumber = 6;
coil.distanceBetween2Wire = 50*10^-3;
coil.rateIncreasingWire = 2;

%% Set the variable for the Field calculation

% X position in meter
coil.x_Start =-0.120;%0.136
coil.x_Stop = 0.120;
coil.x_Step = 0.01;
coil.x_Value = coil.x_Start:coil.x_Step:coil.x_Stop;
% Y position in meter
coil.y_Start = -0.120;
coil.y_Stop = 0.120;
coil.y_Step = 0.01;
coil.y_Value = coil.y_Start:coil.y_Step:coil.y_Stop;
% Z position in meter
coil.z_Start =-0.002;
coil.z_Stop = 0.002;
coil.z_Step = 0.002;
coil.z_Value = coil.z_Start:coil.z_Step:coil.z_Stop;

coil.current = 1;
coil.sphere_radius = 0.025;
coil.coil_radius = 0.317/2;
coil.coil_length = 0.49;
%clear all;
%close all;

% Coding rules :
% Function start with a Capital letter
% variable start with a lower-case letter
%    but : only the MATRIX used for matrix calculation are completly named in CAPITAL


%% Import meshing
% If reduction = 1, the matrix system will be reduced, in order to set all
% the border node to the same value, in order to respect the divergence free criteria
% i.e. the node vector will be reorganised in order to have all the border
% node on top, in order to facilitate the reduction of all the system
coil.reduction = 1;

[coil.listNode,coil.listTriangle,coil.tri] = importMeshWavefront(fullfile('.','data','20x20_R135_H400.obj'));
coil.center = [0 0 0];

%% Target point definition

degreeMax = 5;
orderMax = 5;
rhoReference = 0.08/2; % radius of reference
rk = createTargetPointGaussLegendreAndRectangle7(rhoReference,degreeMax,orderMax);

%% Then we habe to calculate the target field

% Initialize the target ampltiude
bc(1).coefficient = zeros(degreeMax+1,orderMax+1);
bs(1).coefficient = zeros(degreeMax+1,orderMax+1);
bc(2).coefficient = zeros(degreeMax+1,orderMax+1);
bs(2).coefficient = zeros(degreeMax+1,orderMax+1);
bc(3).coefficient = zeros(degreeMax+1,orderMax+1);
bs(3).coefficient = zeros(degreeMax+1,orderMax+1);

% Drive Y
targetCoil = 'DriveY';
bc(2).coefficient(1,1) = 0.015; % MPI - Drive Y - Homegeneous Field in Y direction

B  = RebuildField7bis(bc,bs,rhoReference, rk,'sch');
coil.btarget = [B(1,:) B(2,:) B(3,:)];

%DisplayFieldOnSphere( B,rk,'TargetField' )

clear('B');
 
 %% In order to calculate the resistance matrix, we have to providfe some data :

% For the coil
coil.wireThickness = 0.0078; % (meter) Thickness of the conductor
coil.wireWidth = 0.0078; % (meter) Thickness of the conductor
coil.wireSurface = coil.wireThickness*coil.wireWidth; % in meter %5mmx5mm is equivalent to the number used in Timo's coil or the 7.5*7.5 litz wire
coil.fillFactor = 0.5;
coil.rhoCopper = 1.68*10^-8; % (Ohm*m) resistivity of the copper
coil.wireResistivity = coil.rhoCopper/coil.fillFactor;  % (Ohm*m) resistivity of the wire

%% Information required for the optimisation
%optimizationType = 'QP';
optimizationType = 'standardTikhonov';
reg=4.0*10^-6;
calculateR = 0;
calculateL = 0;
calculateLwp = 0;

%% Details for the wire discretisation

coil.startingWireNumber = 24;
coil.rateIncreasingWire = 2;
coil.distanceBetween2Wire = 30*10^-3; % 2 wire of the same stream function level

%% Set the variable for the Field calculation

coil.current = 1;
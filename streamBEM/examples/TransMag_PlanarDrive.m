%clear all;
%close all;

% Coding rules :
% Function start with a Capital letter
% variable start with a lower-case letter
%    but : only the MATRIX used for matrix calculation are completly named in CAPITAL


%% Import meshing
% If reduction = 1, the matrix system will be reduced, in order to set all
% the border node to the same value, in order to respect the xxx criteria
% i.e. the node vector will be reorgenized in order to have all the border
% node on top, in order to facilitate the reduction of all the system

%[coil.node,coil.triangle,coil.tri] = importMeshBlender('data\20x20_R100_H280.obj');
[coil.listNode,coil.listTriangle,coil.tri] = importMeshBlender('./data/30x13_R145_D314.obj');
%[coil.listNode,coil.listTriangle,coil.tri] = importMeshBlender('./data/2500_R145_D250.obj');
%[coil.listNode,coil.listTriangle,coil.tri] = importMeshBlender('../data/BEM/mesh/blender/20x20_R400_H1600_meshRandom.obj');
%[coil.node,coil.triangle,coil.tri] = importMeshBlender('data\75x28_R119_H280.obj');
coil.center = [0 0 0];
coil.reduction = 1;
coil.rateIncreasingWire = 1;

%% Target point definition
addpath('..\SphericalHarmonics\')

degreeMax = 15;
orderMax = 15;
rhoReference = 0.025; % radius of reference
rk = createTargetPointGaussLegendreAndRectangle7(rhoReference,degreeMax,orderMax);
%% Then we habe to calculate the field in a given direction

% define the target ampltiude
bc(1).coefficient = zeros(degreeMax+1,orderMax+1);
bs(1).coefficient = zeros(degreeMax+1,orderMax+1);
bc(2).coefficient = zeros(degreeMax+1,orderMax+1);
bs(2).coefficient = zeros(degreeMax+1,orderMax+1);
bc(3).coefficient = zeros(degreeMax+1,orderMax+1);
bs(3).coefficient = zeros(degreeMax+1,orderMax+1);

% Quadrupole
%bc(1).coefficient(2,2) = targetAmplitude; % Quadrupole
%bs(2).coefficient(2,2) = -targetAmplitude; % Quadrupole

% Drive
bc(3).coefficient(1,1) = 0.003; % Drive Y
%bc(1).coefficient(1,1) = 0.003; % Drive X

B  = RebuildField7bis(bc,bs,rhoReference, rk,'sch');

%targetCoil = 'dBzdx';
%targetCoil = 'dBzdy';
%targetCoil = 'dBzdz';
%coil.btarget = [B(3,:)];

%targetCoil = 'Quad';
%coil.btarget = [B(1,:) B(2,:) B(3,:)];

targetCoil = 'DriveZ';
coil.btarget = [B(1,:) B(2,:) B(3,:)];

%targetCoil = 'DriveX';
%coil.btarget = [B(1,:) B(2,:) B(3,:)];

%coil.btarget = [B(2,:)];
%coil.error = 0.1;
%DisplayFieldOnSphere( B,rk,'TargetField' )

clear('B');
 
 %% In order to calculate the resistance matrix, we have to providfe some data :

% For the coil
coil.wireThickness = 0.00395; % (meter) Thickness of the conductor
coil.wireWidth = 0.00395; % (meter) Thickness of the conductor
coil.wireSurface = coil.wireThickness*coil.wireWidth; % in meter %5mmx5mm is equivalent to the number used in Timo's coil or the 7.5*7.5 litz wire
coil.fillFactor = 0.5;
coil.rhoCopper = 1.68*10^-8; % (Ohm*m) resistivity of the copper
coil.rho = coil.rhoCopper*coil.fillFactor;
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
coil.coil_radius = 0.1;
coil.coil_length = 0.2;
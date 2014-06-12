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
[coil.listNode,coil.listTriangle,coil.tri] = importMeshBlender('./data/41x61_R100_H400.obj');
%[coil.listNode,coil.listTriangle,coil.tri] = importMeshBlender('../data/BEM/mesh/blender/20x20_R400_H1600_meshRandom.obj');
%[coil.node,coil.triangle,coil.tri] = importMeshBlender('data\75x28_R119_H280.obj');
coil.center = [0 0 0];
coil.reduction = 1;
coil.rateIncreasingWire = 1;

%% Target point definition
addpath('..\SphericalHarmonics\')

degreeMax = 5;
orderMax = 5;
rhoReference = 0.04; % radius of reference
rk = createTargetPointGaussLegendreAndRectangle7(rhoReference,degreeMax,orderMax);
%% Then we habe to calculate the field in a given direction

% define the target ampltiude
targetAmplitude = 0.5*rhoReference;
for i=1:7
    bc(1).coefficient(i,:) = zeros(1,7);
    bs(1).coefficient(i,:) = zeros(1,7);
    bc(2).coefficient(i,:) = zeros(1,7);
    bs(2).coefficient(i,:) = zeros(1,7);
    bc(3).coefficient(i,:) = zeros(1,7);
    bs(3).coefficient(i,:) = zeros(1,7);
end
%bs(2).coefficient(3,3) = 10^-5; % shimming
%bc(2).coefficient(1,1) = 15*10^-3; % Drive Y
%bs(3).coefficient(2,2) = 1*rhoReference; % MRI gradient Y
%bc(3).coefficient(2,2) = 1*rhoReference; % MRI gradient X
%bc(3).coefficient(2,1) = targetAmplitude; % MRI gradient X
% Quadrupole
%bc(1).coefficient(2,2) = targetAmplitude; % Quadrupole
%bs(2).coefficient(2,2) = -targetAmplitude; % Quadrupole

% Drive
bc(2).coefficient(1,1) = 0.015; % Drive

B  = RebuildField7bis(bc,bs,rhoReference, rk,'sch');

%targetCoil = 'dBzdx';
%targetCoil = 'dBzdy';
%targetCoil = 'dBzdz';
%coil.btarget = [B(3,:)];

%targetCoil = 'Quad';
%coil.btarget = [B(1,:) B(2,:) B(3,:)];

targetCoil = 'DriveY';
coil.btarget = [B(2,:)];
coil.btarget = [B(1,:) B(2,:) B(3,:)];

%coil.btarget = [B(2,:)];
%coil.error = 0.1;
%DisplayFieldOnSphere( B,rk,'TargetField' )

clear('B');
 
 %% In order to calculate the resistance matrix, we have to providfe some data :

% For the coil
coil.wireThickness = 0.0078; % (meter) Thickness of the conductor
coil.wireWidth = 0.0078; % (meter) Thickness of the conductor
coil.wireSurface = coil.wireThickness*coil.wireWidth; % in meter %5mmx5mm is equivalent to the number used in Timo's coil or the 7.5*7.5 litz wire
coil.fillFactor = 0.5;
coil.rhoCopper = 1.68*10^-8; % (Ohm*m) resistivity of the copper
coil.rho = coil.rhoCopper*coil.fillFactor;
coil.wireResistivity = coil.rhoCopper/coil.fillFactor;  % (Ohm*m) resistivity of the wire

%% Part to calculate
%optimizationType = 'QP';
optimizationType = 'standardTikhonov';
reg=4.0*10^-6;
calculateR = 0;
calculateL = 0;
calculateLwp = 0;

coil.startingWireNumber = 16;
coil.distanceBetween2Wire = 30*10^-3;
coil.rateIncreasingWire = 2;

%% Set the variable for the Field calculation

% X position in meter
x_Start =-0.120;%0.136
x_Stop = 0.120;
x_Step = 0.001;
x_Value = x_Start:x_Step:x_Stop;
% Y position in meter
y_Start = -0.120;
y_Stop = 0.120;
y_Step = 0.001;
y_Value = y_Start:y_Step:y_Stop;
% Z position in meter
z_Start =-0.002;
z_Stop = 0.002;
z_Step = 0.002;
z_Value = z_Start:z_Step:z_Stop;

current = 1;
sphere_radius = 0.025;
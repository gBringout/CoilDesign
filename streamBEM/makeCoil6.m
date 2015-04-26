% Coil optimization based on a BEM methods using stream functions

clear all;
close all;

%% 1. Add the relevant path
addpath(genpath(fullfile('..')))
addpath(genpath(fullfile('..','..','SphericalHarmonics')))

%% 2. Load the information specific to each the coil design and 
% add the relevant path for the optimization package
disp('2. Load the design details')

%TransMag_Quadrupole
%TransMag_Drive
%TransMag_PlanarDrive
%MRI_GY
IWMPI_Drive

if strcmp(optimizationType,'standardTikhonov') || strcmp(optimizationType,'generalizedTikhonov')
    % Please download the regularization tools of: http://www.imm.dtu.dk/~pcha/Regutools/
    % and add the folder to matlab' path:
    addpath(genpath(fullfile('..','..','regu')))
    %addpath(genpath(fullfile('..','..','..','Software','regu')))
elseif strcmp(optimizationType,'QP')
    % Please download the YALMIP TOOLBOX from
    % http://users.isy.liu.se/johanl/yalmip/
    % and add the folder to matlab' path :
    addpath(genpath(fullfile('..','..','YALMIP')))
    %addpath(genpath(fullfile('..','..','..','Software','YALMIP')))
end

tStart=tic;

%% 3. Sort the mesh nodes to have all boundaries on top of the mesh. Then
% the mesh is pre-processed to calculate general properties as
% the air of triangle, the barycenter, which triangle is connected to which node,
% and the basis function associated to each triangle. The mesh and the target
% points are displayed to be able to check them
disp('3.Process the mesh')

disp('  Re-order the node to have all the border on the top of the node vector')
if coil.reduction == 1
    [coil.listNode,coil.listTriangle,coil.subBoundaries] = NodeSorting2(coil.listNode,coil.listTriangle);
    fprintf('%i border detected on the coil \n',size(coil.subBoundaries,1))
else
    coil.subBoundaries = [];
end
 
disp('  Process the triangle')
[coil.triangle,coil.node] = processMesh(coil.listTriangle,coil.listNode,orderQuadrature);

% Calculate the basis function of the mesh.
disp('Calculate the basis function')
coil.basis = basisFunction4(coil.node, coil.triangle,coil.center,orderQuadrature);


disp('  Ploting the mesh and the target points');
figure('Name','Data verification')
subplot(3,3,1)
trimesh(coil.listTriangle,coil.listNode(:,1),coil.listNode(:,2),coil.listNode(:,3),ones(size(coil.listNode,1),1)); 
axis square
subplot(3,3,2)
plot3(rk(:,1),rk(:,2),rk(:,3),'*');
axis square

%% 4. The different matrices are calculted, when required. And displayed
disp('4.Calculate the requires matrices')
% The Laplacian

coil.Lwp = zeros(size(coil.node,2),size(coil.node,2));
if calculateLwp
    fprintf(1,'  Calculating the Laplacian Operator.\n');
    if exist('hLaplace')
        coil.Lwp = Laplacian3(coil.node, coil.triangle,hLaplace);
    else
        coil.Lwp = Laplacian3(coil.node, coil.triangle);
    end
end

subplot(3,3,3)
imagesc(log(abs(coil.Lwp)));
title('Lwp')
axis square
colormap(gray)


%The inductance matrix

coil.L = zeros(size(coil.node,2),size(coil.node,2));
if calculateL
    disp('  Calculating the Lmn matrix.');
    coil.L = Lmn10(coil.node, coil.triangle,coil.basis,orderQuadrature);
end

subplot(3,3,4)
Llog = real(log(coil.L));
imagesc(Llog);
title('coil.L (log scale)')
axis square
colormap(gray)


% The resistance matrix

coil.R = zeros(size(coil.node,2),size(coil.node,2));

if calculateR
    disp('  Calculating the Rmn matrix.');
    coil.R = Rmn9(coil.node,coil.triangle,coil.basis,coil.wireResistivity,coil.wireThickness);
end

subplot(3,3,5)
imagesc(coil.R);
title('coil.R')
axis square
colormap('gray')


% The magnetic flux density matrix

disp('  Calculating the Cn matrix.');
[coil.Cx,coil.Cy,coil.Cz] = Cn7(coil.node,coil.triangle,coil.basis,rk,orderQuadrature);
coil.C = [coil.Cx;coil.Cy;coil.Cz];

subplot(3,3,6)
imagesc(coil.C);
title('coil.C')
caxis([min(min(coil.C)) max(max(coil.C))])
colormap('gray')

% The target magnetic flux density matrix is build, accoring to 
% the type of coil which has to be designed

if strcmp(targetCoil,'dBzdx') || strcmp(targetCoil,'dBzdy') || strcmp(targetCoil,'dBzdz')
    coil.Ctarget = coil.Cz;
elseif strcmp(targetCoil,'DriveZ') || strcmp(targetCoil,'DriveY') || strcmp(targetCoil,'DriveZ')
    coil.Ctarget = [coil.Cx;coil.Cy;coil.Cz];
elseif strcmp(targetCoil, 'Quad')
    coil.Ctarget = [coil.Cx;coil.Cy;coil.Cz];
end

%% 5. The matrix are reduced to used the normalized set 
% of basis function (if required). Some are then displayed
disp('5. Normalising the set of basis function')

if coil.reduction && size(coil.subBoundaries,1)>0
    disp('  Reduction of the matrix')
    coil.nonReducedSize = size(coil.R,1);
    coil.Rfull = coil.R;
    coil.R = ReduceSquareMatrix4(coil.Rfull,coil.subBoundaries);
    coil.Lfull = coil.L;
    coil.L = ReduceSquareMatrix4(coil.Lfull,coil.subBoundaries);
    coil.CtargetFull = coil.Ctarget;
    coil.Ctarget = ReduceCMatrix4(coil.CtargetFull,coil.subBoundaries);
    coil.Cfull = coil.C;
    coil.C = ReduceCMatrix4(coil.Cfull,coil.subBoundaries);
    coil.LwpFull = coil.Lwp;
    coil.Lwp = ReduceSquareMatrix4(coil.LwpFull,coil.subBoundaries);
end

subplot(3,3,7)
Llog = real(log(coil.L));
imagesc(Llog);
title('reduced coil.L')
axis square
colormap(gray)

subplot(3,3,8)
imagesc(coil.R);
title('reduced coil.R')
axis square
colormap(gray)

%% 6. The optimization problem is solved
disp('6.Optimisation')

if strcmp(optimizationType,'standardTikhonov')
    disp('  Standard Tikhonov')
    [U,s,V] = csvd(coil.Ctarget);
    spectrum = sort(s(:,1),'descend');
    subplot(3,3,9)
    semilogy(spectrum)
    %figure; [reg_corner,rho,eta,reg_param] = l_curve(U,s,coil.btarget');
    [coil.x_lambda,~,~] = tikhonov(U,s,V,coil.btarget',reg);
    if coil.reduction
        %coil.s = retrieveCurrentVector2(coil.x_lambda,coil.subBoundaries);
        coil.s = retrieveCurrentVector3(coil.x_lambda,coil.subBoundaries,optimizationType);
    else
        coil.s = coil.x_lambda;
    end

elseif strcmp(optimizationType,'generalizedTikhonov')
    disp('  Generalized Tikhonov')
    [U,sm,X,V,W] = cgsvd(coil.Ctarget,coil.Lwp);
    spectrum = sort(sm(:,1),'descend');
    subplot(3,3,9)
    semilogy(spectrum)
    %figure; [reg_corner,rho,eta,reg_param] = l_curve(U,sm,coil.btarget');
    [coil.x_lambda,~,~] = tikhonov(U,sm,X,coil.btarget',reg);
    if coil.reduction
        coil.s = retrieveCurrentVector3(coil.x_lambda,coil.subBoundaries,optimizationType);
    else
        coil.s = coil.x_lambda;
    end
    
elseif strcmp(optimizationType,'QP')
    disp('  Quadratic Programming')
    % Optimisation parameter
    % YALMIP
    % Here we try to minimize the inductance
    x = sdpvar(size(coil.L,2),1); % define the size of the objective
    Objective = x'*coil.R*x; % set the objective
    
    % build the error constrains on the fields
    lr = zeros(1,size(coil.btarget,2));
    ur = zeros(1,size(coil.btarget,2));
    for j=1:size(coil.Ctarget,1)
        if coil.btarget(j)<0
            % lr and ur stands for "lower" and "upper"
            lr(j) = (1+coil.error)*coil.btarget(j);
            ur(j) = (1-coil.error)*coil.btarget(j);
        elseif coil.btarget(j)>0
            lr(j) = (1-coil.error)*coil.btarget(j);
            ur(j) = (1+coil.error)*coil.btarget(j);
        else
            lr(j) = -0.001;
            ur(j) = 0.001;
        end
    end
    
    % boundary of the solution
    lb = ones(size(coil.L,1),1)*-Inf;
    ub = ones(size(coil.L,1),1)*Inf;
    if coil.reduction
        for i=1:size(coil.subBoundaries,1)
            lb(i) = 0; %boundary of the mesh have to be null
            ub(i) = 0; %boundary of the mesh have to be null
        end
    else
        index = 1;
        for i=1:size(coil.subBoundaries,1)
            for j=size(coil.subBoundaries(i).node,1)
                lb(index) = 0; %boundary of the mesh have to be null
                ub(index) = 0; %boundary of the mesh have to be null
                index = index+1;
            end
        end
    end
    
    %Build the constrains object
    Constraints = [lr' <= coil.Ctarget*x <= ur',...
                   lb<=x<=ub];
               
    
    % Set some options for YALMIP and solver
    options = sdpsettings('verbose',2);
    options = sdpsettings(options,'debug',1);
    %solve the problem
    solvesdp(Constraints,Objective,options)
    
    %retrieve the solution
    coil.s_reduced = double(x);
    
    if coil.reduction
        coil.s = retrieveCurrentVector3(coil.s_reduced,coil.subBoundaries,optimizationType);
    else
        coil.s = coil.s_reduced;
    end
end


displayStreamFunction(coil.listTriangle,coil.s,coil.listNode)

%% 7. The wire centroids are aproximate and the resulting wire are displayed
disp('7. Extract the wire path')
[wire,nbrwire,minimalWireSpacing,res] = NbrWireOptimisation3(coil.listNode,coil.listTriangle,coil.s,coil.wireWidth,coil.startingWireNumber,coil.distanceBetween2Wire,coil.rateIncreasingWire);
displayWire(wire);


%% 8. Different properties are then evaluated
% Calculating the field generated by the wires
[coil.Bx,coil.By,coil.Bz] = Field3(wire,coil.current,rk(:,1)',rk(:,2)',rk(:,3)');

% making a decomposition up to degree and order 5
degreeMax = 5;
orderMax = 5;

% Displaying using the schimdt semi-normalisation
mode = 'sch';
[bc, bs] = getSphericalHarmonicsCoefficientMeasure7(coil.Bx,coil.By,coil.Bz,degreeMax,orderMax,rk,mode);
displaySHC(bc,bs,2,10^6)

% Displaying using the mathematical norme
mode = 'norm';
[bc, bs] = getSphericalHarmonicsCoefficientMeasure7(coil.Bx,coil.By,coil.Bz,degreeMax,orderMax,rk,mode);
displaySHC(bc,bs,2,10^5)

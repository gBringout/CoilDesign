% BEM methods using stream functions
clear all;
close all;
addpath(genpath(fullfile('.')))
addpath(genpath(fullfile('..','..','SphericalHarmonics')))

TransMag_Quadrupole
%TransMag_Drive
%TransMag_PlanarDrive
%MRI_GY

if strcmp(optimizationType,'standardTikhonov') || strcmp(optimizationType,'generalizedTikhonov')
    % Please download the regularization tools of: http://www.imm.dtu.dk/~pcha/Regutools/
    % and add the folder to matlab' path:
    addpath(genpath(fullfile('..','..','regu')))
    %addpath(genpath(fullfile('..','..','..','Software','regu')))
elseif strcmp(optimizationType,'QP')
    % Please download the OPTI TOOLBOX from
    % http://www.i2c2.aut.ac.nz/Wiki/OPTI/, install it
    % and add the folder to matlab' path :
    addpath(genpath(fullfile('..','..','OptiToolbox')))
    %addpath(genpath(fullfile('..','..','..','Software','OptiToolbox')))
end

tStart=tic;

%% Display the Mesh

disp('Re-order the node to have all the border on the top of the node vector')
if coil.reduction == 1
    [coil.listNode,coil.listTriangle,coil.subBoundaries] = NodeSorting2(coil.listNode,coil.listTriangle);
    fprintf('%i border detected on the coil \n',size(coil.subBoundaries,1))
else
    coil.subBoundaries = [];
end


disp('Process the triangle')
[coil.triangle,coil.node] = processMesh(coil.listTriangle,coil.listNode);

% Calculate the basis function of the mesh.
disp('Calculate the basis function')
coil.basis = basisFunction4(coil.node, coil.triangle,coil.center);
%tt = [coil.triangle(:).node];
figure('Name','Data verification')
subplot(3,3,1)
trimesh(coil.listTriangle,coil.listNode(:,1),coil.listNode(:,2),coil.listNode(:,3),ones(size(coil.listNode,1),1)); 
axis square

%% Creating the target points
disp('Ploting the target points');
subplot(3,3,2)
plot3(rk(:,1),rk(:,2),rk(:,3),'*');
axis square

%% Calculating the Laplacian

coil.Lwp = zeros(size(coil.node,2),size(coil.node,2));
if calculateLwp
    fprintf(1,'Calculating the Laplacian Operator.\n');

    coil.Lwp = Laplacian3(coil.node, coil.triangle);
end

subplot(3,3,3)
imagesc(log(abs(coil.Lwp)));
title('Lwp')
axis square
colormap(gray)
%% Attempt to calculate the Lmn matrix

coil.L = zeros(size(coil.node,2),size(coil.node,2));
if calculateL
    disp('Calculating the Lmn matrix.');
    coil.L = Lmn10(coil.node, coil.triangle,coil.basis);
end

subplot(3,3,4)
Llog = real(log(coil.L));
imagesc(Llog);
title('coil.L (log scale)')
axis square
colormap(gray)

%% Attemp to calculate the Rmn matrix

coil.R = zeros(size(coil.node,2),size(coil.node,2));

if calculateR
    disp('Calculating the Rmn matrix.');
    coil.R = Rmn9(coil.node,coil.triangle,coil.basis,coil.wireResistivity,coil.wireThickness);
end

subplot(3,3,5)
imagesc(coil.R);
title('coil.R')
axis square
colormap('gray')

%% Attemp to calculate Cn
disp('Calculating the Cn matrix.');
[coil.Cx,coil.Cy,coil.Cz] = Cn7(coil.node,coil.triangle,coil.basis,rk);
coil.C = [coil.Cx;coil.Cy;coil.Cz];

if strcmp(targetCoil,'dBzdx') || strcmp(targetCoil,'dBzdy') || strcmp(targetCoil,'dBzdz')
    coil.Ctarget = coil.Cz;
elseif strcmp(targetCoil,'DriveZ') || strcmp(targetCoil,'DriveY') || strcmp(targetCoil,'DriveZ')
    coil.Ctarget = [coil.Cx;coil.Cy;coil.Cz];
elseif strcmp(targetCoil, 'Quad')
    coil.Ctarget = [coil.Cx;coil.Cy;coil.Cz];
end

subplot(3,3,6)
imagesc(coil.C);
title('coil.C')
caxis([min(min(coil.C)) max(max(coil.C))])
colormap('gray') 

%% Reduction of the matrix to remove the border influence

if coil.reduction && size(coil.subBoundaries,1)>0
    disp('Reduction of the matrix')
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

%% Solve the optimization problem

if strcmp(optimizationType,'standardTikhonov')
    disp('Standard Tikhonov')
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
    disp('Generalized Tikhonov')
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
    disp('Quadratic Programming')
    % Optimisation parameter

    H = coil.L;
    f = zeros(size(H,1),1);                

    % Linear Constraints
    % or Linear Inequality Constraints (Ax <= b)
    A = coil.Ctarget;
    btarget = [coil.btarget];

    rl = zeros(1,size(btarget,2));
    ru = zeros(1,size(btarget,2));
    for j=1:size(A,1)
        if btarget(j)<0
            % rl and ru stands for "lower" and "upper"
            rl(j) = (1+coil.error)*btarget(j);
            ru(j) = (1-coil.error)*btarget(j);
        elseif btarget(j)>0
            rl(j) = (1-coil.error)*btarget(j);
            ru(j) = (1+coil.error)*btarget(j);
        else
            rl(j) = -0.001;
            ru(j) = 0.001;
        end
    end
    % we do not constrain the amplitude of the solution
    lb = ones(size(H,1),1)*-Inf;
    ub = ones(size(H,1),1)*Inf;

    % but we do it on the boundaries
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

    % Quadratic Constraint
    % or Quadratic Inequality (x'Qx + l'x <= r)
    %Q = {[R] [L]};             
    %l = [zeros(size(R,1),1),zeros(size(L,1),1)];
    %r = [90*10^-3;10e8];
    Q2 = coil.R;
    l = zeros(size(coil.R,1),1);
    r = 10000;


    % change the option.
    %   Solver 'ipopt' is for the linearisation (i.e. Poole method)
    %   Solver 'BONMIN' is not working
    %   Solver 'SCIP' is not wirking
    %   Solver 'Cplex' preview version is too small
    opts = optiset('display','iter',...
            'solver','ipopt',...
            'maxiter',500,...
            'maxtime',100);

    % QPQC
    Opt = opti('qp',H,f,'lin',A,rl,ru,'qc',Q2,l,r,'bounds',lb,ub,'options',opts); %bounded lineare constrains
    % QPLC
    %Opt = opti('qp',H,f,'lin',A,rl,ru,'bounds',lb,ub,'options',opts);
    %Opt = opti('qp',H,f,'ineq',A,b,'bounds',lb,ub,'qc',Q2,l,r,'options',opts);

    % Solve the QP problem
    [coil.s_reduced,fval,exitflag,info] = solve(Opt);
    if coil.reduction
        coil.s = retrieveCurrentVector3(coil.s_reduced,coil.subBoundaries,optimizationType);
    else
        coil.s = coil.s_reduced;
    end
end


displayStreamFunction(coil.listTriangle,coil.s,coil.listNode)

%% Exctract the wire based on a simple method
[wire,nbrwire,minimalWireSpacing,res] = NbrWireOptimisation3(coil.listNode,coil.listTriangle,coil.s,coil.wireWidth,coil.startingWireNumber,coil.distanceBetween2Wire,coil.rateIncreasingWire);
displayWire(wire);


%% Field properties calculation
coil.B = coil.Cfull*coil.s;
coil.p_dis = coil.s'*coil.Rfull*coil.s;
coil.e_stored = 0.5*coil.s'*coil.Lfull*coil.s;

%% Calculation of the field discritized coil
% disp('Field calculation with the discritized model')
% 
[coil.Bx,coil.By,coil.Bz] = Field2(wire,coil.current,coil.x_Value,coil.y_Value,coil.z_Value);
% Babs = sqrt(Bx.^2+By.^2+Bz.^2);
DisplayField(coil.Bx,coil.By,coil.Bz,coil.x_Value,coil.y_Value,coil.z_Value,-40*10^-6,40*10^-6);
if strcmpi(targetCoil,'Quad') || strcmpi(targetCoil,'dBzdy')
    [gradientLinearity,gradientValueCenter] = Linearity(coil.Bx,coil.By,coil.Bz,coil.sphere_radius,coil.x_Value,coil.y_Value,coil.z_Value);
    procentPlane = [0.01,0.03,0.05,0.08,0.10,0.15];
    procentVolume = 0.05;
    DisplayFieldGradientLinearity(coil.By,wire,coil.coil_radius,coil.coil_length,procentPlane,procentVolume,coil.x_Value,coil.y_Value,coil.z_Value)
    %DisplayGradient(By,wire,x_Value,y_Value,z_Value)
    % minimumFieldAmplitude = -0.005;
    % maximumFieldAmplitude = 0.005;
    % DisplayField(Babs,By,Bz,x_Value,y_Value,z_Value,minimumFieldAmplitude,maximumFieldAmplitude)
elseif strcmpi(targetCoil,'DriveX')
    [coil.maxHomoX,coil.maxHomoY,coil.maxHomoZ,coil.meanHomoX,coil.meanHomoY,coil.meanHomoZ,coil.Bx0,coil.By0,coil.Bz0] = Homogeneity(coil.Bx,coil.By,coil.Bz,coil.sphere_radius,coil.x_Value,coil.y_Value,coil.z_Value);
    minimumFieldAmplitude = 0;
    maximumFieldAmplitude = 1;
    DisplayField(sqrt((coil.Bx./coil.By0).^2),sqrt((coil.By./coil.By0).^2),sqrt((coil.Bz./coil.By0).^2),coil.x_Value,coil.y_Value,coil.z_Value,minimumFieldAmplitude,maximumFieldAmplitude)
    procentPlane = [0.05,0.10,0.20];
    procentVolume = 0.05;
    DisplayFieldHomogeneity(coil.Bx,wire,coil.coil_radius,coil.coil_length,procentPlane,procentVolume,coil.x_Value,coil.y_Value,coil.z_Value)
elseif strcmpi(targetCoil,'DriveY')
    [coil.maxHomoX,coil.maxHomoY,coil.maxHomoZ,coil.meanHomoX,coil.meanHomoY,coil.meanHomoZ,coil.Bx0,coil.By0,coil.Bz0] = Homogeneity(coil.Bx,coil.By,coil.Bz,coil.sphere_radius,coil.x_Value,coil.y_Value,coil.z_Value);
    minimumFieldAmplitude = 0;
    maximumFieldAmplitude = 1;
    DisplayField(sqrt((coil.Bx./coil.By0).^2),sqrt((coil.By./coil.By0).^2),sqrt((coil.Bz./coil.By0).^2),coil.x_Value,coil.y_Value,coil.z_Value,minimumFieldAmplitude,maximumFieldAmplitude)
    procentPlane = [0.05,0.10,0.20];
    procentVolume = 0.05;
    DisplayFieldHomogeneity(coil.By,wire,coil.coil_radius,coil.coil_length,procentPlane,procentVolume,coil.x_Value,coil.y_Value,coil.z_Value)
elseif strcmpi(targetCoil,'DriveZ')
    [coil.maxHomoX,coil.maxHomoY,coil.maxHomoZ,coil.meanHomoX,coil.meanHomoY,coil.meanHomoZ,coil.Bx0,coil.By0,coil.Bz0] = Homogeneity(coil.Bx,coil.By,coil.Bz,coil.sphere_radius,coil.x_Value,coil.y_Value,coil.z_Value);
    minimumFieldAmplitude = 0;
    maximumFieldAmplitude = 1;
    DisplayField(sqrt((coil.Bx./coil.By0).^2),sqrt((coil.By./coil.By0).^2),sqrt((coil.Bz./coil.By0).^2),coil.x_Value,coil.y_Value,coil.z_Value,minimumFieldAmplitude,maximumFieldAmplitude)
    procentPlane = [0.05,0.10,0.20];
    procentVolume = 0.05;
    DisplayFieldHomogeneity(coil.Bz,wire,coil.coil_radius,coil.coil_length,procentPlane,procentVolume,coil.x_Value,coil.y_Value,coil.z_Value)
else
    disp('Error : unknown coil type')
end
%% generate the input file for FastHenry
WriteCoordWireFastHenry(wire,'DesignedCoil',coil.wireWidth,coil.wireThickness,25000,25000,1)

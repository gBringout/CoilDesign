% BEM methods - with shield
clear all;
close all;
addpath(genpath(fullfile('..')))
addpath(genpath(fullfile('..','..','SphericalHarmonics')))


DY_shielding_TransMagConf % the script used to make the transmag calcalation
%Q0_shielding %the script to calculate the shielding of the Quadrupole by  the solenoid

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
if shield.reduction == 1
    [shield.listNode,shield.listTriangle,shield.subBoundaries] = NodeSorting2(shield.listNode,shield.listTriangle);
    fprintf('%i border detected on the shield \n',size(shield.subBoundaries,1))
else
    shield.subBoundaries = [];
end

disp('Process the triangle')
[coil.triangle,coil.node] = processMesh(coil.listTriangle,coil.listNode);
[shield.triangle,shield.node] = processMesh(shield.listTriangle,shield.listNode);

disp('Calculate the basis function')
coil.basis = basisFunction4(coil.node, coil.triangle,coil.center);
shield.basis = basisFunction4(shield.node, shield.triangle,shield.center);

figure('Name','Data verification')
subplot(3,3,1)
trimesh(coil.listTriangle,coil.listNode(:,1),coil.listNode(:,2),coil.listNode(:,3),ones(size(coil.listNode,1),1)); 
hold all
trimesh(shield.listTriangle,shield.listNode(:,1),shield.listNode(:,2),shield.listNode(:,3),ones(size(shield.listNode,1),1)); 
axis square

%% Creating the target points
disp('Ploting the target points');
subplot(3,3,2)
plot3(rk(:,1),rk(:,2),rk(:,3),'*');
axis square

%% Calculating the Laplacian

coil.Lwp = zeros(size(coil.node,2),size(coil.node,2));
shield.Lwp = zeros(size(shield.node,2),size(shield.node,2));

if calculateLwp
    fprintf(1,'Calculating the Laplacian Operator.\n');

    coil.Lwp = Laplacian3(coil.node, coil.triangle);
    shield.Lwp = Laplacian3(shield.node, shield.triangle);
end

coupling.LwpUp = zeros(size(coil.node,2),size(shield.node,2));
coupling.LwpDown = zeros(size(shield.node,2),size(coil.node,2));
coupling.Lwp = [[coil.Lwp coupling.LwpUp];[coupling.LwpDown shield.Lwp]];
subplot(3,3,3)
imagesc(log(abs(coupling.Lwp)));
title('Lwp')
axis square
colormap(gray)
%% Attempt to calculate the Lmn matrix

coil.L = zeros(size(coil.node,2),size(coil.node,2));
shield.L = zeros(size(shield.node,2),size(shield.node,2));
coupling.LUp = zeros(size(coil.node,2),size(shield.node,2));
coupling.LDown = coupling.LUp';
coupling.L = [[coil.L coupling.LUp];[coupling.LDown shield.L]];

if calculateL
    disp('Calculating the Lmn matrix.');
    coil.L = Lmn10(coil.node, coil.triangle,coil.basis);
    shield.L = Lmn10(shield.node, shield.triangle,shield.basis);

    % Coupling [    coil    0
    %               0       shield]
    coupling.LUp = Lmn10(coil.node, coil.triangle,coil.basis,shield.node,shield.triangle,shield.basis);
    coupling.LDown = coupling.LUp';
    coupling.L = [[coil.L coupling.LUp];[coupling.LDown shield.L]];
end

subplot(3,3,4)
Llog = real(log(coupling.L));
Lilog = imag(log(coupling.L)); % I don't know why their is an imaginary part 
imagesc(Llog);
title('coupling.L (log scale)')
axis square
colormap(gray)

%% A
shield.Ax = zeros(size(shield.node,2),size(coil.node,2));
shield.Ay = zeros(size(shield.node,2),size(coil.node,2));
shield.Az = zeros(size(shield.node,2),size(coil.node,2));

if calculateA
    disp('Calculating the Ax,Ay and Az matrix.');
    [shield.Ax,shield.Ay,shield.Az] = Amn2(shield.node, shield.triangle,shield.basis,coil.node, coil.triangle,coil.basis);
end


subplot(3,3,3)
imagesc([shield.Ax,shield.Ay,shield.Az]')
title('shield.A (log scale)')
%% Attemp to calculate the Rmn matrix

coil.R = zeros(size(coil.node,2),size(coil.node,2));
shield.R = zeros(size(shield.node,2),size(shield.node,2));
coupling.RUp = zeros(size(coil.node,2),size(shield.node,2));
coupling.RDown = zeros(size(shield.node,2),size(coil.node,2));

if calculateR
    disp('Calculating the Rmn matrix.');
    coil.R = Rmn9(coil.node,coil.triangle,coil.basis,coil.wireResistivity,coil.wireThickness);
    shield.R = Rmn9(shield.node,shield.triangle,shield.basis,shield.wireResistivity,shield.wireThickness);
end

coupling.R = [[coil.R coupling.RUp];[coupling.RDown shield.R]];

subplot(3,3,5)
imagesc(coupling.R);
title('coupling.R')
axis square
colormap('gray')

%% Attemp to calculate Cn
disp('Calculating the Cn matrix.');
[coil.Cx,coil.Cy,coil.Cz] = Cn7(coil.node,coil.triangle,coil.basis,rk);
[shield.Cx,shield.Cy,shield.Cz] = Cn7(shield.node,shield.triangle,shield.basis,rk);
coil.C = [coil.Cx;coil.Cy;coil.Cz];
shield.C = [shield.Cx;shield.Cy;shield.Cz];

if strcmp(targetCoil,'dBzdx') || strcmp(targetCoil,'dBzdy') || strcmp(targetCoil,'dBzdz')
    coil.Ctarget = coil.Cz;
elseif strcmp(targetCoil,'DriveZ') || strcmp(targetCoil,'DriveY') || strcmp(targetCoil,'DriveZ')
    coil.Ctarget = [coil.Cx;coil.Cy;coil.Cz];
elseif strcmp(targetCoil, 'Quad')
    coil.Ctarget = [coil.Cx;coil.Cy;coil.Cz];
end

% Coupling [    coil    0
%               0       shield]
coupling.CUp = zeros(3*size(rk,1),size(shield.node,2));
coupling.CDown = zeros(3*size(rk,1),size(coil.node,2));
coupling.C = [[coil.C coupling.CUp];[coupling.CDown shield.C]];
subplot(3,3,6)
imagesc(coupling.C);
title('coupling.C')
caxis([min(min(coil.C)) max(max(coil.C))])
colormap('gray')

%% Reduction of the matrix to remove the border stuff

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
if shield.reduction && size(shield.subBoundaries,1)>0
    shield.nonReducedSize = size(shield.R,1);
    shield.Rfull = shield.R;
    shield.R = ReduceSquareMatrix4(shield.Rfull,shield.subBoundaries);
    shield.Lfull = shield.L;
    shield.L = ReduceSquareMatrix4(shield.Lfull,shield.subBoundaries);
    shield.Cfull = shield.C;
    shield.C = ReduceCMatrix4(shield.Cfull,shield.subBoundaries);
    shield.LwpFull = shield.Lwp;
    shield.Lwp = ReduceSquareMatrix4(shield.LwpFull,shield.subBoundaries);
else
    shield.nonReducedSize = size(shield.R,1);
    shield.Rfull = shield.R;
    shield.Lfull = shield.L;
    shield.Cfull = shield.C;
    shield.LwpFull = shield.Lwp;
end
if (coil.reduction && size(coil.subBoundaries,1)>0) || (shield.reduction && size(shield.subBoundaries,1)>0)
    coupling.LDownFull = coupling.LDown;
    coupling.LDown = ReduceLUpMatrix4(coupling.LDownFull,shield.subBoundaries,coil.subBoundaries);
    coupling.LUpFull = coupling.LUp;
    coupling.LUp = coupling.LDown';
    coupling.L = [[coil.L coupling.LUp];[coupling.LDown shield.L]];
end

subplot(3,3,7)
Llog = real(log(coupling.L));
Lilog = imag(log(coupling.L)); % I don't know why theiur is an imaginary part 
imagesc(Llog);
title('reduced coupling.L')
axis square
colormap(gray)

subplot(3,3,8)
imagesc(coupling.R);
title('reduced coupling.R')
axis square
colormap(gray)

subplot(3,3,9)
imagesc(coupling.LDownFull);
title('coupling.LDownFull')
axis square
colormap(gray)

%% Calculating the coupling factor

if calculateL && calculateR
    Mis = coupling.LDown;
    Mii = shield.L;
    Rii = shield.R;

    [Q,Lambda] = eig(Mii\Rii); % use the eigenvalue decomposition 
    % Why do we get sometime negative values? This is bad! This is because
    % some matrix are not positive semidefinite
    % we have to use the reduce version of the matrix to avoid that
    %Lambda = abs(Lambda);

    % Calculation of ai
    nbrLambda = size(Lambda,2);
    a = zeros(nbrLambda,nbrLambda);
    asin = zeros(nbrLambda,nbrLambda);
    acos = zeros(nbrLambda,nbrLambda);
    aexp = zeros(nbrLambda,nbrLambda);
    aAmplitude = zeros(nbrLambda,nbrLambda);

    t=100000/(coil.freq)+(1/4)/coil.freq; %Sin part after 100 000 oscillation
    for j=1:nbrLambda
        a(j,j) = (-2*pi*coil.freq)*(Lambda(j,j)*exp(-Lambda(j,j)*t)...
                -Lambda(j,j)*cos(2*pi*coil.freq*t)...
                -2*pi*coil.freq*sin(2*pi*coil.freq*t))/((2*pi*coil.freq)^2+Lambda(j,j)^2);
        aexp(j,j) = (-2*pi*coil.freq)*(Lambda(j,j)*exp(-Lambda(j,j)*t))/((2*pi*coil.freq)^2+Lambda(j,j)^2);
        asin(j,j) = (-2*pi*coil.freq)*(-2*pi*coil.freq)/((2*pi*coil.freq)^2+Lambda(j,j)^2);
        acos(j,j) = (-2*pi*coil.freq)*(-Lambda(j,j))/((2*pi*coil.freq)^2+Lambda(j,j)^2);
        aAmplitude(j,j) = sqrt(asin(j,j)^2+acos(j,j)^2);
    end
    coupling.Factor =           Q*a*            (-inv(Q)*inv(Mii)*Mis);
    coupling.FactorSin =        Q*asin*         (-inv(Q)*inv(Mii)*Mis);
    coupling.FactorCos =        Q*acos*         (-inv(Q)*inv(Mii)*Mis);
    coupling.FactorAmplitude =  Q*aAmplitude*   (-inv(Q)*inv(Mii)*Mis);
end

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
        coil.s_reduced = coil.x_lambda;
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
    %Opt = opti('qp',H,f,'lin',A,rl,ru,'qc',Q2,l,r,'bounds',lb,ub,'options',opts); %bounded lineare constrains
    % QPLC
    Opt = opti('qp',H,f,'lin',A,rl,ru,'bounds',lb,ub,'options',opts);
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

%% Calculation of the current distribution on the shield
if calculateL && calculateR
    if shield.reduction && coil.reduction
        shield.s_induce_reduced = coupling.Factor*coil.s_reduced;
        %s_induce_reduced = Q*a*(-inv(Q)*inv(Mii)*Mis)*s_reduced;
        shield.s_induce_reducedCos = coupling.FactorCos*coil.s_reduced;
        shield.s_induce_reducedSin = coupling.FactorSin*coil.s_reduced;
        shield.s_induce_reducedAmplitude =coupling.FactorAmplitude*coil.s_reduced;
        %s_induce_reduced = coupling.Factor*s_reduced;
        shield.s_induce = retrieveCurrentVector3(shield.s_induce_reduced,shield.subBoundaries,'noQP');
        shield.s_induceCos = retrieveCurrentVector3(shield.s_induce_reducedCos,shield.subBoundaries,'noQP');
        shield.s_induceSin = retrieveCurrentVector3(shield.s_induce_reducedSin,shield.subBoundaries,'noQP');
        shield.s_induceAmplitude = retrieveCurrentVector3(shield.s_induce_reducedAmplitude,shield.subBoundaries,'noQP');
    else
        disp('Please check that this is coherent!')
        shield.s_induce = coupling.Factor*coil.s;
        %shield.s_induce = coupling.Factor*coil.s_reduced;
        %s_induce_reduced = Q*a*(-inv(Q)*inv(Mii)*Mis)*s_reduced;
        shield.s_induceCos = coupling.FactorCos*coil.s;
        shield.s_induceSin = coupling.FactorSin*coil.s;
        shield.s_induceAmplitude =coupling.FactorAmplitude*coil.s;
        s_induce_reduced = coupling.Factor*s_reduced;
    end
    
    displayStreamFunction(shield.listTriangle,shield.s_induce,shield.listNode)
    %displayStreamFunction(shield.listTriangle,shield.s_induceAmplitude,shield.listNode)
end

%% Properties calculation

shield.B = shield.Cfull*shield.s_induce;
shield.p_dis = shield.s_induce'*shield.Rfull*shield.s_induce;
shield.e_stored = 0.5*shield.s_induce'*shield.Lfull*shield.s_induce;

coil.B = coil.Cfull*coil.s;
coil.p_dis = coil.s'*coil.Rfull*coil.s;
coil.e_stored = 0.5*coil.s'*coil.Lfull*coil.s;

%% Field loss calculation - Drive

coil.B = coil.Cfull*coil.s;
%Whe have to split the field
nbrPoint = size(rk,1);
Bx = coil.B(1:nbrPoint);
By = coil.B(nbrPoint+1:2*nbrPoint);
Bz = coil.B(2*nbrPoint+1:end);
[coil.bc,coil.bs] = getSphericalHarmonicsCoefficientMeasure7(Bx,By,Bz,degreeMax,orderMax,rk,'sch');
coil.amplitude = coil.bc(2).coefficient(1,1);
coil.phase = 0;
fprintf('Coil.\n Generated field: %2.4g T.\n Phase: %0.4g radian.\n Power Loss: %0.4g W\n',coil.amplitude,coil.phase,coil.p_dis)

shield.Bcos = shield.Cfull*shield.s_induceCos;
shield.Bsin = shield.Cfull*shield.s_induceSin;

nbrPoint = size(rk,1);
Bx = shield.Bcos(1:nbrPoint);
By = shield.Bcos(nbrPoint+1:2*nbrPoint);
Bz = shield.Bcos(2*nbrPoint+1:end);
[shield.bcCos,shield.bsCos] = getSphericalHarmonicsCoefficientMeasure7(Bx,By,Bz,degreeMax,orderMax,rk,'sch');

Bx = shield.Bsin(1:nbrPoint);
By = shield.Bsin(nbrPoint+1:2*nbrPoint);
Bz = shield.Bsin(2*nbrPoint+1:end);
[shield.bcSin,shield.bsSin] = getSphericalHarmonicsCoefficientMeasure7(Bx,By,Bz,degreeMax,orderMax,rk,'sch');

shield.phase = atan2(shield.bcCos(2).coefficient(1,1),shield.bcSin(2).coefficient(1,1)); % form A*sin(wt+phi)
shield.amplitudebis = sqrt(shield.bcSin(2).coefficient(1,1)^2+shield.bcCos(2).coefficient(1,1)^2);

shield.B = shield.Cfull*shield.s_induce;
%We have to split the field
nbrPoint = size(rk,1);
Bx = shield.B(1:nbrPoint);
By = shield.B(nbrPoint+1:2*nbrPoint);
Bz = shield.B(2*nbrPoint+1:end);
[shield.bc,shield.bs] = getSphericalHarmonicsCoefficientMeasure7(Bx,By,Bz,degreeMax,orderMax,rk,'sch');
shield.amplitude = shield.bc(2).coefficient(1,1);
fprintf('Shield.\n Generated field: %2.4g / %2.4g T.\n Phase: %0.4g radian.\n Power Loss: %0.4g W\n',shield.amplitude,shield.amplitudebis,shield.phase,shield.p_dis)
%displaySHC(shield.bc,shield.bs,2)

fprintf('Factor DC/AC: %i\n',(coil.amplitude-shield.amplitudebis)/coil.amplitude)

%% Field loss calculation - Gradient

coil.B = coil.Cfull*coil.s;
%Whe have to split the field
nbrPoint = size(rk,1);
Bx = coil.B(1:nbrPoint);
By = coil.B(nbrPoint+1:2*nbrPoint);
Bz = coil.B(2*nbrPoint+1:end);
[coil.bc,coil.bs] = getSphericalHarmonicsCoefficientMeasure7(Bx,By,Bz,degreeMax,orderMax,rk,'sch');
coil.amplitude = coil.bc(1).coefficient(2,2)/rhoReference;
coil.phase = 0;
fprintf('Coil.\n Generated field: %2.4g T/m.\n Phase: %0.4g radian.\n Power Loss: %0.4g W\n',coil.amplitude,coil.phase,coil.p_dis)

shield.Bcos = shield.Cfull*shield.s_induceCos;
shield.Bsin = shield.Cfull*shield.s_induceSin;

nbrPoint = size(rk,1);
Bx = shield.Bcos(1:nbrPoint);
By = shield.Bcos(nbrPoint+1:2*nbrPoint);
Bz = shield.Bcos(2*nbrPoint+1:end);
[shield.bcCos,shield.bsCos] = getSphericalHarmonicsCoefficientMeasure7(Bx,By,Bz,degreeMax,orderMax,rk,'sch');

Bx = shield.Bsin(1:nbrPoint);
By = shield.Bsin(nbrPoint+1:2*nbrPoint);
Bz = shield.Bsin(2*nbrPoint+1:end);
[shield.bcSin,shield.bsSin] = getSphericalHarmonicsCoefficientMeasure7(Bx,By,Bz,degreeMax,orderMax,rk,'sch');

shield.phase = atan2(shield.bcCos(2).coefficient(1,1),shield.bcSin(2).coefficient(1,1)); % form A*sin(wt+phi)
shield.amplitudebis = sqrt(shield.bcSin(2).coefficient(1,1)^2+shield.bcCos(2).coefficient(1,1)^2);

shield.B = shield.Cfull*shield.s_induce;
%We have to split the field
nbrPoint = size(rk,1);
Bx = shield.B(1:nbrPoint);
By = shield.B(nbrPoint+1:2*nbrPoint);
Bz = shield.B(2*nbrPoint+1:end);
[shield.bc,shield.bs] = getSphericalHarmonicsCoefficientMeasure7(Bx,By,Bz,degreeMax,orderMax,rk,'sch');
shield.amplitude = shield.bc(1).coefficient(2,2)/rhoReference;
fprintf('Shield.\n Generated field: %2.4g / %2.4g T/m.\n Phase: %0.4g radian.\n Power Loss: %0.4g W\n',shield.amplitude,shield.amplitudebis,shield.phase,shield.p_dis)
%displaySHC(shield.bc,shield.bs,2)

fprintf('Factor DC/AC: %i\n',(coil.amplitude-shield.amplitudebis)/coil.amplitude)


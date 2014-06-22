% BEM methods - with shield
clear all;
close all;
addpath('.\examples');
addpath('.\data');
addpath('.\coreFunctions');
addpath('.\displayFunctions');
addpath('.\otherFunctions');
addpath(genpath('..\..\SphericalHarmonics\'));

%fieldByMPI % Drive coil FFL
%fieldByPig
%fieldCTMPI
%Receive_Mandy_fieldByMPI
%gradientBzyMRI
%gradientByyMPI % used to make the quadrupole
%fieldByOPENMPI % Used to make the flat coils
%gradientByyHercules
%FFLstatic
%ShimmingBy22MPI%Shimming for FFL
%ShieldMPI
%DY_p_ShieldMPI% shield coil with a cylinder
DY_shielding_TransMagConf % the script used to make the transmag calcalation
%DriveY_Human_Higor % induced cirrent density in the body
%GY_MRI_Human % Validation with the experiment of Sanchez
%DY_shielding_TransMagConf % the script used to make the transmag calcalation
%InducedVoltage_Test1
%InducedVoltage_Test2

if strcmp(optimizationType,'standardTikhonov') || strcmp(optimizationType,'generalizedTikhonov')
    % Please download the regularization tools of: http://www.imm.dtu.dk/~pcha/Regutools/
    % and add the folder to matlab' path:
    addpath('C:\Users\bringout.IMT\Desktop\Dropbox\Ph.D\Software\regu');
elseif strcmp(optimizationType,'QP')
    % Please download the OPTI TOOLBOX from
    % http://www.i2c2.aut.ac.nz/Wiki/OPTI/, install it
    % and add the folder to matlab' path :
    addpath('C:\Users\bringout.IMT\Desktop\Dropbox\Ph.D\Software\OptiToolbox')
end

tStart=tic;
%% Display the Mesh

%shield.listNode(:,1) = shield.listNode(:,1)*(0.08/0.3);
%shield.listNode(:,2) = shield.listNode(:,2)*(0.08/0.3);
%shield.listNode(:,3) = shield.listNode(:,3)*(0.08/0.3);
%shield.listNode(:,3) = shield.listNode(:,3)-0.08;
%shield.listNode(:,1) = shield.listNode(:,1)-0.13;
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

% Calculate the basis function of the mesh.
%node = coil.node;
%triangle=coil.triangle;

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
%[coil.Cx,coil.Cy,coil.Cz] = Cn5(coil.node,coil.basis,rk);
[coil.Cx,coil.Cy,coil.Cz] = Cn7(coil.node,coil.triangle,coil.basis,rk);
[shield.Cx,shield.Cy,shield.Cz] = Cn7(shield.node,shield.triangle,shield.basis,rk);
%[shield.Cx,shield.Cy,shield.Cz] = Cn5(shield.node,shield.basis,rk);
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
    %Mis = coupling.LDown11;
    %Mii = shield.L11;
%     Mis = coupling.LDownFull;
%     Mii = shield.Lfull;
%     Rii = shield.Rfull;
    Mis = coupling.LDown; % test with not the full matrix
    Mii = shield.L;
    Rii = shield.R;

    %mean(mean(inv(Mii)*Mii - eye(440))); % this should be very small (10^-18). ok this is working
    %mean(mean(Mii*inv(Mii) - eye(440))); % this should be very small (10^-18). ok this is working
    %mean(mean(inv(Rii)*Rii - eye(440))); % this should be very small (10^-18). ok this is not working with the non-Reduced matrix !
    [Q,Lambda] = eig(Mii\Rii); % use the eigenvalue decomposition % Why do we get sometime negative values? This is bad!
    %[Q,Lambda] = eig(Rii,Mii); % use the eigenvalue decomposition
    Lambda = abs(Lambda);
    %X = Q*Lambda*(inv(Q)); % A' is the linear algebraic transpose of A
    %Mii\Rii - X;
    %mean(mean(Mii\Rii - X)); % should be zero. ok
    %Q/(Q); % This should be the identity matrix. A' is the linear algebraic transpose of A

    % Calculation of ai
    nbrLambda = size(Lambda,2);
    a = zeros(nbrLambda,nbrLambda);
    asin = zeros(nbrLambda,nbrLambda);
    acos = zeros(nbrLambda,nbrLambda);
    aexp = zeros(nbrLambda,nbrLambda);
    aAmplitude = zeros(nbrLambda,nbrLambda);
    % for a big t
    % the first part (exp(-lamdai...) is zero
    %t=100000/(freq)+0; %Cos part after 100.000 oscillation
    t=100000/(coil.freq)+(1/4)/coil.freq; %Sin part after 100 000 oscillation
    % Before was: 
    %t=100000/(2*pi*freq)+(pi/2)/(2*pi*freq); %Sin part after 100 000 oscillation
    for j=1:size(Lambda,2)
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
    %reg = 10^-6;
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

    H = coil.R;
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

%i=x;
displayStreamFunction(coil.listTriangle,coil.s,coil.listNode)
%[jTri,jAbsTri,j,jAbs] = calculateCurrentDensity(coil.s,coil.node,coil.triangle,coil.basis);
%displayCurrentDensity(coil.listTriangle,coil.listNode,jTri,jAbsTri)

% displayStreamFunction(coil.listTriangle,j,coil.listNode)
%caxis([-1000 1000])
%figure;semilogy(results(2,1:20))
%0.5*coil.s'*coil.Lfull*coil.s;
% i'*-LwpFull*i
% i'*Rfull*i
% s_reduced'*coil.R*s_reduced
% min(eig(R))

% %% Exctract the wire based on a simple method
% %wireThickness = 0.01;
% %nbrContour = startingWireNumber;
% startingWireNumber = 9;
% disp('Wire extraction')
% clear('wire');
% [wire,nbrwire,minimalWireSpacing] = NbrWireOptimisation3(coil.node,coil.triangle,coil.s,coil.wireWidth,startingWireNumber,coil.distanceBetween2Wire,coil.rateIncreasingWire);
% wire = cleanWire(wire);
% %minimalWireSpacing = DistanceBetweenWire(wire);
% %radius(wire,'xy');
% displayWire(wire);
% Length2(wire);
% %displayWireAsPoint(wire);
% %displayWireSubFigure(wire);
% %[wire2] = WireRotation(wire,pi/4,'z');
% %displayWire(wire2);
% % wire(19) = [];

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
    displayStreamFunction(shield.listTriangle,shield.s_induceAmplitude,shield.listNode)
    %[jTri,jAbsTri,j,jAbs] = calculateCurrentDensity(shield.s_induce,shield.node,shield.triangle,shield.basis);
    %displayCurrentDensity(shield.listTriangle,shield.listNode,jTri,jAbsTri)
    %[jTri,jAbsTri,j,jAbs] = calculateCurrentDensity(coil.s,coil.node,coil.triangle,coil.basis);
    %displayCurrentDensity(coil.listTriangle,coil.listNode,jTri,jAbsTri)
    %displayStreamFunction(shield.listTriangle,j,shield.listNode)
    toc(tStart)

    % plot(node(row,1),
    %hold all
    %plot3(node(row,1),node(row,2),node(row,3),'r*')
end

%% Test voltage 1 (do not work)
% if shield.reduction
%     res = shield.Rfull*shield.s_induce; % resistance associated to each node
%     %[~,jabs] = calculateCurrentDensity(shield.triangle,shield.s_induce,shield.node,shield.nodeLinkToTriangle,shield.basis);
%     %dcVol = res.*jabs; %voltage associated to each node
%     dcVol = res.*shield.s_induce; %voltage associated to each node
%     indu = shield.Lfull*shield.s_induce; %inductance of the element
%     acVolSelf = indu*2*pi*freq.*shield.s_induce;
%     acVolIndu = (coupling.LDownFull*coil.s)*2*pi*freq.*shield.s_induce;
% else 
%     res = shield.R*shield.s_induceAmplitude; % resistance associated to each node
%     %[~,jabs] = calculateCurrentDensity(shield.triangle,shield.s_induce,shield.node,shield.nodeLinkToTriangle,shield.basis);
%     %dcVol = res.*jabs; %voltage associated to each node
%     dcVol = res.*shield.s_induceAmplitude; %voltage associated to each node
%     indu = shield.L*shield.s_induceAmplitude; %inductance of the element
%     acVolSelf = indu*2*pi*freq.*shield.s_induceAmplitude;
%     acVolIndu = (coupling.LDown*coil.s)*2*pi*freq.*shield.s_induceAmplitude;
% end
% vol = dcVol+acVolSelf+acVolIndu;
% displayVoltage(shield.listTriangle,vol,shield.listNode)
%% Test Voltage 2 (do not seems to work either)
% 
% res = shield.Rfull*shield.s_induce; % resistance associated to each node
% %[~,~,j,jabs] = calculateCurrentDensity(shield.s_induce,shield.node,shield.triangle,shield.basis);
% % don't we have to just take the self resistance and inductance into
% % account?
% dcVol(:,1) = shield.Rfull*(j(:,1)); %DC voltage x associated to each node
% dcVol(:,2) = shield.Rfull*(j(:,2)); %DC voltage y associated to each node
% dcVol(:,3) = shield.Rfull*(j(:,3)); %DC voltage z associated to each node
% 
% acVol(:,1) = shield.Lfull*j(:,1)*2*pi*freq; %AC voltage x  of the element
% acVol(:,2) = shield.Lfull*j(:,2)*2*pi*freq; %AC voltage x  of the element
% acVol(:,3) = shield.Lfull*j(:,3)*2*pi*freq; %AC voltage x  of the element
% 
% 
% %vol = sqrt((acVol(:,1)).^2+(acVol(:,2)).^2+(acVol(:,3)).^2);
% vol = sqrt((dcVol(:,1)).^2+(dcVol(:,2)).^2+(dcVol(:,3)).^2);
% 
% %displayCurrentDensity(shield.listTriangle,shield.listNode,jTri,jAbsTri)
% displayVoltage(shield.listTriangle,vol,shield.listNode); title('dc voltage')
% 
% vol = sqrt((dcVol(:,1)+acVol(:,1)).^2+(dcVol(:,2)+acVol(:,2)).^2+(dcVol(:,3)+acVol(:,3)).^2);
% displayVoltage(shield.listTriangle,vol,shield.listNode); title('dc + ac voltage')

%% Test voltage 3
% scale the results to get the correct gradient :

coil.B = coil.Cfull*coil.s;
%Whe have to split the field
nbrPoint = size(rk,1);
Bx = coil.B(1:nbrPoint);
By = coil.B(nbrPoint+1:2*nbrPoint);
Bz = coil.B(2*nbrPoint+1:end);
[coil.bc,coil.bs] = getSphericalHarmonicsCoefficientMeasure7(Bx,By,Bz,degreeMax,orderMax,rk,'sch');


if strcmp(targetCoil,'dBzdx') 
    coil.amplitude = coil.bc(3).coefficient(2,2); % X gradient
elseif strcmp(targetCoil,'dBzdy')
    coil.amplitude = coil.bs(3).coefficient(2,2); % Y gradient
elseif strcmp(targetCoil,'dBzdz')
    coil.amplitude = coil.bc(3).coefficient(2,1); % Z gradient
end

coil.correction = 0.03*rhoReference/coil.amplitude;

Ex = 2*pi*freq*shield.Ax*coil.s*coil.correction;
Ey = 2*pi*freq*shield.Ay*coil.s*coil.correction;
Ez = 2*pi*freq*shield.Az*coil.s*coil.correction;
E = sqrt(Ex.^2+Ey.^2+Ez.^2);
displayVoltage(shield.listTriangle,E,shield.listNode); title('ac voltage')
toc(tStart)

%% Test 4

[~,~,j,jabs] = calculateCurrentDensity(shield.s_induce,shield.node,shield.triangle,shield.basis);
displayCurrentDensity(shield.listTriangle,shield.listNode,jTri,jAbsTri)

if false
    %% Field calculation

    shield.B = shield.Cfull*shield.s_induce;
    shield.p_dis = shield.s_induce'*shield.Rfull*shield.s_induce;
    shield.e_stored = 0.5*shield.s_induce'*shield.Lfull*shield.s_induce;

    coil.B = coil.Cfull*coil.s;
    coil.p_dis = coil.s'*coil.Rfull*coil.s;
    coil.e_stored = 0.5*coil.s'*coil.Lfull*coil.s;

    %% Try the spherical harmonics decomposition

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
    %displaySHC(coil.bc,coil.bs,2)


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
    %Whe have to split the field
    nbrPoint = size(rk,1);
    Bx = shield.B(1:nbrPoint);
    By = shield.B(nbrPoint+1:2*nbrPoint);
    Bz = shield.B(2*nbrPoint+1:end);
    [shield.bc,shield.bs] = getSphericalHarmonicsCoefficientMeasure7(Bx,By,Bz,degreeMax,orderMax,rk,'sch');
    shield.amplitude = shield.bc(2).coefficient(1,1);
    fprintf('Shield.\n Generated field: %2.4g / %2.4g T.\n Phase: %0.4g radian.\n Power Loss: %0.4g W\n',shield.amplitude,shield.amplitudebis,shield.phase,shield.p_dis)
    %displaySHC(shield.bc,shield.bs,2)

    fprintf('Factor DC/AC: %i\n',(coil.amplitude+shield.amplitudebis)/coil.amplitude)

    %% Compare the shiel diameter
    shieldDiameter = max(shield.node(:,2));
    filename = sprintf('Coil_0.119_Shield_%4.5g.mat',shieldDiameter);
    save(filename,'-struct','shield','bc')
    save(filename,'-struct','shield','bs','-append')
    %save('Coil_0.119.mat','-struct','coil','bc')
    %save('Coil_0.119.mat','-struct','coil','bs','-append')
    disp('Saved');

    %% Analysing the results
    % clear all
    % close all
    %displaySHC(coil.bc,coil.bs,2)
    addpath('C:\Users\bringout.IMT\Desktop\Dropbox\Ph.D\Software\matlab2tikz_0.4.1\src')

    coil = load('Coil_0.119.mat');
    shieldDiameter = [0.135 0.14 0.145 0.150 0.16 0.17 0.18];
    %coefToKeep = [[1 1] ;[3 1];[3 3];[5 1]];% Shield 160
    %coefToKeep = [[1 1] ;[2 2];[3 1];[3 3];[4 4];[5 1];[5 5]];% First version of the article
    coefToKeep = [[1 1] ;[3 1];[3 3];[5 1];[5 3];[7 1]];% Coil

    % Initialize with the coil data
    i=1;
    for j=1:size(coefToKeep,1)
        if coefToKeep(j,2)>0
            SHC(i,j+1) = coil.bc(2).coefficient(coefToKeep(j,1),coefToKeep(j,2));
        else
            SHC(i,j+1) = coil.bs(2).coefficient(coefToKeep(j,1),abs(coefToKeep(j,2)));
        end
    end


    for i=1:size(shieldDiameter,2)
        filename = sprintf('Coil_0.119_Shield_%4.5g.mat',shieldDiameter(i));
        load(filename)
        for j=1:size(coefToKeep,1)
            SHC(i+1,1) = shieldDiameter(i);
            if coefToKeep(j,2)>0
                SHC(i+1,j+1) = bc(2).coefficient(coefToKeep(j,1),coefToKeep(j,2));
            else
                SHC(i+1,j+1) = bs(2).coefficient(coefToKeep(j,1),abs(coefToKeep(j,2)));
            end
            SHC2(i+1,j+1) = (SHC(1,j+1)+SHC(i+1,j+1))/SHC(1,j+1);
        end
        %displaySHC(bc,bs,2)
    end

    figure
    strLegend = '';
    for i=1:size(coefToKeep,1)
        %subplot(2,3,i)
        %plot(SHC(:,1),SHC(:,i+1)/max(SHC(:,i+1))) The normalized
        %progression of harmonics
        strLegend = [strLegend;sprintf('%gth order-%gth degree',coefToKeep(i,1)-1,coefToKeep(i,2)-1)];
        plot(SHC(2:end,1),(SHC(1,i+1)+SHC(2:end,i+1))/SHC(1,i+1))
        %plot(SHC(2:end,1),(-SHC(2:end,i+1))/SHC(1,i+1))
        %plot(SHC(2:end,1),SHC(2:end,i+1))
        hold all
    end
    legend(strLegend)

    %matlab2tikz

    %% Calculation of the field discritized coil
    % disp('Field calculation with the discritized model')
    % 
    % [Bx,By,Bz] = Field2(wire,current,x_Value,y_Value,z_Value);
    % Babs = sqrt(Bx.^2+By.^2+Bz.^2);
    % %DisplayField(Bx,By,Bz,x_Value,y_Value,z_Value,-1*10^-3,1*10^-3);
    % if strcmpi(targetFieldType,'gradient')
    %     [gradientLinearity,gradientValueCenter] = Linearity(Bx,By,Bz,sphere_radius,x_Value,y_Value,z_Value);
    %     coil_radius = Radius(wire,'xy');
    %     coil_length = CoilLength(wire);
    %     procentPlane = [0.01,0.03,0.05,0.08,0.10,0.15];
    %     procentVolume = 0.05;
    %     DisplayFieldGradientLinearity(By,wire,coil_radius,coil_length,procentPlane,procentVolume,x_Value,y_Value,z_Value)
    %     %DisplayGradient(By,wire,x_Value,y_Value,z_Value)
    %     % minimumFieldAmplitude = -0.005;
    %     % maximumFieldAmplitude = 0.005;
    %     % DisplayField(Babs,By,Bz,x_Value,y_Value,z_Value,minimumFieldAmplitude,maximumFieldAmplitude)
    % elseif strcmpi(targetFieldType,'field')
    %     [maxHomoX,maxHomoY,maxHomoZ,meanHomoX,meanHomoY,meanHomoZ,Bx0,By0,Bz0] = Homogeneity(Bx,By,Bz,sphere_radius,x_Value,y_Value,z_Value);
    %     minimumFieldAmplitude = 0;
    %     maximumFieldAmplitude = 1;
    %     DisplayField(sqrt((Bx./By0).^2),sqrt((By./By0).^2),sqrt((Bz./By0).^2),x_Value,y_Value,z_Value,minimumFieldAmplitude,maximumFieldAmplitude)
    %     
    %     coil_radius = 0.1;
    %     coil_length = 0.2;
    %     procentPlane = [0.01,0.03,0.05,0.08,0.10,0.15];
    %     procentVolume = 0.05;
    %     DisplayFieldHomogeneity(By,wire,coil_radius,coil_length,procentPlane,procentVolume,x_Value,y_Value,z_Value)
    % 
    % else
    %     disp('Error : unknown targetFieldType')
    % end
    % 
    % rCoils = Resistance(Length2(wire),wireSurface,wireConductivity); %5mmx5mm is equivalent to the number used in Timo's coil or the 7.5*7.5 litz wire
    %
    %% Calculation of Value for the Ideal coil
    % disp('Calculation of Value according to the IBEM method')
    % %lTotale = 0;
    % %rTotale = 0;
    % if reduction == 1
    %     [~,bxIdeal] = IdealField(i,Cx,position);
    %     [~,byIdeal] = IdealField(i,Cy,position);
    %     [~,bzIdeal] = IdealField(i,Cz,position);
    %     rTotale = IdealResistance(x_lambda,R);
    %     lTotale = IdealInductance(x_lambda,L);
    % else
    %     BxIdeal = IdealField(I./max(I),Cx,position);
    %     ByIdeal = IdealField(I./max(I),Cy,position);
    %     BzIdeal = IdealField(I./max(I),Cz,position);
    %     rTotale = IdealResistance(x_lambda,R);
    %     lTotale = IdealInductance(x_lambda,L);
    % end
    % 
    % x0 = (size(x_Value,2)-1)/2 +1;
    % y0 = (size(y_Value,2)-1)/2 +1;
    % z0 = (size(z_Value,2)-1)/2 +1;
    % idealBy0 = B(x0,y0,z0);
    % if strcmpi(targetFieldType,'gradient')
    %     disp('TODO');
    % elseif strcmpi(targetFieldType,'field')
    %     maxIdealHomogeneity = max(max(max(B./idealBy0)))-1;
    % else
    %     disp('Error : unknown targetFieldType')
    % end
    % 
    % %% Efficiency opti
    % minRegu = -12;
    % maxRegu = -6;
    % nbrPoint = 10;
    % reguParameter = logspace(minRegu,maxRegu,nbrPoint);
    % Results = zeros(nbrPoint,5);
    % 
    % for k = 1:nbrPoint
    %     tic;
    %     reg = reguParameter(k);
    %     [x_lambda,rho,eta] = tikhonov(U2,s,V,btarget',reg); 
    %     if reduction
    %         i = retrieveCurrentVector(x_lambda,nonReducedSize,nbrBorderNode);
    %     else
    %         I = x_lambda;
    %         i = I;
    %         i(1:nbrBorderNode) = 0; 
    %     end
    %     
    %     clear('wire');
    %     [wire,nbrwire,minimalWireSpacing] = NbrWireOptimisation3(node,triangle,-i,wireThickness,startingWireNumber,distanceBetween2Wire,rateIncreasingWire);
    %     wire = cleanWire(wire);
    %     
    %     [Bx,By,Bz] = Field2(wire,current,x_Value,y_Value,z_Value);
    %     Babs = sqrt(Bx.^2+By.^2+Bz.^2);
    %     [maxHomoX,maxHomoY,maxHomoZ,meanHomoX,meanHomoY,meanHomoZ,Bx0,By0,Bz0] = Homogeneity(Bx,By,Bz,sphere_radius,x_Value,y_Value,z_Value);
    %     Results(k,1) = reg;
    %     Results(k,2) = By0;
    %     Results(k,3) = current;
    %     Results(k,4) = By0/current;
    %     Results(k,5) = size(wire,2);
    %     fprintf('Loop %i on %i - Done in %5.0f sec.\n',k,nbrPoint,toc);
    % end
    % 
    %% Results scanning
    % minRegu = 1*10^-10;
    % maxRegu = 5*10^-3;
    % nbrPoint = 10;
    % current = 1;
    % x_Value = -0.005:0.001:0.005;
    % y_Value = -0.005:0.001:0.005;
    % z_Value = -0.005:0.001:0.005;
    % sphere_radius = 0.0250;
    % reguParameter = linspace(minRegu,maxRegu,nbrPoint);
    % Results = zeros(nbrPoint,5);
    % figure
    % for j=1:size(reguParameter,2)
    %     reg = reguParameter(j);
    %     % Standard Tikho (Identity)
    %     %[x_lambda,rho,eta] = tikhonov(U2,s,V2,btarget',reg);
    %     % Genralized Tikho
    %     [x_lambda,rho,eta] = tikhonov(U,sm,X,btarget',reg);
    %     
    %     if reduction
    %         i = retrieveCurrentVector(x_lambda,nonReducedSize,nbrBorderNode);
    %     else
    %         I = x_lambda;
    %         i = I;
    %         i(1:nbrBorderNode) = 0; 
    %     end
    %     subplot(2,5,j)
    %     displayCurrentDensity(triangle,i,node)
    %     
    % %    disp('Wire extraction')
    % %    clear('wire');
    % %    [wire,nbrwire,minimalWireSpacing] = NbrWireOptimisation3(node,triangle,-i,wireThickness,startingWireNumber,distanceBetween2Wire,rateIncreasingWire);
    % %    wire = cleanWire(wire);
    % %    
    % % 
    % %     disp('Field calculation with the discritized model')
    % %     [Bx,By,Bz] = Field2(wire,current,x_Value,y_Value,z_Value);
    % %     %DisplayField(Bx,By,Bz,x_Value,y_Value,z_Value,-1*10^-3,1*10^-3);
    % %     [gradientLinearity,gradientValueCenter] = Linearity(Bx,By,Bz,sphere_radius,x_Value,y_Value,z_Value);
    % %     
    % %     Results(j,1) = reg;
    % %     Results(j,2) = gradientValueCenter(1,1);
    % %     Results(j,3) = Length2(wire);
    % %     Results(j,4) = nbrwire;
    % %     Results(j,5) = minimalWireSpacing.value;
    % end

    %% Look for the time signal of the stuff
    % 
    % tstart=100000/(freq)+0; %Cos part after 100.000 oscillation
    % time = tstart+linspace(0,1/freq,100);
    % Pprime = (-inv(Q)*inv(Mii)*Mis)*coil.s_reduced;
    % amplitude = zeros(10,100);
    % 
    % for index=1:size(time,2)
    %     %test(index) = 
    %     %t=100000/(2*pi*freq)+(pi/2)/(2*pi*freq); %Sin part after 100 000 oscillation
    %     t=time(index);
    %     for j=1:size(Lambda,2)
    %         a(j,j) = (-2*pi*freq)*(Lambda(j,j)*exp(-Lambda(j,j)*t)...
    %                 -Lambda(j,j)*cos(2*pi*freq*t)...
    %                 -2*pi*freq*sin(2*pi*freq*t))/((2*pi*freq)^2+Lambda(j,j)^2);
    %         
    %         asin(j,j) = (-2*pi*freq)*(-2*pi*freq)*sin(2*pi*freq*t)/((2*pi*freq)^2+Lambda(j,j)^2);
    %     end
    % 
    %     s_induce_reduced = Q*a*Pprime;
    %     s_induce_reducedSin = Q*asin*Pprime;
    %     %s_induce_reduced = coupling.Factor*s_reduced;
    %     s_induce = retrieveCurrentVector(s_induce_reduced,shield.nonReducedSize,shield.nbrBorderNode);
    %     s_induceSin = retrieveCurrentVector(s_induce_reducedSin,shield.nonReducedSize,shield.nbrBorderNode);
    %     shield.s = s_induce;
    %     % Field calculation
    % 
    %     shield.B = shield.CwholeField*shield.s;
    %     BSin = shield.CwholeField*s_induceSin;
    %     %Whe have to split the field
    %     nbrPoint = size(rk,1);
    %     Bx = shield.B(1:nbrPoint);
    %     By = shield.B(nbrPoint+1:2*nbrPoint);
    %     Bz = shield.B(2*nbrPoint+1:end);
    %     [shield.bc,shield.bs] = getSphericalHarmonicsCoefficientMeasure2(Bx,By,Bz,sqrt(nbrPoint),rk);
    %     %fprintf('Shield. Generated field: %2.2g T. Power Loss: %0.4g W\n',shield.bc(2).coefficient(1,1),shield.p_dis)
    %     %displaySHC(shield.bc,shield.bs,2)
    % 
    %     amplitude(1,index) = shield.bc(2).coefficient(1,1); % Time dependan amplitude gene by the shield
    %     amplitude(2,index) = cos(2*pi*freq*time(index)); %Cosinus function
    %     amplitude(3,index) = sin(2*pi*freq*time(index)); %Sinus function
    %     amplitude(4,index) = a(15,15); %Exemple a factor a
    %     
    %     Bx = BSin(1:nbrPoint);
    %     By = BSin(nbrPoint+1:2*nbrPoint);
    %     Bz = BSin(2*nbrPoint+1:end);
    %     [shield.bc,shield.bs] = getSphericalHarmonicsCoefficientMeasure2(Bx,By,Bz,sqrt(nbrPoint),rk);
    %     
    %     amplitude(7,index) = shield.bc(2).coefficient(1,1); % Field just with the sinus function
    % end
    % 
    % for index=1:size(time,2)
    %     amplitude(5,index) = coil.bc(2).coefficient(1,1)*sin(2*pi*freq*time(index)); %Field generated by the coil
    %     amplitude(6,index) = amplitude(5,index) + amplitude(1,index);%Both field together (shield + coil), corrected with the phase
    % end
    % fitting= amplitude(1,:);
    % 
    % 
    % figure
    % hold all;
    % plot(amplitude(1,:))
    % plot(amplitude(5,:))
    % plot(amplitude(6,:))
    % plot(0.02*amplitude(2,:))
    % plot(0.02*amplitude(3,:))
    % plot(amplitude(7,:))
    % legend('Shield','Coil','Both','Cosine','Sine','sine part of BEM')
    % grid on
    % 
    % max(amplitude(5,:))
    % max(amplitude(1,:))
    % max(amplitude(6,:))
    % amplitude(6,26)

    %% Test to get the sinus cosinus stuff
    % 
    % shield.bcCos(2).coefficient(1,1)
    % 
    % shield.bcSin(2).coefficient(1,1)

    %% Test solution to calculate the amplitude of the shield signal
    % for j=1:size(Lambda,2)
    %     %a(j,j) = (-2*pi*freq)*(Lambda(j,j)*exp(-Lambda(j,j)*t)...
    %     %        -Lambda(j,j)*cos(2*pi*freq*t)...
    %     %        -2*pi*freq*sin(2*pi*freq*t))/((2*pi*freq)^2+Lambda(j,j)^2);
    %     %asin(j,j) = (-2*pi*freq)*(-2*pi*freq)/((2*pi*freq)^2+Lambda(j,j)^2);
    %     %acos(j,j) = (-2*pi*freq)*(-Lambda(j,j))/((2*pi*freq)^2+Lambda(j,j)^2);
    %     aAmplitude(j,j) = sqrt(asin(j,j)^2+acos(j,j)^2);
    %     aPhase(j,j) = atan(acos(j,j)/asin(j,j));
    % end
    % 
    % s_induce_reducedCos = Q*acos*(-inv(Q)*inv(Mii)*Mis)*s_reduced;
    % s_induce_reducedSin = Q*asin*(-inv(Q)*inv(Mii)*Mis)*s_reduced;
    % s_induce_reducedAmplitude = Q*aAmplitude*(-inv(Q)*inv(Mii)*Mis)*s_reduced;
    % s_induce_reducedPhase = Q*aPhase*(-inv(Q)*inv(Mii)*Mis)*s_reduced;
    % 
    % 
    % s_induce = retrieveCurrentVector(s_induce_reduced,shield.nonReducedSize,shield.nbrBorderNode);
    % s_induceCos = retrieveCurrentVector(s_induce_reducedCos,shield.nonReducedSize,shield.nbrBorderNode);
    % s_induceSin = retrieveCurrentVector(s_induce_reducedSin,shield.nonReducedSize,shield.nbrBorderNode);
    % s_induceAmplitude = retrieveCurrentVector(s_induce_reducedAmplitude,shield.nonReducedSize,shield.nbrBorderNode);
    % s_inducePhase = retrieveCurrentVector(s_induce_reducedPhase,shield.nonReducedSize,shield.nbrBorderNode);
    % 
    % 
    % shield.Bcos = shield.CwholeField*s_induceCos;
    % shield.Bsin = shield.CwholeField*s_induceSin;
    % shield.BAmplitude = shield.CwholeField*s_induceAmplitude;
    % shield.BPhase = shield.CwholeField*s_inducePhase;
    % 
    % nbrPoint = size(rk,1);
    % 
    % Bx = shield.Bcos(1:nbrPoint);
    % By = shield.Bcos(nbrPoint+1:2*nbrPoint);
    % Bz = shield.Bcos(2*nbrPoint+1:end);
    % [shield.bcCos,shield.bsCos] = getSphericalHarmonicsCoefficientMeasure2(Bx,By,Bz,sqrt(nbrPoint),rk);
    % 
    % Bx = shield.Bsin(1:nbrPoint);
    % By = shield.Bsin(nbrPoint+1:2*nbrPoint);
    % Bz = shield.Bsin(2*nbrPoint+1:end);
    % [shield.bcSin,shield.bsSin] = getSphericalHarmonicsCoefficientMeasure2(Bx,By,Bz,sqrt(nbrPoint),rk);
    % 
    % Bx = shield.BAmplitude(1:nbrPoint);
    % By = shield.BAmplitude(nbrPoint+1:2*nbrPoint);
    % Bz = shield.BAmplitude(2*nbrPoint+1:end);
    % [shield.bcAmplitude,shield.bsAmplitude] = getSphericalHarmonicsCoefficientMeasure2(Bx,By,Bz,sqrt(nbrPoint),rk);
    % shield.bcAmplitude(2).coefficient(1,1)
    % 
    % Bx = shield.BPhase(1:nbrPoint);
    % By = shield.BPhase(nbrPoint+1:2*nbrPoint);
    % Bz = shield.BPhase(2*nbrPoint+1:end);
    % [shield.bcPhase,shield.bsPhase] = getSphericalHarmonicsCoefficientMeasure2(Bx,By,Bz,sqrt(nbrPoint),rk);
    % shield.bcPhase(2).coefficient(1,1) % This dont works!
    % 
    % shield.amplitude = sqrt((shield.bcCos(2).coefficient(1,1))^2+(shield.bcSin(2).coefficient(1,1))^2);
    % shield.phase = atan(shield.bcCos(2).coefficient(1,1)/shield.bcSin(2).coefficient(1,1));

    %% Efficiency of the coil function of the frequency
    % frequencyArray = [100:500:1000 2000:5000:20000 30000:50000:100000];
    % tic
    % for ii = 1:size(frequencyArray,2)
    %     freq = frequencyArray(ii);
    %     nbrLambda = size(Lambda,2);
    % %    a = zeros(nbrLambda,nbrLambda);
    %     asin = zeros(nbrLambda,nbrLambda);
    %     acos = zeros(nbrLambda,nbrLambda);
    % %    aexp = zeros(nbrLambda,nbrLambda);
    % %    aAmplitude = zeros(nbrLambda,nbrLambda);
    % % we have to correct the resistance for the skin effect
    %     %shield.skinDepth = sqrt(2*shield.rhoCopper/(shield.muCopper*2*pi*freq));
    %     %shield.wireThickness = 4*shield.skinDepth; % (meter) Thickness of the shield. Should be 4*skin depth
    %     newSkinDeph =  sqrt(2*shield.rhoCopper/(shield.muCopper*2*pi*freq));
    %     newThickness = 4*newSkinDeph;
    %     ratioChange(ii) = 1;%newThickness/shield.wireThickness;
    %     [Q,Lambda] = eig(Mii\(Rii*ratioChange(ii))); % use the eigenvalue decomposition
    %     
    %     % for a big t
    %     % the first part (exp(-lamdai...) is zero
    %     %t=100000/(freq)+0; %Cos part after 100.000 oscillation
    %     t=100000/(freq)+(1/4)/freq; %Sin part after 100 000 oscillation
    %     % Before was: 
    %     %t=100000/(2*pi*freq)+(pi/2)/(2*pi*freq); %Sin part after 100 000 oscillation
    %     for j=1:size(Lambda,2)
    %         a(j,j) = (-2*pi*freq)*(Lambda(j,j)*exp(-Lambda(j,j)*t)...
    %                 -Lambda(j,j)*cos(2*pi*freq*t)...
    %                 -2*pi*freq*sin(2*pi*freq*t))/((2*pi*freq)^2+Lambda(j,j)^2);
    % %        aexp(j,j) = (-2*pi*freq)*(Lambda(j,j)*exp(-Lambda(j,j)*t))/((2*pi*freq)^2+Lambda(j,j)^2);
    %         asin(j,j) = (-2*pi*freq)*(-2*pi*freq)/((2*pi*freq)^2+Lambda(j,j)^2);
    %         acos(j,j) = (-2*pi*freq)*(-Lambda(j,j))/((2*pi*freq)^2+Lambda(j,j)^2);
    % %        aAmplitude(j,j) = sqrt(asin(j,j)^2+acos(j,j)^2);
    %     end
    %     coupling.Factor = Q*a*(-inv(Q)*inv(Mii)*Mis);
    %     coupling.FactorSin = Q*asin*(-inv(Q)*inv(Mii)*Mis);
    %     coupling.FactorCos = Q*acos*(-inv(Q)*inv(Mii)*Mis);
    % %    coupling.FactorAmplitude = Q*aAmplitude*(-inv(Q)*inv(Mii)*Mis);
    %     shield.Coupling = shield.C*coupling.Factor;
    %     shield.CouplingSin = shield.C*coupling.FactorSin;
    %     shield.CouplingCos = shield.C*coupling.FactorCos;
    % %    shield.CouplingAmplitude = shield.C*coupling.FactorAmplitude;
    % 
    %     s_induce_reduced = coupling.Factor*coil.s_reduced;
    %     s_induce_reducedCos = coupling.FactorCos*coil.s_reduced;
    %     s_induce_reducedSin = coupling.FactorSin*coil.s_reduced;
    % %    s_induce_reducedAmplitude =coupling.FactorAmplitude*coil.s_reduced;
    % 
    %     shield.s_induce = retrieveCurrentVector(s_induce_reduced,shield.nonReducedSize,shield.nbrBorderNode);
    %     shield.s_induceCos = retrieveCurrentVector(s_induce_reducedCos,shield.nonReducedSize,shield.nbrBorderNode);
    %     shield.s_induceSin = retrieveCurrentVector(s_induce_reducedSin,shield.nonReducedSize,shield.nbrBorderNode);
    % %    shield.s_induceAmplitude = retrieveCurrentVector(s_induce_reducedAmplitude,shield.nonReducedSize,shield.nbrBorderNode);
    % 
    % %     if size(coil.btarget,2) == size(coil.Cfull,1)
    % %         % When this is true, we just one field components, so we have to use
    % %         % another C
    % %         shield.CwholeField = [shield.Cx;shield.Cy;shield.Cz];
    % %     end
    % 
    %     shield.B = shield.CwholeField*shield.s_induce;
    %     shield.Bcos = shield.CwholeField*shield.s_induceCos;
    %     shield.Bsin = shield.CwholeField*shield.s_induceSin;
    % 
    % %Whe have to split the field
    %     nbrPoint = size(rk,1);
    %     Bx = shield.B(1:nbrPoint);
    %     By = shield.B(nbrPoint+1:2*nbrPoint);
    %     Bz = shield.B(2*nbrPoint+1:end);
    %     [shield.bc,shield.bs] = getSphericalHarmonicsCoefficientMeasure2(Bx,By,Bz,sqrt(nbrPoint),rk);
    %     shield.amplitude = shield.bc(2).coefficient(1,1);
    % 
    %     nbrPoint = size(rk,1);
    %     Bx = shield.Bcos(1:nbrPoint);
    %     By = shield.Bcos(nbrPoint+1:2*nbrPoint);
    %     Bz = shield.Bcos(2*nbrPoint+1:end);
    %     [shield.bcCos,shield.bsCos] = getSphericalHarmonicsCoefficientMeasure2(Bx,By,Bz,sqrt(nbrPoint),rk);
    % 
    %     Bx = shield.Bsin(1:nbrPoint);
    %     By = shield.Bsin(nbrPoint+1:2*nbrPoint);
    %     Bz = shield.Bsin(2*nbrPoint+1:end);
    %     [shield.bcSin,shield.bsSin] = getSphericalHarmonicsCoefficientMeasure2(Bx,By,Bz,sqrt(nbrPoint),rk);
    % 
    %     shield.phase = atan2(shield.bcCos(2).coefficient(1,1),shield.bcSin(2).coefficient(1,1)); % form A*sin(wt+phi)
    %     shield.amplitudebis = sqrt(shield.bcSin(2).coefficient(1,1)^2+shield.bcCos(2).coefficient(1,1)^2);
    % 
    %     amplitudeShield(ii) = shield.amplitudebis;
    %     amp(ii) = shield.amplitude;
    %     amplitudeCoil(ii) = coil.amplitude;
    %     phase(ii) = shield.phase;
    % end
    % toc
    % figure
    % plot(frequencyArray,amplitudeShield)

    %% Check the basis function
    % 
    % % 
    shield.basis = basisFunction3(shield.node, shield.triangle, shield.nodeLinkToTriangle,shield.center);
    displayBasis(shield.node,shield.basis,shield.center)
    xlim([-0.3 0.3])
    ylim([-0.3 0.3])
    zlim([-0.3 0.3])
    displayBasis(coil.node,coil.basis,coil.center)
    % 
    for i=1:size(shield.basis,1)
        if shield.node(i).nbrTriangle==32
            i
        end
    end

    %% test with the stream function and the sphere
    shield.s_test = zeros(2557,1);

    shield.s_test(1) = 1;

    displayStreamFunction(shield.listTriangle,shield.s_test,shield.listNode)
    [jTri,jAbsTri,j,jAbs] = calculateCurrentDensity(shield.s_test,shield.node,shield.triangle,shield.basis);
    displayCurrentDensity(shield.listTriangle,shield.listNode,jTri,jAbsTri)


     tt = 0; tt2=0;
     for i=1:size(shield.L,1)
         tt(i) = shield.L(i,i);
         tt2(i) = shield.R(i,i);
     end
     for i=1:size(shield.node,2)
         tt3(i) = sum([shield.triangle(shield.node(i).linkToTriangle).air]);
     end
      figure;plot(tt3);title('air of triangle linked to the node')
      figure;plot([shield.triangle(:).air]');title('air of triangle')
     figure;plot(tt);title('Inductance');figure;plot(tt2);title('resistance')
end
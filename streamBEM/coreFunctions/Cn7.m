function [Cx,Cy,Cz] = Cn7(node,triangle,basis,r)

% This file aim at calculating the so called "Cn matrix"
%
%
% node : a matrix with the 3d position of each node (in meter)
% triangle : a matrix linking 3 node together to form a triangle
% basis : structure with the basis vector and the related informations
% r : the target point

%ChangeLog
% v5: adapted to the new structure of node and triangle
% v7: clean the code

tic;


[u,v,ck] = triGaussPoints(2);
for i=1:size(u,1)
    w(i) = 1-u(i)-v(i);
end


mu_0 = 4*pi*10^-7;
coef = mu_0/(4*pi);
nbrIntegrationPoints = size(ck,1);

Cx = zeros(size(r,1),size(node,2));
Cy = zeros(size(r,1),size(node,2));
Cz = zeros(size(r,1),size(node,2));

%activate the parallel function if available
if license('test','Distrib_Computing_Toolbox')
    matlabVersion = version;
    matlabVersion = str2num(matlabVersion(1:3));
    if matlabVersion < 8.2
        [TF,~] = license('checkout', 'Distrib_Computing_Toolbox');
        if TF
            schd = findResource('scheduler', 'configuration', 'local');
            numWorkers = schd.ClusterSize;
        end

        if matlabpool('size') == 0  && TF && numWorkers >1
            % checking to see if the pool is already open and of we have the licence
            % and at least 2 cores
            matlabpool open
        end
    elseif matlabVersion >= 8.2 
        poolobj = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(poolobj)
            parpool;
        end
    end
end

fprintf('Loop n=%5.0i',0);
for n=1:size(node,2) % For a given node
    tempCx = zeros(1,size(r,1));
    tempCy = zeros(1,size(r,1));
    tempCz = zeros(1,size(r,1));
    for k=1:size(r,1); % And for a given target point
        bigSumCx = 0;
        bigSumCy = 0;
        bigSumCz = 0;
        for i=1:node(n).nbrTriangle; %We select the triangle that link to this node
            currentTriangle_i = node(n).linkToTriangle(i);
            r_x = r(k,1); %X coordinate of the provided point
            r_y = r(k,2); %Y coordinate of the provided point
            r_z = r(k,3); %Z coordinate of the provided point

            vni_x = basis(n).triangle(i).value(1);
            vni_y = basis(n).triangle(i).value(2);
            vni_z = basis(n).triangle(i).value(3);
            
            rpk_x(:,1) = basis(n).triangle(i).r_o(:,1);
            rpk_y(:,1) = basis(n).triangle(i).r_o(:,2);
            rpk_z(:,1) = basis(n).triangle(i).r_o(:,3);
            for l=1:nbrIntegrationPoints

                normRRp = sqrt((r_x-rpk_x(l))^2+(r_y-rpk_y(l))^2+(r_z-rpk_z(l))^2);
                %fC = ck(l)*cross(-(r(k,1:3)-rp(l,:))/normRRp^3,vni);% To
                %optimise the executzion time, we do not use those functions
                fCx(l) = ck(l)*((-vni_z*(r_y-rpk_y(l))+(vni_y*(r_z-rpk_z(l))))/(normRRp^3));
                fCy(l) = ck(l)*((-vni_x*(r_z-rpk_z(l))+(vni_z*(r_x-rpk_x(l))))/(normRRp^3));
                fCz(l) = ck(l)*((-vni_y*(r_x-rpk_x(l))+(vni_x*(r_y-rpk_y(l))))/(normRRp^3));
            end

            sum1Cx = sum(fCx);
            bigSumCx = bigSumCx + sum1Cx*(2*triangle(currentTriangle_i).air);
            %bigSumCx = bigSumCx + scallingToCompensateTheTransform1*sum1Cx;
            
            sum1Cy = sum(fCy);
            bigSumCy = bigSumCy + sum1Cy*(2*triangle(currentTriangle_i).air);
            %bigSumCy = bigSumCy + scallingToCompensateTheTransform1*sum1Cy;
            
            sum1Cz = sum(fCz);
            bigSumCz = bigSumCz + sum1Cz*(2*triangle(currentTriangle_i).air);
            %bigSumCz = bigSumCz + scallingToCompensateTheTransform1*sum1Cz;
        end
        tempCx(k) = coef*bigSumCx;
        tempCy(k) = coef*bigSumCy;
        tempCz(k) = coef*bigSumCz;
        %C(k,n) = coef*bigSum;
    end
    Cx(:,n) = tempCx;
    Cy(:,n) = tempCy;
    Cz(:,n) = tempCz;
    fprintf(1,'\b\b\b\b\b%5.0i',n);
end
%fprintf(1,' - Done \n');
fprintf(' - Done in %5.0f sec.\n',toc);
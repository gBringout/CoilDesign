function [AX,AY,AZ] = Amn2(node_1, triangle_1, basis_1,node_2, triangle_2, basis_2)

% This file aim at calculating the magnetic vector potential
%
% node : a matrix with the 3d position of each node (in meter)
% triangle : a matrix linking 3 node together to form a triangle
% basis : a matrix with the description of the basis function
%
% changelog : based on Lmn11 and Cn6, for a cleaner version using the basis structure



[~,~,ck] = triGaussPoints(2);

%%
mu_0 = 4*pi*10^-7;
coef = mu_0/(4*pi);
nbrIntegrationPoints = size(ck,1);

AX = zeros(size(node_1,2),size(node_2,2));
AY = zeros(size(node_1,2),size(node_2,2));
AZ = zeros(size(node_1,2),size(node_2,2));

[TF,~] = license('checkout', 'Distrib_Computing_Toolbox');
if matlabpool('size') == 0  && TF
    % checking to see if the pool is already open and of we have the licence
    matlabpool open
end

tic

parfor m=1:size(node_1,2); %For every node on the induced surface
    tempX = zeros(1,size(node_2,2));
    tempY = zeros(1,size(node_2,2));
    tempZ = zeros(1,size(node_2,2));
    r = node_1(m).coord(:);

    for n=1:size(node_2,2); %and every node on the coil surface
        superBigSumX = 0;
        superBigSumY = 0;
        superBigSumZ = 0;
        
        for i=1:node_2(n).nbrTriangle; %Number of face in this face
            currentTriangle_i = node_2(n).linkToTriangle(i);
            integral = 0;
            vmi = basis_2(n).triangle(i).value;
            r_o = basis_2(n).triangle(i).r_o;

            for k = 1:nbrIntegrationPoints
                integral = integral + ck(k)/norm(r'-r_o(k,:));
            end

            superBigSumX = superBigSumX+vmi(1)*integral*2*triangle_2(currentTriangle_i).air;
            superBigSumY = superBigSumY+vmi(2)*integral*2*triangle_2(currentTriangle_i).air;
            superBigSumZ = superBigSumZ+vmi(3)*integral*2*triangle_2(currentTriangle_i).air;
        end
        tempX(n) = coef*superBigSumX;
        tempY(n) = coef*superBigSumY;
        tempZ(n) = coef*superBigSumZ;
    end
    AX(m,:) = tempX;
    AY(m,:) = tempY;
    AZ(m,:) = tempZ;
end

fprintf(' - Done in %5.0f sec.\n',toc);
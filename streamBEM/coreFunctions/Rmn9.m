function [R] = Rmn9(node, triangle,basis,rho,t)

% This file aim at calculating the so called "Rmn matrix"
%
% Equation come from the thesis of Michael Poole "Improved Equipment and
% Techniques for Dynamic Shimming in High Field MRI
% page 66
%
% Node : a matrix with the 3d position of each Node (in meter)
% Triangle : a matrix linking 3 Node together to form a triangle
% nodeLinkToTriangle : structure with the triangle linked to the same Node
%
% Changelog :
% Based on Rmn6 with an exterbal basis function
% Based on Rmn7 with the use of the new structure of node and triangle
% Based on Rmn8 with the use of a new structure for the basis

tic;
R = zeros(size(node,2),size(node,2));

conductorSurfaceResistance = rho/t;
if matlabpool('size') == 0 % checking to see if my pool is already open
   matlabpool open 4
end

fprintf('Loop n=%5.0i',0);
%for m=1:size(node,1);
parfor m=1:size(node,2); % m is a node
    temp = zeros(1,size(node,2));
    for n=m:size(node,2); %n is a node
        bigSum = 0;
        for i=1:node(m).nbrTriangle
            currentTriangle_m = node(m).linkToTriangle(i);
            for j=1:node(n).nbrTriangle
                currentTriangle_n = node(n).linkToTriangle(j);
				
                % if the triangle share a common node
                if currentTriangle_m == currentTriangle_n
					vmi = basis(m).triangle(i).value;
					vnj = basis(n).triangle(j).value; %#ok<PFBNS>
                    bigSum = bigSum + dot(vmi,vnj)*triangle(m).air;
                end
            end
        end

        temp(n) = conductorSurfaceResistance*bigSum;
    end
    R(m,:) = temp;
    %fprintf('Loop m=%i \n',m);
    %fprintf(1,'\b\b\b\b\b%5.0i',m);
end


% the Rmn operator is symetric
for m=2:size(node,2) % m is a node
    for n=1:m-1 %n is a node
        R(m,n) = R(n,m);
    end
end
%fprintf(1,' - Done \n');
fprintf(' - Done in %5.0f sec.\n',toc);
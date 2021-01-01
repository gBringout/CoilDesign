function [Lwp] = Laplacian3(node, triangle,h)   
% This file aim at calculating the Laplacian of a mesh
% It has been implemented only for a flat mesh in the XZ plan
%
%
% node : a matrix with the 3d position of each node (in meter)
% triangle : a matrix linking 3 node together to form a triangle
% h a optional parameter linked to the mesh properties. example: 0.0088


% if h is note provided, it is estimated the average length between to point of the mesh
if nargin<3
    nbrSegment = 0;
    for t=1:size(triangle,2)
        nbrSegment = nbrSegment + 1;
        a = node(triangle(t).node(1)).coord; %first point of the triangle
        nbrSegment = nbrSegment + 1;
        b = node(triangle(t).node(2)).coord; %second point of the triangle
        nbrSegment = nbrSegment + 1;
        c = node(triangle(t).node(3)).coord; %third point of the triangle
        lengthSide(nbrSegment-2) = norm(a-b);
        lengthSide(nbrSegment-1) = norm(b-c);
        lengthSide(nbrSegment) = norm(c-a);
    end
    h = mean(lengthSide)/3;
end
constante = 1/(4*pi*h^2);

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

%%
dim1Lwp = size(node,2);
dim2Lwp = size(node,2);
Lwp = zeros(dim1Lwp,dim2Lwp);
tic
parfor w=1:dim1Lwp % for each node
    temp = zeros(1,dim2Lwp);
    sumW = 0;
    for p=1:dim2Lwp % for each node
        sumP =0;
        factor = (exp(-norm(node(p).coord-node(w).coord)^2/(4*h)));
        if node(p).coord(:,2) == node(w).coord(:,2) %Only for open coils; It means when we have the same sub surface
            for t=1:size(node(p).linkToTriangle,1)
                sumP = sumP + triangle(t).air/3;%we use triangular mesh, so there is always 3 point for each face
            end
        end
        temp(p) = sumP*factor;
    end
    for t=1:size(triangle,2)
        p1 = node(triangle(t).node(1)).coord;%Triangular mesh, there is just 3 face, 3 point
        p2 = node(triangle(t).node(2)).coord;
        p3 = node(triangle(t).node(3)).coord;
        factor1 = (exp(-norm(p1-node(w).coord)^2/(4*h)));
        factor2 = (exp(-norm(p2-node(w).coord)^2/(4*h)));
        factor3 = (exp(-norm(p3-node(w).coord)^2/(4*h)));
        if p1(2) == node(w).coord(2) %Only for open coils; It means when we have the same sub surface
            sumW = sumW + triangle(t).air*(factor1+factor2+factor3);
        end
    end
    temp(w) = temp(w)-sumW;
    Lwp(w,:) = constante*temp(:);
end
fprintf(' - Done in %5.0f sec.\n',toc);
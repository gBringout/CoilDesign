function [triangleNew,nodeNew] = processMesh(triangle,node)
% Process the mesh to calculate the air of triangle, the barycenter, which
% triangle is connected to which node, etc.

tic
%% Value useful for the Gauss-Legendre integration
[u,v,~] = triGaussPoints(2);
w = zeros(size(u,1),1);
for i=1:size(u,1)
    w(i) = 1-u(i)-v(i);
end

%% Process the triangle
for i=1:size(triangle,1)
    % linked to node
    triangleNew(i).node = triangle(i,:);
    Anode = triangle(i,1);
	A = node(Anode,:);
    
    Bnode = triangle(i,2);
	B = node(Bnode,:);
    
    Cnode = triangle(i,3);
	C = node(Cnode,:);
    
    %calculate the air, according to: http://geomalgorithms.com/a01-_area.html
	triangleNew(i).air = (1/2)*norm(cross(A-C,B-C)); 
    
    %use the signed triangle air to find the correct normal direction
    % and correct the order of the points
    if triangleNew(i).air>0
        triangleNew(i).normal = cross(A-C,B-C)/norm(cross(A-C,B-C));
        triangleNew(i).node(1) = Cnode; % Sort the points in the right order
        triangleNew(i).node(2) = Anode;
        triangleNew(i).node(3) = Bnode;
    elseif  triangleNew(i).air<0
        triangleNew(i).air = -triangleNew(i).air;
        triangleNew(i).normal = cross(B-C,A-C)/norm(cross(B-C,A-C));
        triangleNew(i).node(1) = Cnode; % Sort the points in the right order
        triangleNew(i).node(2) = Bnode;
        triangleNew(i).node(3) = Anode;
        Atemp = A;
        A = B;
        B = Atemp;
    else
        disp('triangle with an area of zero, this is wierd')
    end

    % calculate the center as the barycenter
    %see: http://mathforum.org/library/drmath/view/54899.html
	triangleNew(i).center = (A+B+C)/3;
    triangleNew(i).r_o = changtRef(A, B, C,u, v, w);
end

%% Process the node
for i=1:size(node,1)
    %copy the coordinate
    nodeNew(i).coord = node(i,:);
    
    [nodeNew(i).linkToTriangle, ~] = find(triangle==i); % Find all the face contaigning the i_th node
    nodeNew(i).nbrTriangle = size(nodeNew(i).linkToTriangle,1);
end
fprintf(' - Done in %5.0f sec.\n',toc);
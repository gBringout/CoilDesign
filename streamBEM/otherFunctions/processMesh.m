function [triangleNew,nodeNew] = processMesh(triangle,node)


for i=1:size(triangle,1)
    % linked to node
    triangleNew(i).node = triangle(i,:);
    
%     % air
%     % according to wikipedia page : http://fr.wikipedia.org/wiki/Aire_d%27un_triangle
%     Ax = node(triangle(i,1),1);
% 	Ay = node(triangle(i,1),2);
% 	Az = node(triangle(i,1),3);
% 	Bx = node(triangle(i,2),1);
% 	By = node(triangle(i,2),2);
% 	Bz = node(triangle(i,2),3);
% 	Cx = node(triangle(i,3),1);
% 	Cy = node(triangle(i,3),2);
% 	Cz = node(triangle(i,3),3);
%     
% 	triangleNew(i).air = (1/2)*sqrt((det([[Cx Ax Bx];[Cy Ay By];[1 1 1]]))^2+(det([[Cy Ay By];[Cz Az Bz];[1 1 1]]))^2+(det(det([[Cz Az Bz];[Cx Ax Bx];[1 1 1]])))^2); 
%     
    
	A = node(triangle(i,1),:);
	B = node(triangle(i,2),:);
	C = node(triangle(i,3),:);
    
    %calculate the air, according to: http://geomalgorithms.com/a01-_area.html
	triangleNew(i).air = abs((1/2)*norm(cross(A-C,B-C))); 
    
    
    % calculate the center
    %see: http://mathforum.org/library/drmath/view/54899.html
	triangleNew(i).center = (A+B+C)/3;
end

for i=1:size(node,1)
    %copy the coordinate
    nodeNew(i).coord = node(i,:);
    
    [nodeNew(i).linkToTriangle, ~] = find(triangle==i); % Find all the face contaigning the i_th node
    nodeNew(i).nbrTriangle = size(nodeNew(i).linkToTriangle,1);
end
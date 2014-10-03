function [P] = changtRef(P1, P2, P3,u, v, w)

% This file aim at calculating the position in the cartesians coordinate
% (x,y,z) over a triangle,
% from point in a generalized coordinate system (u,v,w). Note that w is
% genereted using the shape function associated with a linear interpolation
% over a unit isocele triangle.
% The point P3 is the point with coordinate (0,0)

%
%changelog : 

P=zeros(size(u,1),3);
for i=1:size(u,1)
    P(i,:) = P1*u(i)+P2*v(i)+P3*w(i); % We change the coordinate from (x,y,z) to (u,v) (page 71)
end
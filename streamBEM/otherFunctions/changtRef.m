function [P] = changtRef(P1, P2, P3,u, v, w)
% This file aim at calculating the changement of coordinates for a triangle.

% u,v,w: coordinate on the unit triangle
% P1, P2, P3: coordinate in the global coordinate system. Note that the
% point P3 is the one maped with the point (0,0)

P=zeros(size(u,1),3);
for i=1:size(u,1)
    P(i,:) = P1*u(i)+P2*v(i)+P3*w(i);
end
function [r_o] = changtRef(A, B, C,u, v, w)

% This file aim at calculating the rmi for the parfor loop of Lmn9
% That is, from the point u,v on the unitary triangle, we can calculate
% back the coordinate (x,y,z) in the other referential

%
%changelog : 

r_o=zeros(size(u,1),3);
for i=1:size(u,1)
    r_o(i,:) = B*u(i)+A*v(i)+C*w(i); % We change the coordinate from (x,y,z) to (u,v) (page 71)
end
function [listNode2] = translateMesh(listNode,translateX,translateY,translateZ)
% translate the list of point using the 3 distance in meter
% Only the first three coordinate are taken into account
% the order is taken as
%1: x
%2: y
%3: z

listNode2(:,1) = listNode(:,1)+translateX;
listNode2(:,2) = listNode(:,2)+translateY;
listNode2(:,3) = listNode(:,3)+translateZ;
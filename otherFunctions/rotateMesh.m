function [listNode2] = rotateMesh(listNode,angleX,angleY,angleZ)
% rotate the list of point using the 3 angle in degree
% the rotation is always around the defined axe
% Start with the z axe, then the y and finally the x one.
% for example: angleX around the X axe
% Only the first three coordinate are taken into account
% the order is taken as
%1: x
%2: y
%3: z

angleRadianX = angleX*(2*pi)/360;
angleRadianY = angleY*(2*pi)/360;
angleRadianZ = angleZ*(2*pi)/360;

RX =    [[1 0 0]...
        ;[0 cos(angleRadianX) sin(angleRadianX)]...
        ;[0 -sin(angleRadianX) cos(angleRadianX)]];
    
RY =    [[cos(angleRadianY) 0 -sin(angleRadianY)]...
        ;[0 1 0]...
        ;[sin(angleRadianY) 0 cos(angleRadianY)]];
    
RZ =    [[cos(angleRadianZ) sin(angleRadianZ) 0]...
        ;[-sin(angleRadianZ) cos(angleRadianZ) 0]...
        ;[0 0 1]];
R = RX*RY*RZ;
% R = [[cos(angleRadianY)*cos(angleRadianZ) 0 0]...
%     ;[0 cos(angleRadianX)*cos(angleRadianZ) 0]...
%     ;[0 0 cos(angleRadianY)*cos(angleRadianX)]];
listNode2 = zeros(size(listNode,1),3);
for i=1:size(listNode,1)
    listNode2(i,:) = R*(listNode(i,:)');
end
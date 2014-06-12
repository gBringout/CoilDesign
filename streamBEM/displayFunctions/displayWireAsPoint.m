function [] = displayWireAsPoint(Wire)


figure
hold all;
for i = 1 :size(Wire,2)
    if Wire(i).currentDirection == 1
        plot3(Wire(i).Coord(1,:),Wire(i).Coord(2,:),Wire(i).Coord(3,:),'b*');
    elseif Wire(i).currentDirection == -1
        plot3(Wire(i).Coord(1,:),Wire(i).Coord(2,:),Wire(i).Coord(3,:),'r*');
    end
end
%view([0 0]);
xlabel('X axes');
ylabel('Y axes');
zlabel('Z axes');
axis equal;
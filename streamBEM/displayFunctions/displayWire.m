function [] = displayWire(wire,addToActualFigure)

if nargin < 2
  addToActualFigure = 0;
end

if addToActualFigure == 0
    figure
end

hold all;
for i = 1 :size(wire,2)
    if wire(i).currentDirection == 1
        plot3(wire(i).Coord(1,:),wire(i).Coord(2,:),wire(i).Coord(3,:),'blue');
    else %if wire(i).currentDirection == -1
        plot3(wire(i).Coord(1,:),wire(i).Coord(2,:),wire(i).Coord(3,:),'red');
    end
end
view([0 0]);
xlabel('X axis /m');
ylabel('Y axis /m');
zlabel('Z axis /m');
axis equal;
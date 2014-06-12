function [] = displayWireSubFigure(wire)


figure
index = ceil(sqrt(size(wire,2)));

for i = 1 :size(wire,2)
    subplot(index,index,i)
    if wire(i).currentDirection == 1
        plot3(wire(i).Coord(1,:),wire(i).Coord(2,:),wire(i).Coord(3,:),'b*');
    else %if wire(i).currentDirection == -1
        plot3(wire(i).Coord(1,:),wire(i).Coord(2,:),wire(i).Coord(3,:),'r*');
    end
    xlim([-0.3 0.3]);
    ylim([-0.3 0.3]);
    zlim([-0.3 0.3]);
    title(sprintf('%i',i))
    axis square
end

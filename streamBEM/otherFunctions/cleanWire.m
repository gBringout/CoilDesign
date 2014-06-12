function [wire] = cleanWire(wire)


for i = 1 :size(wire,2)
    j = 1;
    while j<=size(wire(i).Coord,2)-1
        if (wire(i).Coord(1,j) == wire(i).Coord(1,j+1) &&...
                wire(i).Coord(2,j) == wire(i).Coord(2,j+1) &&...
                wire(i).Coord(3,j) == wire(i).Coord(3,j+1)) || ...
                ((abs( norm(wire(i).Coord(:,j)-wire(i).Coord(:,j+1),2))<0.003))
            %if two point are the same, one after the other
            %or if the distance between to adjacente point is lower than 3
            %mm
            wire(i).Coord(:,j+1) = []; %delete the second point
        else
            j = j+1;
        end
    end
end

% Then we look if we have wire with NaN as coordinate
for i = size(wire,2):-1:1
    if isnan(wire(i).Coord(1,1))
        wire(i) = [];
    end
end
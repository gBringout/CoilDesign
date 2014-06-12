function [wire] = closeWire(wire)
% Close the loop by copying the first point to the last position


for i = 1 :size(wire,2)
    t = size(wire(i).Coord,2);
    wire(i).Coord(:,t+1) =  wire(i).Coord(:,1);
end
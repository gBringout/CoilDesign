function [wire] = makeSolenoid(center, radius, nbrSegment,direction,plan)

if nargin < 5
    center = [[0 0 0]' [0 0 0]' [0 0 0]'];
    radius  = [0.2 0.3 0.25];
    nbrSegment  = [20 30 100];
    direction  = [1 1 -1];
    plan = 'xy';
    disp('Error,not enough parameter')
end

% nPoint = 1;
% wire.Coord(1,nPoint) = x;
% wire.Coord(2,nPoint) = y;
% wire.Coord(3,nPoint) = z;
if any([size(center,2) ~= size(radius,2)...
        size(center,2) ~= size(nbrSegment,2)...
        size(center,2) ~= size(direction,2)])
    disp('Error, Center, Radius, nbrSegment and direction have to have the same size')
else
    if strcmp(plan,'xy')
        for i=1:size(center,2)
            nbrOfPoint = 1;
            wire(i).currentDirection = direction(i);
            for j=1:nbrSegment(i)+1
                angle = (j-1)*2*pi/nbrSegment(i);
                wire(i).Coord(1,nbrOfPoint) = radius(i)*cos(angle) + center(1,i);
                wire(i).Coord(2,nbrOfPoint) = radius(i)*sin(angle) + center(2,i);
                wire(i).Coord(3,nbrOfPoint) = center(3,i);
                nbrOfPoint = nbrOfPoint+1;
            end
        end
    end
end


%displayWire(wire);

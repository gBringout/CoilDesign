function [wire] = makeLoop(center, radius, nbrSegment,direction,plan)
% Implement a function to make a planar wire loop at a given 'center', with a given 'radius', made of a given number of straight segment 'nbrSegment', with the current flowing in a direction. Many loops can be create with a single function call by increasing the dimension of all inputs accordingly.

% We have to provide
% center : [m m m] the center position in 3D, as x y z
% radius : [m] the radius of the loop
% nbrSegment : [ ] the number of straight segment used to make each loop
% direction : [ ] (-1 or 1) direction of the flowing current. 1 = counter clock wise
% plan : [] plan on which the loop is made.

% output
% wire : structure with the wire x,y,z position which form a loop and the direction in which the current should flow
%
% created by Gael Bringout - 04.02.2015


if nargin < 5
    center = [[0 0 0]' [0 0 0]' [0 0 0]'];
    radius  = [0.2 0.3 0.25];
    nbrSegment  = [20 30 100];
    direction  = [1 1 -1];
    plan = 'xy';
    error('Error,not enough parameter')
end

% nPoint = 1;
% wire.Coord(1,nPoint) = x;
% wire.Coord(2,nPoint) = y;
% wire.Coord(3,nPoint) = z;
if any([size(center,2) ~= size(radius,2)...
        size(center,2) ~= size(nbrSegment,2)...
        size(center,2) ~= size(direction,2)])
    error('Error, Center, Radius, nbrSegment and direction have to have the same size')
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

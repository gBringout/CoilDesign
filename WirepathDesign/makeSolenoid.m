function [wire] = makeSolenoid(center, radius, nbrSegment,direction,plan)
% Implement a function to make a planar wire loop at a given 'center', with a given 'radius', made of a given number of straight segment 'nbrSegment', with the current flowing in a direction. Many loops can be create with a single function call by increasing the dimension of all inputs accordingly. This is the old name of the function. This will stay like that to ensure comptatibility

% We have to provide
% center : [m m m] the center position in 3D, as x y z
% radius : [m] the radius of the loop
% nbrSegment : [ ] the number of straight segment used to make each loop
% direction : [ ] (-1 or 1) direction of the flowing current. 1 = counter clock wise
% plan : [] plan on which the loop is made.

% output
% wire : structure with the wire x,y,z position which form a loop and the direction in which the current should flow
%
% created by Gael Bringout - 03.10.2013
% updated according to Ankit comments by Gael Bringout - 04.02.2015

disp('function is depreciated, please use makeLoop instead')

wire = makeLoop(center, radius, nbrSegment,direction,plan);
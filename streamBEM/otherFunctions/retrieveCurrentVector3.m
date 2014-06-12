function [i] = retrieveCurrentVector3(x_lambda,subBoundaries1,optimizationType)


nbrBorderNode = size(subBoundaries1,1);

if nbrBorderNode == 1
    nbrNodeOnBoundary1 = size(subBoundaries1(1).node,1);
    Ip = x_lambda;
    I0 = Ip(2:end);
    if strcmp(optimizationType,'QP')
        I1(1:nbrNodeOnBoundary1,1) = Ip(1);
    else
        I1(1:nbrNodeOnBoundary1,1) = 0;
    end
    I = [I1;I0];
    i = I;
elseif nbrBorderNode == 2
    nbrNodeOnBoundary1 = size(subBoundaries1(1).node,1);
    nbrNodeOnBoundary2 = size(subBoundaries1(2).node,1);
    Ip = x_lambda;
    I0 = Ip(3:end);
    if strcmp(optimizationType,'QP')
        I1(1:nbrNodeOnBoundary2,1) = Ip(2);
        I2(1:nbrNodeOnBoundary1,1) = Ip(1);
    else
        I1(1:nbrNodeOnBoundary2,1) = 0;
        I2(1:nbrNodeOnBoundary1,1) = 0;
    end
    I = [I2;I1;I0];
    i = I;
else
    disp('error : no function to reduce matrix with 3 sub-boundaries');
end
function [MR] = ReduceCMatrix4(M,subBoundaries)

nbrBorderNode = size(subBoundaries,1);

if nbrBorderNode == 1
    nonReducedSize = size(M,2);
    nbrNodeOnBoundary = size(subBoundaries(1).node,1);
    
    M0 = M(:,nbrNodeOnBoundary+1:nonReducedSize);
    M1 = M(:,1);
    MR = [M1 M0];
elseif nbrBorderNode == 2
    nonReducedSize = size(M,2);
    nbrNodeOnBoundary1 = size(subBoundaries(1).node,1);
    nbrNodeOnBoundary2 = size(subBoundaries(2).node,1);
    totalNbrNodeOnBoundaries = nbrNodeOnBoundary1+nbrNodeOnBoundary2;
    
    M0 = M(:,totalNbrNodeOnBoundaries+1:end);
    M1 = M(:,nbrNodeOnBoundary1+nbrNodeOnBoundary2);
    M2 = M(:,nbrNodeOnBoundary1);
    MR = [M2 M1 M0];     
else
    disp('error : no function to reduce matrix with 3 sub-boundaries');
end
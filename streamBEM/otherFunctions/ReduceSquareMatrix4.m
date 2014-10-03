function [MR] = ReduceSquareMatrix4(M,subBoundaries1)
%we suppose that all the border are on the top of the matrix
% we use the new version of the border
% subBoundaries1 = coil.subBoundaries1 ;
% and correct the limitation to squarre Matrix

nbrBorderNode = size(subBoundaries1,1);

if nbrBorderNode == 1
    nonReducedSize = size(M);
    nbrNodeOnBoundary = size(subBoundaries1(1).node,1);
    M0 = M(nbrNodeOnBoundary+1:nonReducedSize(1),nbrNodeOnBoundary+1:nonReducedSize(2));
    M1 = M(1,nbrNodeOnBoundary+1:nonReducedSize(2));
    M2 = M(1,1);
    M3 = M(nbrNodeOnBoundary+1:nonReducedSize(1),1);

    MR = [[M2 M1];[M3 M0]];
elseif nbrBorderNode == 2
    nonReducedSize = size(M);
    nbrNodeOnBoundary1 = size(subBoundaries1(1).node,1);
    nbrNodeOnBoundary2 = size(subBoundaries1(2).node,1);
    totalNbrNodeOnBoundaries = nbrNodeOnBoundary1+nbrNodeOnBoundary2;
    % Here we replace the whole equation just by one. As stated by Peeren
    % on page 63, 'all the value are the same' so we can just keep one
    % equation per boundary condition
    
    M1 = M(nbrNodeOnBoundary1,totalNbrNodeOnBoundaries+1:end);
    M2 = M(nbrNodeOnBoundary1,nbrNodeOnBoundary1+nbrNodeOnBoundary2);
    M3 = M(nbrNodeOnBoundary1,nbrNodeOnBoundary1);
    
    M4 = M(nbrNodeOnBoundary1+nbrNodeOnBoundary2,totalNbrNodeOnBoundaries+1:end);
    M5 = M(nbrNodeOnBoundary1+nbrNodeOnBoundary2,nbrNodeOnBoundary1+nbrNodeOnBoundary2);
    M6 = M(nbrNodeOnBoundary1+nbrNodeOnBoundary2,nbrNodeOnBoundary1);
    
    M0 = M(totalNbrNodeOnBoundaries+1:end,totalNbrNodeOnBoundaries+1:end);
    M7 = M(totalNbrNodeOnBoundaries+1:end,nbrNodeOnBoundary1+nbrNodeOnBoundary2);
    M8 = M(totalNbrNodeOnBoundaries+1:end,nbrNodeOnBoundary1);
    
    MR = [[M3 M2 M1];[M6 M5 M4];[M8 M7 M0]];
else
    disp('error : no function to reduce matrix with 3 or more sub-boundaries');
end
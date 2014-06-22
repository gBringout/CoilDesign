function [MR] = ReduceLUpMatrix4(M,subBoundariesVertical,subBoundariesHorizontal)
% used as ReduceLUpMatrix2(coupling.LDownFull,shield.subBoundaries,coil.subBoundaries)

% nbrBorderHorizontal is the column direction %20 of 320
% nbrBorderVertical is the line direction %20 of 420
% M = L
nbrBorderHorizontal = size(subBoundariesHorizontal,1);
nbrBorderVertical = size(subBoundariesVertical,1);
if nbrBorderVertical == 2 && nbrBorderHorizontal == 2 %case of the cylinder
    %if nbrBorderVertical~=nbrBorderVertical || nbrBorderHorizontal(1)~=nbrBorderHorizontal(2)
    %    disp('error: the boundary have to have the same number of elements')
    %end
        
    nonReducedSize = size(M);
    nbrNodeOnBoundaryV1 = size(subBoundariesVertical(1).node,1);
    nbrNodeOnBoundaryV2 = size(subBoundariesVertical(2).node,1);
    nbrNodeOnBoundaryV = nbrNodeOnBoundaryV1 + nbrNodeOnBoundaryV2;
    nbrNodeOnBoundaryH1 = size(subBoundariesHorizontal(1).node,1);
    nbrNodeOnBoundaryH2 = size(subBoundariesHorizontal(2).node,1);
    nbrNodeOnBoundaryH = nbrNodeOnBoundaryH1 + nbrNodeOnBoundaryH2;
    totalNbrNodeOnBoundaries = nbrNodeOnBoundaryV+nbrNodeOnBoundaryH;
    % Here we replace the whole equation just by one. As stated by Peeren
    % on page 63, 'all the value are the same' so we can just keep one
    % equation per boundary condition

    %M1 = sum(M(1:95,190+1:2470),1);
    M1 = M(nbrNodeOnBoundaryV1,nbrNodeOnBoundaryH1 + nbrNodeOnBoundaryH2 + 1:end); %ok
    M2 = M(nbrNodeOnBoundaryV1,nbrNodeOnBoundaryH1 + nbrNodeOnBoundaryH2); %ok
    M3 = M(nbrNodeOnBoundaryV1,nbrNodeOnBoundaryH1); %ok
    
    M4 = M(nbrNodeOnBoundaryV1+nbrNodeOnBoundaryV2,nbrNodeOnBoundaryH1 + nbrNodeOnBoundaryH2 + 1:end); %ok
    M5 = M(nbrNodeOnBoundaryV1+nbrNodeOnBoundaryV2,nbrNodeOnBoundaryH1 + nbrNodeOnBoundaryH2); %ok
    M6 = M(nbrNodeOnBoundaryV1+nbrNodeOnBoundaryV2,nbrNodeOnBoundaryH1); %ok
    
    
    M0 = M(nbrNodeOnBoundaryV1+nbrNodeOnBoundaryV2+1:end,nbrNodeOnBoundaryH1 + nbrNodeOnBoundaryH2 + 1:end);
    M7 = M(nbrNodeOnBoundaryV1+nbrNodeOnBoundaryV2+1:end,nbrNodeOnBoundaryH1 + nbrNodeOnBoundaryH2); %ok
    M8 = M(nbrNodeOnBoundaryV1+nbrNodeOnBoundaryV2+1:end,nbrNodeOnBoundaryH1); %ok
    
    MR = [[M3 M2 M1];[M6 M5 M4];[M8 M7 M0]];
elseif nbrBorderVertical == 0 && nbrBorderHorizontal == 2
    
    nonReducedSize = size(M,2);
    nbrNodeOnBoundary1 = size(subBoundariesHorizontal(1).node,1);
    nbrNodeOnBoundary2 = size(subBoundariesHorizontal(2).node,1);
    totalNbrNodeOnBoundaries = nbrNodeOnBoundary1+nbrNodeOnBoundary2;
    
    M0 = M(:,nbrNodeOnBoundary1+nbrNodeOnBoundary2+1:end);
    M1 = M(:,nbrNodeOnBoundary1+nbrNodeOnBoundary2);
    M2 = M(:,nbrNodeOnBoundary1);
    MR = [M2 M1 M0];     
else
    disp('error : no function to reduce matrix with 1,3,... or inequalsub-boundaries');
    MR = 0;
end
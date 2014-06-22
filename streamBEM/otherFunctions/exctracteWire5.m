function [wire2] = exctracteWire5(node,triangle,I,nbrContour,distanceBetween2Wire)

%%
% nbrnodeSide = 60;
% nbrnodeLength = 41;
% distanceBetween2Wire = 0.01;
% I = MatrixResults(3).I;
% nbrContour = 6;16;
%clear('wire','phiLevel');
% Lines Amplitude for the contour fonction
nContour = nbrContour; % Number of contour for the coil
phiMin = min(I);
phiMax = max(I);
phiContour = zeros(nContour,1);
for i=1:nContour
    phiContour(i) = phiMin + ((2*i-1)*(phiMax-phiMin))/(4*(nContour/2));
end

%figure
%clear 'phiLevel'
for i = 1:size(phiContour,1)
    phiCn=phiContour(i);
    k=1;
    if phiCn>=0
        phiLevel(i).currentDirection = 1;
    else
        phiLevel(i).currentDirection = -1;
    end
    % On each triangle we are checking if there is a level
    for j= 1 : size(triangle,1)
        % Let A,B and C be the point on the unit Triangle
        node1 = triangle(j,1);
        node2 = triangle(j,2);
        node3 = triangle(j,3);
        
        A = [0;0;I(node1)];
        B = [1;0;I(node2)];
        C = [0;1;I(node3)];
        % if t>=0 and t<1, then the plan with equation z = phiCn
        % Intersection the vector AB on point P:
        % P = A+t*(B-A)
        % Does the plan cross AB,BC or CA ?
        tBA = (phiCn - A(3))/(B(3)-A(3));
        tCB = (phiCn - B(3))/(C(3)-B(3));
        tAC = (phiCn - C(3))/(A(3)-C(3));
        nbIntersection = 0;
        if (tBA<=1 && tBA>= 0) 
            nbIntersection = nbIntersection +1;
            P(:,nbIntersection) = A+tBA*(B-A);
        end
        if (tCB<=1 && tCB>= 0)
            nbIntersection = nbIntersection +1;
            P(:,nbIntersection) = B+tCB*(C-B);
        end
        if  (tAC<=1 && tAC>= 0)
            nbIntersection = nbIntersection +1;
            P(:,nbIntersection) = C+tAC*(A-C);
        end
        if nbIntersection >= 3
            disp('Error : you can not have 3 intersectio between a triangle and a plane!')
        end
            
        % P are point on the unit triangle
        % p1 = [0 (phi3-phiCn)/(phi3-phi2) phiCn]+0.0*vectS;
        
        
        if nbIntersection >= 1
            % We can now interpolate between the two point P
            % for x comprised between x0 and x1
            % y = y0 + (y1-y0)*(x-x0)/(x1-x0);
            x = linspace(P(1,1),P(1,2),4);
            %We reset the number of point to zero, as the vector x now
            %include also the two previously found points.
            nbIntersection = 0;
            for l=1:size(x,2)
                nbIntersection = nbIntersection +1;
                Pinterpo(1,nbIntersection) = x(l);
                Pinterpo(2,nbIntersection) = P(2,1) + (P(2,2)-P(2,1))*(x(l)-P(1,1))/(P(1,2)-P(1,1));
                Pinterpo(3,nbIntersection) = P(3,1);
            end
            % We then remove the last and first point, as they are
            % present on all the triangle
            Pinterpo(:,4) = [];
            Pinterpo(:,1) = [];
            nbIntersection = nbIntersection-2;

            % LetD,E and F be the point on the real triangle 
            
            D = node(triangle(j,1),:); % A wich correspond to x3 for Poole
            E = node(triangle(j,2),:); % B which correspond to x1 for Poole
            F = node(triangle(j,3),:); % C which correspond to x2 for Poole
            P = Pinterpo;
            for l=1:nbIntersection
                w = 1-P(1,l)-P(2,l);
                A1 = [E(1)*P(1,l)+F(1)*P(2,l)+D(1)*w...
                      E(2)*P(1,l)+F(2)*P(2,l)+D(2)*w...
                      E(3)*P(1,l)+F(3)*P(2,l)+D(3)*w];
                  
                phiLevel(i).Coord(1,k) = A1(1);
                phiLevel(i).Coord(2,k) = A1(2);
                phiLevel(i).Coord(3,k) = A1(3);
                k=k+1;
            end

        end
    end
end
%displayWireAsPoint(phiLevel);
%%
% Sort the point into wire
%
%clear 'wire'
% phiLevel = phiLevelHand;
ii=0;
for i=1:size(phiLevel,2)
    ii = ii+1;% we increase the number of wire
    a = phiLevel(i).Coord(:,1); % We start with the first point, arbitrary
    wire(ii).Coord(:,1) = a; % We save the first point in the new wire Matrix
    wire(ii).currentDirection = phiLevel(i).currentDirection; % We also save the current direction, which is given by the current density amplitude
    phiLevel(i).Coord(:,1) = []; %We then remove it
    remainingPoint = size(phiLevel(i).Coord,2);
    nbrPointWire = 1;
    
    for k=1:remainingPoint % for all the point
        length1 = zeros(1,remainingPoint); % initialisation
        length2 = zeros(1,remainingPoint); % initialisation
        for j=1:remainingPoint
            % we calculate the norm of each vector
            c = phiLevel(i).Coord(1:3,j);
            length1(j) = norm(a-c); %norm
            length2(j) = dot(a,c); %norm
        end
%        [~,IX] = sort(length,'ascend'); % sort with the smaller distance first

        %[c1,index1] = min(length1);
        [~,index] = min(length1);
         % if the next point is close enough but not too close
        if length1(index)<=distanceBetween2Wire && length1(index) >0.0001
            %we make it the next point
            nbrPointWire = nbrPointWire+1;
            % and we change the reference point
            a = phiLevel(i).Coord(:,index);
            wire(ii).Coord(:,nbrPointWire) = a;
            %Delete it from the list
            phiLevel(i).Coord(:,index) = [];
        elseif length1(index)<=distanceBetween2Wire && length1(index) <0.0001
            % we just remove the point
            phiLevel(i).Coord(:,index) = [];            
        else
            % if the next point is too far, then it is another wire
            % and we change the reference point
            a = phiLevel(i).Coord(:,index);
            ii=ii+1;% Make a new wire
            wire(ii).Coord(:,1) = a; % We save the first point in the new wire Matrix
            wire(ii).currentDirection = phiLevel(i).currentDirection; % We also save the current direction, which is given by the current density amplitude
            phiLevel(i).Coord(:,index) = []; %We then remove it
            %nbrPoint = size(phiLevel(i).Coord,2);
            %set the new number of point
            nbrPointWire = 1;
        end
        remainingPoint = remainingPoint-1;
    end
end
%displayWireAsPoint(wire);
% displayWire(wire);
%% So we have to select the two first point in order to choose the rotation
% direction
% wireSave = wire;
%wire = wireSave;
%clear 'wire2'
%close all
%displayWire(wire(8));
ii=1;
%uu = 8;
for i=1:size(wire,2)
    % We start with the first point, arbitrary
    A = wire(i).Coord(:,1);
    % We save the first point in the new wire Matrix
    wire2(ii).Coord(:,1) = wire(i).Coord(:,1);
    %We then remove it
    wire(i).Coord(:,1) = []; 
    % We also save the current direction, which is given 
    % by the current density amplitude
    wire2(ii).currentDirection = wire(i).currentDirection;
    
    % lets concidere the referentiel with origine O (0;0;0)
    % Let G be the gravity point of the actual wire
    G = [mean(wire(i).Coord(1,:));mean(wire(i).Coord(2,:));mean(wire(i).Coord(3,:))];
    % our first vector is the one perpendicular to z
    % let the point pG be the orthogonal projection of G on z
    %pG = [0;0;G(3)];
    pG = [0;0;0];
    % Vect Z is thus the vector pGG
    vectZ = (G-pG)/norm(G-pG);% dot(vectZ,vectZ);

    
    nbrPointLoop = size(wire(i).Coord,2); % We save the number of point
    % we have to sort, because we are going to delete element of
    % the matrix. We still have to remove the first point thus
    for k=1:nbrPointLoop-1 % for all the point, exepted the last one
        
        % the vectX is the one pointing toward A
        % we first have to find the origine of the new repere
        % to do this we first have to find the plan equation
        % A is a point on the plan
        % vectZ is the normal of the plane
        % we have vectZ.(r-A) = 0, with r = [x;y;z]
        % the plane has for equation
        % dot(vectZ,[x;y;z]) = dot(vectZ,A)
        % vectZ(1)*x + vectZ(2)*y + vectZ(3)*z = dot(vectZ,A);
        % the origine is the intersection of the vectZ with the plane.
        % Any point P belonging to the vectZ can be written as
        % P = pG + t*vectZ
        % any any point belonging to the plane is defined respect the equation
        % vectZ(1)*x + vectZ(2)*y + vectZ(3)*z = dot(vectZ,A);
        % that is 
        % vectZ(1)*(0+t*vectZ(1)) + vectZ(2)*(0+t*vectZ(2)) +
        % vectZ(3)*(G(3)+t*vectZ(3)) = dot(vectZ,A);
        % t*(vectZ(1)^2+vectZ(2)^2+vectZ(3)^2) = dot(vectZ,A) -  vectZ(3)*G(3)
        t = (dot(vectZ,A) -  vectZ(3)*G(3))/(vectZ(1)^2+vectZ(2)^2+vectZ(3)^2);
        P0 = pG + t*vectZ;
        % now the vectX can be defined as P0A
        vectX = (A-P0)/norm(A-P0);% dot(vectX,vectX);
        % now, to define vectY, we use the vectoriel product
        % we now that, in a right hand coordinate system
        vectY = cross(-vectX,vectZ);% dot(vectY,vectY);
    
        %We calculate the distance between all the point remaining in this
        %contour with the actual point (A)
        
        %Find the approximate line of the 12 point
        % we use a maximum of 12 points, as the point where already sorted
        % and as their is 2 point per triangle, it means that the next
        % point should be in the next 6 triangle. This seems realistic.
        %nbrPoint = min(12,size(wire(i).Coord,2));
        nbrPoint = size(wire(i).Coord,2);
        
        length = zeros(1,nbrPoint); % initialisation
        angle = zeros(1,nbrPoint); % initialisation
        for j=1:nbrPoint
            B = wire(i).Coord(:,j);
            %we projet our point in our new coordinate system
            pBx = dot(B-P0,vectX);
            pBy = dot(B-P0,vectY);
            %we calculate the norm of AB in the new coordinate sys.
            length(j) = norm(B-A);
            % then we calculate the in plane angle
            angle(j) = atan2(pBy,pBx);%atan2(-pBy,pBx);
        end
        [~,IX] = sort(length,'ascend');
        % sort with the smaller distance first
        % Then, we go along the sorted area, as long as we do not find a
        % point which go in the clockwise direction (angle positif)
        % we continue
        kk=1;
        while angle(IX(kk))< 0 && kk<nbrPoint
            kk=kk+1;
        end
        
        if kk==nbrPoint
            fprintf('Error during the ordering of the point belonging to the wire \nno good match founded, i=%i , k=%i \n',i,k)
            break
        else
            % Then save this point to the matrix
            B = wire(i).Coord(:,IX(kk));
            wire2(ii).Coord(:,k+1) = B;
            % And use it as the new reference point
            % but it will just change the length calculation, not the
            % referentiel ! THIS DOESN'T WORK
            % Aold = A;
            A = B;
            wire(i).Coord(:,IX(kk)) = []; % Delete the used point
        end

    end
    ii=ii+1;
end
%displayWire(wire2);
%%
% displayWire(wire2);
%displayWireSubFigure(wire2)
%displayWire(wire);
%displayWireAsPoint(wire);
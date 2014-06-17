function [basis] = basisFunction4(node, triangle, center)
% Calculation of the basis function vector for each triangle
%
% node : a matrix with the 3d position of each node (in meter)
% triangle : a matrix linking 3 node together to form a triangle
% center : a vector with the coordinate of the center of the mesh




%% Value useful for the Gauss-Legendre integration
[u,v,ck] = triGaussPoints(2);
for i=1:size(u,1)
    w(i) = 1-u(i)-v(i);
end

%%
tic;
for m=1:size(node,2); %For every node
    for i=1:node(m).nbrTriangle; %Elements linked to the current node
        actualTriangle_i = node(m).linkToTriangle(i);

        % Set the 2 other nodes connected to the actual triangle,
        % and trying to keep the order. I think it is vain because
        % of the node reorganization made to sort the boundaries
        mNodeIndex = find(triangle(actualTriangle_i).node == m);
        aNodeIndex = mNodeIndex+1;
        bNodeIndex = mNodeIndex+2;
        if aNodeIndex>3; aNodeIndex = mod(aNodeIndex,3); end;
        if bNodeIndex>3; bNodeIndex = mod(bNodeIndex,3); end;

        Ami = node(triangle(actualTriangle_i).node(aNodeIndex)).coord; %A - coordinate of the second node of the i_th face of the n_th node
        Bmi = node(triangle(actualTriangle_i).node(bNodeIndex)).coord; %B - coordinate of the third node of the i_th face of the n_th node
        Cmi = node(m).coord; %C - coordinate of the first node of the i_th face of the n_th node

            
        
        norme_AC = norm(Cmi-Ami);
        norme_BC = norm(Cmi-Bmi);
        norme_AB = norm(Bmi-Ami);
        %norme_dmi = sin(acos((norme_BC^2+norme_AC^2-norme_AB^2)/(2*norme_BC*norme_AC)))*norme_AC; % Using Al Kashi's angle and then SohCahToa stuff
        % dmi is going from C to the middle of AB
        norme_dmi = sin(acos((-norme_BC^2+norme_AC^2+norme_AB^2)/(2*norme_AB*norme_AC)))*norme_AC;% Using Al Kashi's angle and then SohCahToa stuff

        
         air = (1/2)*norm(cross(Ami-Cmi,Bmi-Cmi));
%         if air<0 % use the signed air, to determine in which direction we rotate
%             % http://geomalgorithms.com/a01-_area.html
%             basis(m).triangle(i).A = Bmi;
%             basis(m).triangle(i).B = Ami;
%             basis(m).triangle(i).C = Cmi;
%             vmi = (Ami-Bmi)/(norme_AB*norme_dmi); %why do we divide by 2 times the air of the triangle? To match Poole definition. It comes from the normal. see peeren pdf 93
%         else
%             basis(m).triangle(i).A = Ami;
%             basis(m).triangle(i).B = Bmi;
%             basis(m).triangle(i).C = Cmi;
%             vmi = (Bmi-Ami)/(norme_AB*norme_dmi); %why do we divide by 2 times the air of the triangle? To match Poole definition
%         end
%         
        
        %we calculate the cross product of 2 vector of the
        %triangle, i.e. the normal of the triangle. This is a
        %vector
        scaling = 10^8;
        vectorNormal = cross(scaling*(Ami-Cmi),scaling*(Bmi-Cmi));
        % we want to compare this vector to the one going from the
        % center of the volume to on point of the surface.
        % we suppose that the center is at (0,0,0).
        % To avoid numeric issue, we scale the vector with a factor 10^8
        vectorOA = scaling*(Ami-center);
        % we now caclulate the dot product. If the dot product is
        % positive, the 2 vectors are going in the same direction.
        % This mean that the normal of the triangle ABC is going in
        % the right direction:
        % Otherwise we have to change that.
        dotProduct = dot(vectorNormal,vectorOA);
        % if the dot product is negativ,
        % it means that the normal vector is not pointing in the same
        % direction as OA
        % which happens if the vector CA is not at the right of the vector CB
        % thus then we have to change the orientation of AB
        % But this seems to be unprecise, probably due to the small norm of
        % the normal vector
        % we may also do that with the normed air http://geomalgorithms.com/a01-_area.html
        if dotProduct(1) < 0
            %vector BA
            %vmi = (Ami-Bmi)/(norme_AB*norme_dmi); %why do we divide by 2 times the air of the triangle? To match Poole definition
            vmi = (Ami-Bmi)/(2*air); %why do we divide by 2 times the air of the triangle? To match Poole definition
            vectorNormal = cross(Bmi-Cmi,Ami-Cmi)/norm(cross((Bmi-Cmi),(Ami-Cmi)));
            basis(m).triangle(i).A = Bmi;
            basis(m).triangle(i).B = Ami;
        else
            %vector AB
            %vmi = (Bmi-Ami)/(norme_AB*norme_dmi); %why do we divide by 2 times the air of the triangle? To match Poole definition
            vmi = (Bmi-Ami)/(2*air);
            vectorNormal = cross((Ami-Cmi),(Bmi-Cmi))/norm(cross((Ami-Cmi),(Bmi-Cmi)));
            basis(m).triangle(i).A = Ami;
            basis(m).triangle(i).B = Bmi;
        end

        basis(m).triangle(i).C = Cmi;
        basis(m).triangle(i).value = vmi;
        basis(m).triangle(i).vectorNormal = vectorNormal;
%         
        % calculate the rmi for the inductance
        % which is a referential change
        % It give us the position (x,y,z) of the integration point (u,v).
        basis(m).triangle(i).r_o = changtRef(basis(m).triangle(i).A, basis(m).triangle(i).B, basis(m).triangle(i).C,u, v, w);
    end
end

fprintf(' - Done in %5.0f sec.\n',toc);
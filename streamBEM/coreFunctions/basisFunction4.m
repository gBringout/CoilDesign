function [basis] = basisFunction4(node, triangle, center)


% Calculation of the basis function vector for each traingle
% in order to have for each point the corresponding v vector, as defined in
% the thesis of Poole page 63

% node = shield.node;
% triangle = shield.triangle;
% nodeLinkToTriangle = shield.nodeLinkToTriangle;
% m=375 %with 32 nodes
% i=19
% i=18

% m=547 %with 7 nodes
% i=5
% i=2
%
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
        %basis(m).triangle(i).scaling = sqrt(((Ami(2)-Cmi(2))*(Bmi(3)-Cmi(3))-(Ami(3)-Cmi(3))*(Bmi(2)-Cmi(2)))^2+...
        %                      ((Ami(3)-Cmi(3))*(Bmi(1)-Cmi(1))-(Ami(1)-Cmi(1))*(Bmi(3)-Cmi(3)))^2+...
        %                      ((Ami(1)-Cmi(1))*(Bmi(2)-Cmi(2))-(Ami(2)-Cmi(2))*(Bmi(1)-Cmi(1)))^2);
        basis(m).triangle(i).C = Cmi;
        basis(m).triangle(i).value = vmi;
        basis(m).triangle(i).vectorNormal = vectorNormal;
%         
        % calculate the rmi for the inductance
        % which is a referentiel changement
        % It give us the position (x,y,z) of the integration point (u,v).
        basis(m).triangle(i).r_o = changtRef(basis(m).triangle(i).A, basis(m).triangle(i).B, basis(m).triangle(i).C,u, v, w);
    end
end

fprintf(' - Done in %5.0f sec.\n',toc);
 


% %% Test to see the calculation of the basis function: succed
% % %tested with blender mesh 20x20_R100_H280.obj
% % 
% colorIndex = 1;
% %colorChoice = ['y' 'm' 'c' 'r' 'g' 'b' 'k'];
% colorChoice = ['y' 'y' 'y' 'y' 'y' 'y' 'y' 'r' 'y' 'y' 'y' 'y' 'y' 'y' 'y'];
% figure
% % Calculation of the basis function vector for each traingle
% for m=1:size(node,1); %For every node
%     colorIndex = mod(colorIndex+1,size(colorChoice,2))+1;
%     for i=1:size(nodeLinkToTriangle(m).triangle,1); %Elements linked to the current node
%         actualTriangle_i = nodeLinkToTriangle(m).triangle(i);
%         
%         
%         % Set the 2 other nodes connected to the actual triangle,
%         % and trying to keep the order. I think it is vain because
%         % of the node reorganization made to sort the boundaries
%         mNodeIndex = find(triangle(actualTriangle_i,:) == m);
%         aNodeIndex = mNodeIndex+1;
%         bNodeIndex = mNodeIndex+2;
%         if aNodeIndex>3; aNodeIndex = mod(aNodeIndex,3); end;
%         if bNodeIndex>3; bNodeIndex = mod(bNodeIndex,3); end;
% 
%         Ami = node(triangle(actualTriangle_i,aNodeIndex),:); %A - coordinate of the second node of the i_th face of the n_th node
%         Bmi = node(triangle(actualTriangle_i,bNodeIndex),:); %B - coordinate of the third node of the i_th face of the n_th node
%         Cmi = node(m,:); %C - coordinate of the first node of the i_th face of the n_th node
%         
%         %plot the triangle
%         hold all
%         %plot3(Ami(1),Ami(2),Ami(3),'*')
%         %plot3(Bmi(1),Bmi(2),Bmi(3),'*')
%         plot3(Cmi(1),Cmi(2),Cmi(3),'*','Color',colorChoice(colorIndex))
%         
%         norme_AC = norm(Cmi-Ami);
%         norme_BC = norm(Cmi-Bmi);
%         norme_AB = norm(Bmi-Ami);
%         %norme_dmi = sin(acos((norme_BC^2+norme_AC^2-norme_AB^2)/(2*norme_BC*norme_AC)))*norme_AC; % Using Al Kashi's angle and then SohCahToa stuff
%         % dmi is the going from C to the middle of AB
%         norme_dmi = sin(acos((-norme_BC^2+norme_AC^2+norme_AB^2)/(2*norme_AB*norme_AC)))*norme_AC;% Using Al Kashi's angle and then SohCahToa stuff
% 
%         %we calculate the cross product of 2 vector of the
%         %triangle, i.e. the normal of the triangle. This is a
%         %vector
%         vectorNormal = cross(Ami-Cmi,Bmi-Cmi);
%         % we want to compare this vector to the one going from the
%         % center of the volume to on point of the surface.
%         % we suppose that the center is at (0,0,0).
%         vectorOA = Ami;
%         % we now caclulate the dot product. If the dot product is
%         % positive, the 2 vectors are going in the same direction.
%         % This mean that the normal of the triangle ABC is going in
%         % the right direction:
%         % Otherwise we have to change that.
%         dotProduct = dot(vectorNormal,vectorOA);
%         % if the cross product is negativ,
%         % it means that the normal vector is not pointing in the right direction
%         % which happens if the vector CA is not at the right of the vector CB
%         % thus then we have to change the orientation of AB
%         if dotProduct > 0
%             %vector AB
%             vmi = (Bmi-Ami)/(norme_AB*norme_dmi);
%             vectorNormal = cross(Bmi-Cmi,Ami-Cmi);
%             if(Ami(2)>0)
%                 h = quiver3(Ami(1),Ami(2),Ami(3),vmi(1),vmi(2),vmi(3),0.00010,'Color',colorChoice(colorIndex));
%                 adjust_quiver_arrowhead_size(h, 10);
%             end
%         else
%             %vector BA
%             vmi = (Ami-Bmi)/(norme_AB*norme_dmi);
%             if(Bmi(2)>0)
%                 h = quiver3(Bmi(1),Bmi(2),Bmi(3),vmi(1),vmi(2),vmi(3),0.00010,'Color',colorChoice(colorIndex),'LineWidth',2);
%                 adjust_quiver_arrowhead_size(h, 10);
%             end
%         end
%         
%         basis(m,i).value = vmi;
%         basis(m,i).vectorNormal = vectorNormal;
%     end
% end
% xlabel('x axis');
% ylabel('y axis');
% zlabel('z axis');

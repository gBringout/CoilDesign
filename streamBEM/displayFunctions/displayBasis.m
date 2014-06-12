function [] = displayBasis(node,basis,center)

%% Test plot basisFunction : ok

%% Test to see the calculation of the basis function: succed
% %tested with blender mesh 20x20_R100_H280.obj
% 
colorIndex = 1;
colorChoice = ['y' 'm' 'c' 'r' 'g' 'b' 'k'];
%colorChoice = ['y' 'y' 'y' 'y' 'y' 'y' 'y' 'r' 'y' 'y' 'y' 'y' 'y' 'y' 'y'];
figure
% Calculation of the basis function vector for each traingle
for m=61:15:100%1:37:size(node,2); %For every node
    colorIndex = mod(colorIndex+1,size(colorChoice,2))+1;
    for i=1:node(m).nbrTriangle; %Elements linked to the current node

        colorIndex = mod(colorIndex+1,size(colorChoice,2))+1;
        Ami = basis(m).triangle(i).A; %A - coordinate of the second node of the i_th face of the n_th node
        Bmi = basis(m).triangle(i).B; %B - coordinate of the third node of the i_th face of the n_th node
        Cmi = basis(m).triangle(i).C; %C - coordinate of the first node of the i_th face of the n_th node
        vmi = basis(m).triangle(i).value;
        nmi = basis(m).triangle(i).vectorNormal;
        
        %plot the triangle
        hold all
        %plot3(Ami(1),Ami(2),Ami(3),'*')
        %plot3(Bmi(1),Bmi(2),Bmi(3),'*')
        plot3(Cmi(1),Cmi(2),Cmi(3),'*','Color',colorChoice(colorIndex))
        centerTriangle = (Ami+Bmi+Cmi)/3;
        text(centerTriangle(1),centerTriangle(2),centerTriangle(3),num2str(i),'HorizontalAlignment','left','FontSize',8);

        %if(Ami(2)>0) % To select just the positive one
            h = quiver3(Ami(1),Ami(2),Ami(3),vmi(1),vmi(2),vmi(3),0.0005,'Color',colorChoice(colorIndex));
            adjust_quiver_arrowhead_size(h, 5);
            h = quiver3(Ami(1),Ami(2),Ami(3),nmi(1),nmi(2),nmi(3),1,'Color',colorChoice(colorIndex));
            adjust_quiver_arrowhead_size(h, 5);
        %end

    end
end

%plot3(center(1),center(2),center(3),'*')
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');

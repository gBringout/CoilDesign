function [L] = Lmn10(node_1, triangle_1, basis_1,node_2, triangle_2, basis_2)

% This file aim at calculating the so called "Lmn matrix"
% This is a newer version as Lmn7, using the new structure of object node
% and triangle
%
% node : a matrix with the 3d position of each node (in meter)
% triangle : a matrix linking 3 node together to form a triangle
% basis : a matrix with the description of the basis function
%
%changelog : try to increase the integration degree
% node_1 = coil.node; triangle_1 = coil.triangle; basis_1 = coil.basis;
if nargin < 4
	node_2 = node_1;
	basis_2 = basis_1; 
     triangle_2 = triangle_1;
	typeCalculation = 'sameSurface';
else
	typeCalculation = 'differentSurfaces';
end

[~,~,ck] = triGaussPoints(2);
cl = ck;

%%
mu_0 = 4*pi*10^-7;
coef = mu_0/(4*pi);
nbrIntegrationPoints = size(ck,1);

L = zeros(size(node_1,2),size(node_2,2)); % 1 = m, 2 = n

[TF,~] = license('checkout', 'Distrib_Computing_Toolbox');
numWorkers = 0;
if TF
    schd = findResource('scheduler', 'configuration', 'local');
    numWorkers = schd.ClusterSize;
end
if matlabpool('size') == 0  && TF && numWorkers >1
    % checking to see if the pool is already open and of we have the licence
    matlabpool open
end

tic

parfor m=1:size(node_1,2); %For every node
    temp = zeros(1,size(node_2,2));
	if strcmp(typeCalculation,'sameSurface')
		startPointSecondLoop = m; %for the mutual inductance of the same surface
		% we can reduce the calculation by using the symetrie of the matrix
	else
		startPointSecondLoop = 1;
		%otherwise, we have to make the whole calculation
	end
    for n=startPointSecondLoop:size(node_2,2); %and every node
    superBigSum = 0;
        for i=1:node_1(m).nbrTriangle; %Elements linked to the current node
            currentTriangle_i = node_1(m).linkToTriangle(i);
            vmi = basis_1(m).triangle(i).value;
            % factor linked to the integration of 1/(rmi-rnj)
            r_o = basis_1(m).triangle(i).r_o;
            
            % We calculate the scaling factor due to the axis change during the integration
            %scallingMi = basis_1(m).triangle(i).scaling;%
            
            for j=1:node_2(n).nbrTriangle; %Number of face in this face
                currentTriangle_j = node_2(n).linkToTriangle(j);
                vnj = basis_2(n).triangle(j).value;  %#ok<PFBNS>
                integral = 0;

                if currentTriangle_i~=currentTriangle_j || strcmp(typeCalculation,'differentSurfaces') %ok
					% If the current triangles are not the same,
					% then we will not share any basis function. We can thus solve that easily

                    % We calculate the scaling factur due to the axis change
                    % during the integration
                    %scallingNj = basis_2(n).triangle(j).scaling;
                    % factor linked to the integration of 1/(rmi-rnj)
                    r_p = basis_2(n).triangle(j).r_o;
                    
                    % Then we calculate the integral (see my lab book)
                    %sum = integral(m).sliced(i,n,j);
%                     for k1 = 1:nbrIntegrationPoints
%                         for k2 = 1:nbrIntegrationPoints
%                             integral = integral + ck(k2)*ck(k1)*1/norm(rmi(k2,:)-rnj(k1,:));
%                         end
%                     end
                    for k = 1:nbrIntegrationPoints
                        for l = 1:nbrIntegrationPoints
                            integral = integral + cl(l)/norm(r_o(k,:)-r_p(l,:));
                        end
                        integral = ck(k)*integral;
                    end

                    % The division here is not logic for me
                    % but it match the results from Sanchez
                    bigSum = integral*(2*triangle_2(currentTriangle_j).air*2*triangle_1(currentTriangle_i).air);

                else
                    % If the 2 triangle are similar, we use an approximate
                    % calculation, according to page 72 equ 5.40 of Poole
                    % thesis
                    Ami = basis_1(m).triangle(i).A;
                    Bmi = basis_1(m).triangle(i).B;
                    Cmi = basis_1(m).triangle(i).C;
                    a = dot(Cmi-Ami,Cmi-Ami);
                    b = dot(Cmi-Ami,Cmi-Bmi);
                    c = dot(Cmi-Bmi,Cmi-Bmi);
                    
                    sa = sqrt(a);
                    sc = sqrt(c);
                    ss = sqrt(a-2*b+c);
                    sac = sqrt(a*c);
                    
                    % The division here is not logic for me
                    % but it match the results from Sanchez
                    % (1/(4*triangle_1(currentTriangle_i).air^2))
                    bigSum = (1*(4*triangle_1(currentTriangle_i).air^2))*(...
                             1/(6*sa) * log( ( (a-b+sa*ss)*(b+sac) )/( (-b+sac)*(-a+b+sa*ss) ))...
                            +1/(6*sc) * log( ( (b+sac)*(-b+c+sc*ss) )/( (b-c+sc*ss)*(-b+sac) ))...
                            +1/(6*ss) * log( ( (a-b+sa*ss)*(-b+c+sc*ss) )/( (b-c+sc*ss)*(-a+b+sa*ss) ))); %#ok<PFBNS>
                end

                vmi_vnj = dot(vmi,vnj);
                superBigSum = superBigSum+vmi_vnj*bigSum;
            end
        end
        %L(m,n) = coef*superBigSum;
        temp(n) = coef*superBigSum;
    end
    L(m,:) = temp;
    %fprintf('Loop m=%i \n',m);
end

% if we calculated only the upper half of the matrix, we have to copy the other half
if strcmp(typeCalculation,'sameSurface')
	for m=2:size(node_1,2) % m is a node
		for n=1:m-1 %n is a node
			L(m,n) = L(n,m);
		end
	end
end

fprintf(' - Done in %5.0f sec.\n',toc);
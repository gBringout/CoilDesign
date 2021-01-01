function [BX_a,BY_a,BZ_a] = Field3(loops,I,position_x,position_y,position_z)
% This file aim at calculating the B field accroding to Biot-savart equation 
% at some given point. This is the same calculation as Field2, just for
% points (and not a grid)
% position vectors have to have the same size

% We have to provide
% loops : [ ] The description of the loops, as segment with a given
% direction
% I : [A] current flowing in the wire
% position_x : [m] x Index component of the point

% output
% BX_a : magnetic flux density in the xIndex direction
%
% created by Gael Bringout - xx.xx.2010

%% calculate the size of the elements
nbrPointX = size(position_x,2);
nbrLoops = size(loops,2);

%% Pre-allocate
BX_a = zeros(nbrPointX,1);
BY_a = zeros(nbrPointX,1);
BZ_a = zeros(nbrPointX,1);

%% Constante
mu0 = 4*pi*10^-7;
mip4 = (mu0*I)/(4*pi);
    
%% Field calculation

disp('Calculation of the Field');

%activate the parallel function if available
if license('test','Distrib_Computing_Toolbox')
    matlabVersion = version;
    matlabVersion = str2num(matlabVersion(1:3));
    if matlabVersion < 8.2
        [TF,~] = license('checkout', 'Distrib_Computing_Toolbox');
        if TF
            schd = findResource('scheduler', 'configuration', 'local');
            numWorkers = schd.ClusterSize;
        end

        if matlabpool('size') == 0  && TF && numWorkers >1
            % checking to see if the pool is already open and of we have the licence
            % and at least 2 cores
            matlabpool open
        end
    elseif matlabVersion >= 8.2 
        poolobj = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(poolobj)
            parpool;
        end
    end
end

for xIndex = 1:nbrPointX;
    xValue = position_x(xIndex);
    yValue = position_y(xIndex);
    zValue = position_z(xIndex);
    sum_X = 0;
    sum_Y = 0;
    sum_Z = 0;

    % for each loops
    for w=1:nbrLoops
        direction = loops(w).currentDirection;
        coord_x = loops(w).Coord(1,:);
        coord_y = loops(w).Coord(2,:);
        coord_z = loops(w).Coord(3,:);

        previous_lz = coord_z(1);
        previous_lx = coord_x(1);
        previous_ly = coord_y(1);
        % for each segment for the actual loops
        for u=2:size(loops(w).Coord,2)

            current_lx = coord_x(u);
            current_ly = coord_y(u);
            current_lz = coord_z(u);

            % x,y and z component of the length of the actual wire element
            lx = (current_lx-previous_lx);
            ly = (current_ly-previous_ly);
            lz = (current_lz - previous_lz);

            % x,y and z component of distance between the wire
            % element middle point and the point at which we want
            % to calculate the field
            rx = xValue - ((current_lx + previous_lx)/2);
            ry = yValue - ((current_ly + previous_ly)/2);
            rz = zValue - ((current_lz + previous_lz)/2);

            % norm at the power of three of the distance
            norm_vector = (sqrt(rx^2+ry^2+rz^2))^3;
            coef_vector = direction*mip4 / norm_vector;

            % Biot-savart cross product and integral
            sum_X = sum_X + coef_vector * (ly*rz-lz*ry);
            sum_Y = sum_Y + coef_vector * (lz*rx-lx*rz);
            sum_Z = sum_Z + coef_vector * (lx*ry-ly*rx);

            previous_lz = current_lz;
            previous_lx = current_lx;
            previous_ly = current_ly;
        end
    end
    BX_a(xIndex) = sum_X;
    BY_a(xIndex) = sum_Y;
    BZ_a(xIndex) = sum_Z;
end

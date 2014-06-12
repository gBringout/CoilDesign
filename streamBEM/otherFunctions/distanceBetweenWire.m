function [minimalWireSpacing] = distanceBetweenWire(Wire)

%%
% Calculate the smaller distance between two wire

minimalWireSpacing.value = 1;
%tic
for i = 1:size(Wire,2) % for each wire
    otherwire = 1:size(Wire,2);
    otherwire(i) = [];% Remove the actual element
    for j=1:size(Wire(i).Coord,2) %For each point in this wire
        Ax = Wire(i).Coord(1,j);
        Ay = Wire(i).Coord(2,j);
        Az = Wire(i).Coord(3,j);
        for k=1:size(otherwire,2) %and for each other wire
            kk = otherwire(k);
            for l=1:size(Wire(kk).Coord,2) %and for each point of this wire
               %Calculate the distance between the two point 
                Bx = Wire(kk).Coord(1,l);
                By = Wire(kk).Coord(2,l);
                Bz = Wire(kk).Coord(3,l);
                distance = sqrt((Ax-Bx)^2+(Ay-By)^2+(Az-Bz)^2);
                if distance<minimalWireSpacing.value
                    minimalWireSpacing.value = distance;
                    minimalWireSpacing.Wire1 = i;
                    minimalWireSpacing.Point1 = j;
                    minimalWireSpacing.Wire2 = kk;
                    minimalWireSpacing.Point2 = l;
                end
            end
        end
    end
end
%tElapsed=toc;
%fprintf('DistanceBetweenWire done in %5.0f sec. \n',toc);
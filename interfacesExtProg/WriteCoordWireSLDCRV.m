function [] = WriteCoordWireSLDCRV(name,wire)

% We supose that Solidworks is set to use mm (default setup)
wire = cleanWire(wire);%remove all duplicate wire.
for i=1:size(wire,2)
    filename = sprintf('%s_Wire_%i.sldcrv',name,i);
    fileID = fopen(filename, 'w');
% with the functin cleanWire, we don't have to clean the wire by hand.
%     %We have first to remove a point which is always two time at the near
%     %end of the data set
%     sizeArray = size(wire(i).Coord,2);
%     wire(i).Coord(:,sizeArray-1) = [];
    fprintf(fileID,'%1.9f \t %1.9f \t %1.9f\n',wire(i).Coord*1000);
    fclose(fileID);
end

function [] = WriteCoordWireFastHenry(wire,wireName,wireWidth,wireHeigth,fmin,fmax,nPointPerDecade)
% 
% wireWidth = 0.0075;
% wireHeigth = 0.0075;
% fmin = 1;
% fmax = 10^6;
% nPointPerDecade = 5;

filename = sprintf('%s_input.fastHenry.inp',wireName);
    
if size(wire,2) ~= 1
    %disp('Error, you have to give a single wire path!')
    %i=1;
    fileID = fopen(filename, 'w');
    %We have first to remove a point which is always two time at the near
    %end of the data set
    fprintf(fileID,'**This is the title line. It will always be ignored**.\n');
    fprintf(fileID,'* Everything is case INsensitive\n');
    fprintf(fileID,'* An asterisk starts a comment line.\n');
    fprintf(fileID,'*\n');
    fprintf(fileID,'* The following line names millimeters as the length units for the rest\n');
    fprintf(fileID,'* of the file.\n');
    fprintf(fileID,'.Units MM\n');
    fprintf(fileID,'* Make z=0 the default z coordinate and copper the default conductivity.\n');
    fprintf(fileID,'* Note that the conductivity is in units 1/(mm*Ohms), not 1/(m*Ohms)\n');
    fprintf(fileID,'* since the default units are millimeters.\n');
    fprintf(fileID,'.Default z=0 sigma=2.6029e4\n');
    fprintf(fileID,'* (Litz wire 10000x0.063)\n');
    fprintf(fileID,'\n');
    fprintf(fileID,'* The nodes of a wire\n');
    
    node = 0;
    for i=1:size(wire,2)
        if 1;%wire(i).currentDirection == 1 % if the current go in one direction
            for j=1:size(wire(i).Coord,2);
                node = node +1;
                fprintf(fileID,'N%i x=%1.3f y=%1.3f z=%1.3f \n',node,wire(i).Coord(1,j)*1000,wire(i).Coord(2,j)*1000,wire(i).Coord(3,j)*1000);
            end
        else %other go in the other direction. Is this completly wrong ? Node don't have direction
            for j=size(wire(i).Coord,2):-1:1
                node = node +1;
                fprintf(fileID,'N%i x=%1.3f y=%1.3f z=%1.3f \n',node,wire(i).Coord(1,j)*1000,wire(i).Coord(2,j)*1000,wire(i).Coord(3,j)*1000);
            end
        end
    end

    fprintf(fileID,'* The segments connecting the nodes\n');
    node = 0;
    segment = 0;
    wirePort = 0;
    for i=1:size(wire,2)
		if 1;%wire(i).currentDirection == 1 % if the current go in one direction
			wirePort(i,1) = node+1;%We save the segment number of the fisrt wire piece
			for j=1:size(wire(i).Coord,2)-1 % we remove the last point because it is not connected to anything.
				node = node +1;
				segment = segment +1;
				fprintf(fileID,'E%i N%i N%i w=%1.3f h=%1.3f \n',segment,node,node+1,wireWidth*1000,wireHeigth*1000);
			end
			node = node +1; % compensate for the last "j" point we remove in the inner for loop
			wirePort(i,2) = node;%We save the segment number of the last wire piece
		else
			%we first start at the end of the wire? NO ! We just change the order of the point
			wirePort(i,1) = node+1;%We save the segment number of the fisrt wire piece
			for j=1:size(wire(i).Coord,2)-1 % we remove the last point because it is not connected to anything.
				node = node +1;
				segment = segment +1;
				fprintf(fileID,'E%i N%i N%i w=%1.3f h=%1.3f \n',segment,node+1,node,wireWidth*1000,wireHeigth*1000); % We just change the sens of the segment
			end
			node = node +1; % compensate for the last "j" point we remove in the inner for loop
			wirePort(i,2) = node;%We save the segment number of the last wire piece
		end
    end
    
    
    fprintf(fileID,'* define one ''port'' of the network\n');
    for i=1:size(wirePort,1)
		if wire(i).currentDirection == 1 % if the current go in one direction
			fprintf(fileID,'.external N%i N%i \n',wirePort(i,1),wirePort(i,2));
		else
			fprintf(fileID,'.external N%i N%i \n',wirePort(i,2),wirePort(i,1));
		end
    end

    fprintf(fileID,'* Frequency range of interest.\n');
    fprintf(fileID,'.freq fmin=%1.1e fmax=%1.1e ndec=%i \n',fmin,fmax,nPointPerDecade);
    
    fprintf(fileID,'* All input files must end with:\n');
    fprintf(fileID,'.End\n');
    fclose(fileID);
else
    i=1;
    fileID = fopen(filename, 'w');
    %We have first to remove a point which is always two time at the near
    %end of the data set
    fprintf(fileID,'**This is the title line. It will always be ignored**.\n');
    fprintf(fileID,'* Everything is case INsensitive\n');
    fprintf(fileID,'* An asterisk starts a comment line.\n');
    fprintf(fileID,'*\n');
    fprintf(fileID,'* The following line names millimeters as the length units for the rest\n');
    fprintf(fileID,'* of the file.\n');
    fprintf(fileID,'.Units MM\n');
    fprintf(fileID,'* Make z=0 the default z coordinate and copper the default conductivity.\n');
    fprintf(fileID,'* Note that the conductivity is in units 1/(mm*Ohms), not 1/(m*Ohms)\n');
    fprintf(fileID,'* since the default units are millimeters.\n');
    fprintf(fileID,'.Default z=0 sigma=2.6029e4\n');
    fprintf(fileID,'* (Litz wire 10000x0.063)\n');
    fprintf(fileID,' \n');
    fprintf(fileID,'* The nodes of a wire\n');
    for i=1:size(wire(1).Coord,2);
        fprintf(fileID,'N%i x=%1.3f y=%1.3f z=%1.3f \n',i,wire(1).Coord(1,i)*1000,wire(1).Coord(2,i)*1000,wire(1).Coord(3,i)*1000);
    end

    fprintf(fileID,'* The segments connecting the nodes\n');
    for i=1:size(wire(1).Coord,2)-1
        fprintf(fileID,'E%i N%i N%i w=%1.3f h=%1.3f \n',i,i,i+1,wireWidth*1000,wireHeigth*1000);
    end
    
    
    fprintf(fileID,'* define one ''port'' of the network\n');
    fprintf(fileID,'.external N1 N%i \n',i+1);

    fprintf(fileID,'* Frequency range of interest.\n');
    fprintf(fileID,'.freq fmin=%1.1e fmax=%1.1e ndec=%i \n',fmin,fmax,nPointPerDecade);
    
    fprintf(fileID,'* All input files must end with:\n');
    fprintf(fileID,'.End\n');
    fclose(fileID);
end

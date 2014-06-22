function [node,triangle,tri] = importMeshWavefront(filename)
% Read an single OBJ object exported by Blender or 3DSmax for exemple
% This has to be a mesh made with triangular elements

fprintf('Using file %s \n',filename)
fid = fopen(filename); %Open the file
tline = fgets(fid); % get the first line of the files
nbrNode = 0; %initialize the number of node
nbrTriangle = 0; %initialize the number of triangle

while ischar(tline) %when we reach the end of the files tline should send back -1
    if tline(1) == 'v' %If we have a "Vertex"
        nbrNode = nbrNode + 1;
        x = str2num(tline(2:end)); %take the rest of the string to extract the coordinates
        node(nbrNode,1) = x(1);
        node(nbrNode,2) = x(2);
        node(nbrNode,3) = x(3);
    elseif tline(1) == 'f' %if we have a "Face"
        nbrTriangle = nbrTriangle + 1;
        x = str2num(tline(2:end)); %take the rest of the string to extract the relation
        triangle(nbrTriangle,1) = x(1);
        triangle(nbrTriangle,2) = x(2);
        triangle(nbrTriangle,3) = x(3);
    end
    tline = fgets(fid); % get the next line
end

fclose(fid);
%disp('Reading ok')
tri = TriRep(triangle,node(:,1),node(:,2),node(:,3));
%figure
%trimesh(tri);
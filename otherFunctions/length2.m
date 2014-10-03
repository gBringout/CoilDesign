function [length] = length2(Wire)

% Approximation of the length of all the loops
%Wires = Wire;
length = 0;
for i = 1:size(Wire,2)
	for j=1:size(Wire(i).Coord,2)-1 %Change the size(x,1) for size(x,2)
		A = Wire(i).Coord(:,j);
		B = Wire(i).Coord(:,j+1);
		norm_AB = sqrt((B(1,1)-A(1,1))^2+(B(2,1)-A(2,1))^2+(B(3,1)-A(3,1))^2);
		length = length + norm_AB;
	end
end
		
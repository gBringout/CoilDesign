function [maxHomoX,maxHomoY,maxHomoZ,meanHomoX,meanHomoY,meanHomoZ,Bx0,By0,Bz0] = Homogeneity(Bx,By,Bz,sphere_radius,X_Value,Y_Value,Z_Value)
% Calculate the maximumum in-homogeneity of a field in a sphere

x0 = (size(X_Value,2)-1)/2 +1;
y0 = (size(Y_Value,2)-1)/2 +1;
z0 = (size(Z_Value,2)-1)/2 +1;
Bx0 = Bx(x0,y0,z0);
By0 = By(x0,y0,z0);
Bz0 = Bz(x0,y0,z0);

homoX = zeros(size(X_Value,2),size(Y_Value,2),size(Z_Value,2));
homoY = zeros(size(X_Value,2),size(Y_Value,2),size(Z_Value,2));
homoZ = zeros(size(X_Value,2),size(Y_Value,2),size(Z_Value,2));

l = 0;
meanHomoX = 0;
meanHomoY = 0;
meanHomoZ = 0;

for i=1:size(X_Value,2)
    for j=1:size(Y_Value,2)
        for k=1:size(Z_Value,2)
            if sqrt(X_Value(i)^2+Y_Value(j)^2+Z_Value(k)^2)<=sphere_radius % in the spehre
                homoX(i,j,k) = abs((Bx(i,j,k)/Bx0))-1;
                homoY(i,j,k) = abs((By(i,j,k)/By0))-1;
                homoZ(i,j,k) = abs((Bz(i,j,k)/Bz0))-1;
                meanHomoX = meanHomoX + homoX(i,j,k);
                meanHomoY = meanHomoY + homoY(i,j,k);
                meanHomoZ = meanHomoZ + homoZ(i,j,k);
                l = l+1;
            end
        end
    end
end

maxHomoX = max(max(max(homoX)));
maxHomoY = max(max(max(homoY)));
maxHomoZ = max(max(max(homoZ)));
function [GradientLinearity,Gradient] = Linearity(Bx,By,Bz,sphere_radius,X_Value,Y_Value,Z_Value)
% Calculate the maximumum non linearity of a field in a sphere

x0 = (size(X_Value,2)-1)/2 +1;
y0 = (size(Y_Value,2)-1)/2 +1;
z0 = (size(Z_Value,2)-1)/2 +1;

LinearityGBx_x = zeros(size(X_Value,2),size(Y_Value,2),size(Z_Value,2));
LinearityGBx_y = zeros(size(X_Value,2),size(Y_Value,2),size(Z_Value,2));
LinearityGBx_z = zeros(size(X_Value,2),size(Y_Value,2),size(Z_Value,2));
LinearityGBy_x = zeros(size(X_Value,2),size(Y_Value,2),size(Z_Value,2));
LinearityGBy_y = zeros(size(X_Value,2),size(Y_Value,2),size(Z_Value,2));
LinearityGBy_z = zeros(size(X_Value,2),size(Y_Value,2),size(Z_Value,2));
LinearityGBz_x = zeros(size(X_Value,2),size(Y_Value,2),size(Z_Value,2));
LinearityGBz_y = zeros(size(X_Value,2),size(Y_Value,2),size(Z_Value,2));
LinearityGBz_z = zeros(size(X_Value,2),size(Y_Value,2),size(Z_Value,2));

l = 0;

X_Step = abs(X_Value(1)-X_Value(2));
Y_Step = abs(Y_Value(1)-Y_Value(2));
if size(Z_Value,2) > 1
    Z_Step = abs(Z_Value(1)-Z_Value(2));
else
    Z_Step = 1;
end
if size(Z_Value,2) > 1
    [GBx_y GBx_x GBx_z] = gradient(Bx,X_Step,Y_Step,Z_Step);
    [GBy_y GBy_x GBy_z] = gradient(By,X_Step,Y_Step,Z_Step);
    [GBz_y GBz_x GBz_z] = gradient(Bz,X_Step,Y_Step,Z_Step);
else
    [GBx_y GBx_x] = gradient(Bx,X_Step,Y_Step);
        GBx_z = zeros(size(X_Value,2),size(Y_Value,2));
    [GBy_y GBy_x] = gradient(By,X_Step,Y_Step);
        GBy_z = zeros(size(X_Value,2),size(Y_Value,2));
    [GBz_y GBz_x] = gradient(Bz,X_Step,Y_Step);
        GBz_z = zeros(size(X_Value,2),size(Y_Value,2));
end

GBx_x0 = GBx_x(x0,y0,z0);
GBx_y0 = GBx_y(x0,y0,z0);
GBx_z0 = GBx_z(x0,y0,z0);
GBy_x0 = GBy_x(x0,y0,z0);
GBy_y0 = GBy_y(x0,y0,z0);
GBy_z0 = GBy_z(x0,y0,z0);
GBz_x0 = GBz_x(x0,y0,z0);
GBz_y0 = GBz_y(x0,y0,z0);
GBz_z0 = GBz_z(x0,y0,z0);


for i=1:size(X_Value,2)
    for j=1:size(Y_Value,2)
        for k=1:size(Z_Value,2)
            if sqrt(X_Value(i)^2+Y_Value(j)^2+Z_Value(k)^2)<=sphere_radius % in the spehre
                LinearityGBx_x(i,j,k) = abs((GBx_x(i,j,k)/GBx_x0))-1;
                LinearityGBx_y(i,j,k) = abs((GBx_y(i,j,k)/GBx_y0))-1;
                LinearityGBx_z(i,j,k) = abs((GBx_z(i,j,k)/GBx_z0))-1;
                
                LinearityGBy_x(i,j,k) = abs((GBy_x(i,j,k)/GBy_x0))-1;
                LinearityGBy_y(i,j,k) = abs((GBy_y(i,j,k)/GBy_y0))-1;
                LinearityGBy_z(i,j,k) = abs((GBy_z(i,j,k)/GBy_z0))-1;
                
                LinearityGBz_x(i,j,k) = abs((GBz_x(i,j,k)/GBz_x0))-1;
                LinearityGBz_y(i,j,k) = abs((GBz_y(i,j,k)/GBz_y0))-1;
                LinearityGBz_z(i,j,k) = abs((GBz_z(i,j,k)/GBz_z0))-1;
                l = l+1;
            end
        end
    end
end

maxLinearityGBx_x = max(max(max(LinearityGBx_x)));
maxLinearityGBx_y = max(max(max(LinearityGBx_y)));
maxLinearityGBx_z = max(max(max(LinearityGBx_z)));

maxLinearityGBy_x = max(max(max(LinearityGBy_x)));
maxLinearityGBy_y = max(max(max(LinearityGBy_y)));
maxLinearityGBy_z = max(max(max(LinearityGBy_z)));

maxLinearityGBz_x = max(max(max(LinearityGBz_x)));
maxLinearityGBz_y = max(max(max(LinearityGBz_y)));
maxLinearityGBz_z = max(max(max(LinearityGBz_z)));

GradientLinearity(1,1) = maxLinearityGBx_x;
GradientLinearity(1,2) = maxLinearityGBx_y;
GradientLinearity(1,3) = maxLinearityGBx_z;
GradientLinearity(2,1) = maxLinearityGBy_x;
GradientLinearity(2,2) = maxLinearityGBy_y;
GradientLinearity(2,3) = maxLinearityGBy_z;
GradientLinearity(3,1) = maxLinearityGBz_x;
GradientLinearity(3,2) = maxLinearityGBz_y;
GradientLinearity(3,3) = maxLinearityGBz_z;

Gradient(1,1) = GBx_x0;
Gradient(1,2) = GBx_y0;
Gradient(1,3) = GBx_z0;
Gradient(2,1) = GBy_x0;
Gradient(2,2) = GBy_y0;
Gradient(2,3) = GBy_z0;
Gradient(3,1) = GBz_x0;
Gradient(3,2) = GBz_y0;
Gradient(3,3) = GBz_z0;
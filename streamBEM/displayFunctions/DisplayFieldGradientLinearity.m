function [] = DisplayFieldGradientLinearity(B,wire,coil_radius,coil_length,procentPlane,procentVolume,x_Value,y_Value,z_Value)
% display the non linearity of a field
% B is a 3d field (So, it has to be Bx, By, Bz or Babs)

x_Start =x_Value(1);
x_Stop = x_Value(size(x_Value,2));
x_Step = abs(x_Value(1)-x_Value(2));
% Y is a position, in meter
y_Start =y_Value(1);
y_Stop = y_Value(size(y_Value,2));
y_Step = abs(y_Value(1)-y_Value(2));
% Z is a position, in meter
z_Start =z_Value(1);
z_Stop = z_Value(size(z_Value,2));
z_Step = abs(z_Value(1)-z_Value(2));

x0 = (size(x_Value,2)-1)/2+1;
y0 = (size(y_Value,2)-1)/2+1;
z0 = (size(z_Value,2)-1)/2+1;

% To draw the coil
xp_Start =-coil_radius;
xp_Stop = coil_radius;
xp_Step = coil_radius/20;
% Y is a position, in meter
yp_Start = -coil_radius;
yp_Stop = coil_radius;
yp_Step = coil_radius/20;
% Z is a position, in meter
zp_Start = -coil_length/2;
zp_Stop = coil_length/2;
zp_Step = coil_length/20;
% X is a position, in meter
xp_Value = xp_Start:xp_Step:xp_Stop;
% Y is a position, in meter
yp_Value = yp_Start:yp_Step:yp_Stop;
% Z is a position, in meter
zp_Value = zp_Start:zp_Step:zp_Stop;

%% Calculation of the linearity in space

LinearityGB_x = zeros(size(x_Value,2),size(y_Value,2),size(z_Value,2));
LinearityGB_y = zeros(size(x_Value,2),size(y_Value,2),size(z_Value,2));
LinearityGB_z = zeros(size(x_Value,2),size(y_Value,2),size(z_Value,2));
l = 0;

[GB_y GB_x GB_z] = gradient(B,x_Step,y_Step,z_Step);

GB_x0 = GB_x(x0,y0,z0);
GB_y0 = GB_y(x0,y0,z0);
GB_z0 = GB_z(x0,y0,z0);


for i=1:size(x_Value,2)
    for j=1:size(y_Value,2)
        for k=1:size(z_Value,2)
            %if sqrt(x_Value(i)^2+y_Value(j)^2+z_Value(k)^2)<=sphere_radius % in the spehre
                LinearityGB_x(i,j,k) = abs((GB_x(i,j,k)/GB_x0)-1);
                LinearityGB_y(i,j,k) = abs((GB_y(i,j,k)/GB_y0)-1);
                LinearityGB_z(i,j,k) = abs((GB_z(i,j,k)/GB_z0)-1);
                
                l = l+1;
            %end
        end
    end
end


%% Display in XY plane Bx
B2 = LinearityGB_x;
figure('Name','Gradient according to x')
subplot(3,2,1)
B3 = squeeze(B2(:,:,z0));
[~,h] = contour(y_Value,x_Value,B3,procentPlane);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)

%plot a circle to show the coil radius
point = 0;
clear ('circleX','circleY');
for i=1:size(xp_Value,2)
    for j=1:size(yp_Value,2)
        if sqrt(xp_Value(i)^2+yp_Value(j)^2) == coil_radius% 
            point = point +1;
            circleX(point) = xp_Value(i);
            circleY(point) = yp_Value(j);
        end
    end
end

hold all;
plot(circleX,circleY,'*');
axis equal
set(gca,'XTick',y_Start:(y_Stop-y_Start)/6:y_Stop)
set(gca,'YTick',x_Start:(x_Stop-x_Start)/6:x_Stop)
grid on
xlabel('Y axis')
ylabel('X axis')
xlim([-coil_radius coil_radius])
ylim([-coil_radius coil_radius])

% YZ plane
%figure
subplot(3,2,[3 4])
B3 = squeeze(B2(x0,:,:));
%B3 = B3'; % To have the Z axis in a more "natural" direction
[~,h] = contour(z_Value,y_Value,B3,procentPlane);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)

%plot a rectangle to show the coil
point = 0;
clear ('circleX','circleY');
for i=1:size(zp_Value,2)
    for j=1:size(yp_Value,2)
        if sqrt(yp_Value(j)^2) == coil_radius% 
            point = point +1;
            circleX(point) = zp_Value(i);
            circleY(point) = yp_Value(j);
        end
    end
end

hold all;
plot(circleX,circleY,'*');
axis equal
set(gca,'XTick',z_Start:(z_Stop-z_Start)/6:z_Stop)
set(gca,'YTick',y_Start:(y_Stop-y_Start)/6:y_Stop)
grid on
xlabel('Z axis')
ylabel('Y axis')
xlim([-coil_length/2 coil_length/2])
ylim([-coil_radius coil_radius])


% XZ plane
%figure
subplot(3,2,[5 6])
B3 = squeeze(B2(:,y0,:));
[~,h] = contour(z_Value,x_Value,B3,procentPlane);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)

%plot a rectangle to show the coil
point = 0;
clear ('circleX','circleY');
for i=1:size(zp_Value,2)
    for j=1:size(yp_Value,2)
        if sqrt(yp_Value(j)^2) == coil_radius% 
            point = point +1;
            circleX(point) = zp_Value(i);
            circleY(point) = yp_Value(j);
        end
    end
end

hold all;
plot(circleX,circleY,'*');
axis equal
set(gca,'XTick',z_Start:(z_Stop-z_Start)/6:z_Stop)
set(gca,'YTick',x_Start:(x_Stop-x_Start)/6:x_Stop)
grid on
xlabel('Z axis')
ylabel('X axis')
xlim([-coil_length/2 coil_length/2])
ylim([-coil_radius coil_radius])

% 3D
%figure
subplot(3,2,2)
isosurface(y_Value,x_Value,z_Value,B2,procentVolume)
if isstruct(wire) % If we provide the wire geometrie, display it.
    hold all;
    axis equal
    for i = 1 :size(wire,2)
        if wire(i).currentDirection == 1
            plot3(wire(i).Coord(2,:),wire(i).Coord(1,:),wire(i).Coord(3,:),'blue');
        elseif wire(i).currentDirection == -1
            plot3(wire(i).Coord(2,:),wire(i).Coord(1,:),wire(i).Coord(3,:),'red');
        end
    end
end
set(gca,'XTick',y_Start:(y_Stop-y_Start)/6:y_Stop)
set(gca,'YTick',x_Start:(x_Stop-x_Start)/6:x_Stop)
set(gca,'ZTick',z_Start:(z_Stop-z_Start)/6:z_Stop)
xlim([-coil_radius coil_radius])
ylim([-coil_radius coil_radius])
axis equal
xlabel('Y axis')
ylabel('X axis')
zlabel('Z axis')

%% Display in XY plane Bx
B2 = LinearityGB_y;
figure('Name','Gradient according to y')
subplot(3,2,1)
B3 = squeeze(B2(:,:,z0));
[~,h] = contour(y_Value,x_Value,B3,procentPlane);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)

%plot a circle to show the coil radius
point = 0;
clear ('circleX','circleY');
for i=1:size(xp_Value,2)
    for j=1:size(yp_Value,2)
        if sqrt(xp_Value(i)^2+yp_Value(j)^2) == coil_radius% 
            point = point +1;
            circleX(point) = xp_Value(i);
            circleY(point) = yp_Value(j);
        end
    end
end

hold all;
plot(circleX,circleY,'*');
axis equal
set(gca,'XTick',y_Start:(y_Stop-y_Start)/6:y_Stop)
set(gca,'YTick',x_Start:(x_Stop-x_Start)/6:x_Stop)
grid on
xlabel('Y axis')
ylabel('X axis')
xlim([-coil_radius coil_radius])
ylim([-coil_radius coil_radius])

% YZ plane
%figure
subplot(3,2,[3 4])
B3 = squeeze(B2(x0,:,:));
%B3 = B3'; % To have the Z axis in a more "natural" direction
[~,h] = contour(z_Value,y_Value,B3,procentPlane);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)

%plot a rectangle to show the coil
point = 0;
clear ('circleX','circleY');
for i=1:size(zp_Value,2)
    for j=1:size(yp_Value,2)
        if sqrt(yp_Value(j)^2) == coil_radius% 
            point = point +1;
            circleX(point) = zp_Value(i);
            circleY(point) = yp_Value(j);
        end
    end
end

hold all;
plot(circleX,circleY,'*');
axis equal
set(gca,'XTick',z_Start:(z_Stop-z_Start)/6:z_Stop)
set(gca,'YTick',y_Start:(y_Stop-y_Start)/6:y_Stop)
grid on
xlabel('Z axis')
ylabel('Y axis')
xlim([-coil_length/2 coil_length/2])
ylim([-coil_radius coil_radius])


% XZ plane
%figure
subplot(3,2,[5 6])
B3 = squeeze(B2(:,y0,:));
[~,h] = contour(z_Value,x_Value,B3,procentPlane);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)

%plot a rectangle to show the coil
point = 0;
clear ('circleX','circleY');
for i=1:size(zp_Value,2)
    for j=1:size(yp_Value,2)
        if sqrt(yp_Value(j)^2) == coil_radius% 
            point = point +1;
            circleX(point) = zp_Value(i);
            circleY(point) = yp_Value(j);
        end
    end
end

hold all;
plot(circleX,circleY,'*');
axis equal
set(gca,'XTick',z_Start:(z_Stop-z_Start)/6:z_Stop)
set(gca,'YTick',x_Start:(x_Stop-x_Start)/6:x_Stop)
grid on
xlabel('Z axis')
ylabel('X axis')
xlim([-coil_length/2 coil_length/2])
ylim([-coil_radius coil_radius])

% 3D
%figure
subplot(3,2,2)
isosurface(y_Value,x_Value,z_Value,B2,procentVolume)
if isstruct(wire) % If we provide the wire geometrie, display it.
    hold all;
    axis equal
    for i = 1 :size(wire,2)
        if wire(i).currentDirection == 1
            plot3(wire(i).Coord(2,:),wire(i).Coord(1,:),wire(i).Coord(3,:),'blue');
        elseif wire(i).currentDirection == -1
            plot3(wire(i).Coord(2,:),wire(i).Coord(1,:),wire(i).Coord(3,:),'red');
        end
    end
end
set(gca,'XTick',y_Start:(y_Stop-y_Start)/6:y_Stop)
set(gca,'YTick',x_Start:(x_Stop-x_Start)/6:x_Stop)
set(gca,'ZTick',z_Start:(z_Stop-z_Start)/6:z_Stop)
xlim([-coil_radius coil_radius])
ylim([-coil_radius coil_radius])
axis equal
xlabel('Y axis')
ylabel('X axis')
zlabel('Z axis')

%% Display in XY plane Bx
B2 = LinearityGB_z;
figure('Name','Gradient according to z')
subplot(3,2,1)
B3 = squeeze(B2(:,:,z0));
[~,h] = contour(y_Value,x_Value,B3,procentPlane);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)

%plot a circle to show the coil radius
point = 0;
clear ('circleX','circleY');
for i=1:size(xp_Value,2)
    for j=1:size(yp_Value,2)
        if sqrt(xp_Value(i)^2+yp_Value(j)^2) == coil_radius% 
            point = point +1;
            circleX(point) = xp_Value(i);
            circleY(point) = yp_Value(j);
        end
    end
end

hold all;
plot(circleX,circleY,'*');
axis equal
set(gca,'XTick',y_Start:(y_Stop-y_Start)/6:y_Stop)
set(gca,'YTick',x_Start:(x_Stop-x_Start)/6:x_Stop)
grid on
xlabel('Y axis')
ylabel('X axis')
xlim([-coil_radius coil_radius])
ylim([-coil_radius coil_radius])

% YZ plane
%figure
subplot(3,2,[3 4])
B3 = squeeze(B2(x0,:,:));
%B3 = B3'; % To have the Z axis in a more "natural" direction
[~,h] = contour(z_Value,y_Value,B3,procentPlane);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)

%plot a rectangle to show the coil
point = 0;
clear ('circleX','circleY');
for i=1:size(zp_Value,2)
    for j=1:size(yp_Value,2)
        if sqrt(yp_Value(j)^2) == coil_radius% 
            point = point +1;
            circleX(point) = zp_Value(i);
            circleY(point) = yp_Value(j);
        end
    end
end

hold all;
plot(circleX,circleY,'*');
axis equal
set(gca,'XTick',z_Start:(z_Stop-z_Start)/6:z_Stop)
set(gca,'YTick',y_Start:(y_Stop-y_Start)/6:y_Stop)
grid on
xlabel('Z axis')
ylabel('Y axis')
xlim([-coil_length/2 coil_length/2])
ylim([-coil_radius coil_radius])


% XZ plane
%figure
subplot(3,2,[5 6])
B3 = squeeze(B2(:,y0,:));
[~,h] = contour(z_Value,x_Value,B3,procentPlane);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)

%plot a rectangle to show the coil
point = 0;
clear ('circleX','circleY');
for i=1:size(zp_Value,2)
    for j=1:size(yp_Value,2)
        if sqrt(yp_Value(j)^2) == coil_radius% 
            point = point +1;
            circleX(point) = zp_Value(i);
            circleY(point) = yp_Value(j);
        end
    end
end

hold all;
plot(circleX,circleY,'*');
axis equal
set(gca,'XTick',z_Start:(z_Stop-z_Start)/6:z_Stop)
set(gca,'YTick',x_Start:(x_Stop-x_Start)/6:x_Stop)
grid on
xlabel('Z axis')
ylabel('X axis')
xlim([-coil_length/2 coil_length/2])
ylim([-coil_radius coil_radius])

% 3D
%figure
subplot(3,2,2)
isosurface(y_Value,x_Value,z_Value,B2,procentVolume)
if isstruct(wire) % If we provide the wire geometrie, display it.
    hold all;
    axis equal
    for i = 1 :size(wire,2)
        if wire(i).currentDirection == 1
            plot3(wire(i).Coord(2,:),wire(i).Coord(1,:),wire(i).Coord(3,:),'blue');
        elseif wire(i).currentDirection == -1
            plot3(wire(i).Coord(2,:),wire(i).Coord(1,:),wire(i).Coord(3,:),'red');
        end
    end
end
set(gca,'XTick',y_Start:(y_Stop-y_Start)/6:y_Stop)
set(gca,'YTick',x_Start:(x_Stop-x_Start)/6:x_Stop)
set(gca,'ZTick',z_Start:(z_Stop-z_Start)/6:z_Stop)
xlim([-coil_radius coil_radius])
ylim([-coil_radius coil_radius])
axis equal
xlabel('Y axis')
ylabel('X axis')
zlabel('Z axis')
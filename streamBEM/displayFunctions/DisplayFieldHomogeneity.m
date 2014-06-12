function [] = DisplayFieldHomogeneity(B,wire,coil_radius,coil_length,procentPlane,procentVolume,x_Value,y_Value,z_Value)


%procent = [0.01,0.03,0.05,0.08,0.10,0.15];
% B = Bx;
% B = By;
% B = Bz;

x_Start =x_Value(1);
x_Stop = x_Value(size(x_Value,2));
% Y is a position, in meter
y_Start =y_Value(1);
y_Stop = y_Value(size(y_Value,2));
% Z is a position, in meter
z_Start =z_Value(1);
z_Stop = z_Value(size(z_Value,2));

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
zp_Start = -coil_length;
zp_Stop = coil_length;
zp_Step = coil_length/20;
% X is a position, in meter
xp_Value = xp_Start:xp_Step:xp_Stop;
% Y is a position, in meter
yp_Value = yp_Start:yp_Step:yp_Stop;
% Z is a position, in meter
zp_Value = zp_Start:zp_Step:zp_Stop;


% XY plane
figure
subplot(3,2,1)
B2 = abs((B/B(x0,y0,z0))-1);
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
set(gca,'XTick',y_Start:0.02:y_Stop)
set(gca,'YTick',x_Start:0.02:x_Stop)
grid on
xlabel('Y axis')
ylabel('X axis')
xlim([-coil_radius coil_radius])
ylim([-coil_radius coil_radius])

% YZ plane
%figure
subplot(3,2,[3 4])
B2 = abs((B/B(x0,y0,z0))-1);
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
%axis equal
set(gca,'XTick',z_Start:0.02:z_Stop)
set(gca,'YTick',y_Start:0.02:y_Stop)
grid on
xlabel('Z axis')
ylabel('Y axis')
xlim([-coil_length coil_length])
ylim([-coil_radius coil_radius])


% XZ plane
%figure
subplot(3,2,[5 6])
B2 = abs((B/B(x0,y0,z0))-1);
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
%axis equal
set(gca,'XTick',z_Start:0.02:z_Stop)
set(gca,'YTick',x_Start:0.02:x_Stop)
grid on
xlabel('Z axis')
ylabel('X axis')
xlim([-coil_length coil_length])
ylim([-coil_radius coil_radius])

% 3D
%figure
subplot(3,2,2)
B2 = abs((B/B(x0,y0,z0))-1);
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
set(gca,'XTick',y_Start:0.05:y_Stop)
set(gca,'YTick',x_Start:0.05:x_Stop)
set(gca,'ZTick',z_Start:0.05:z_Stop)
xlim([-coil_radius coil_radius])
ylim([-coil_radius coil_radius])
axis equal
xlabel('Y axis')
ylabel('X axis')
zlabel('Z axis')
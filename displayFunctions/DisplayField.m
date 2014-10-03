function [] = DisplayField(Bx,By,Bz,x_Value,y_Value,z_Value,minimumFieldAmplitude,maximumFieldAmplitude)

% This file aim at displaying the field in 3 plane, in order to check the value 

sizeX = size(Bx,1);
sizeY = size(Bx,2);
sizeZ = size(Bx,3);


BX_xy_z0 = squeeze(Bx(:,:,(sizeZ-1)/2 +1));
BY_xy_z0 = squeeze(By(:,:,(sizeZ-1)/2 +1));
BZ_xy_z0 = squeeze(Bz(:,:,(sizeZ-1)/2 +1));

BX_xz_y0 = squeeze(Bx(:,(sizeY-1)/2 +1,:));
BY_xz_y0 = squeeze(By(:,(sizeY-1)/2 +1,:));
BZ_xz_y0 = squeeze(Bz(:,(sizeY-1)/2 +1,:));

BX_yz_x0 = squeeze(Bx((sizeX-1)/2 +1,:,:));
BY_yz_x0 = squeeze(By((sizeX-1)/2 +1,:,:));
BZ_yz_x0 = squeeze(Bz((sizeX-1)/2 +1,:,:));

clims = [minimumFieldAmplitude maximumFieldAmplitude]; % Display the graph for +/- FieldAmplitude

%% Please, be cautious with the order of the axis !
% imagesc don't scheck the size of the x and y vector
% It is inversed, as often with matlab

figure
subplot(3,3,1)
imagesc(y_Value,x_Value,BX_xy_z0,clims);
xlabel('Y axes');
ylabel('X axes');
title('BX xy z0');
axis('square')
subplot(3,3,2)
imagesc(z_Value,x_Value,BX_xz_y0,clims);
xlabel('Z axes');
ylabel('X axes');
title('BX xz y0');
axis('square')
subplot(3,3,3)
imagesc(z_Value,y_Value,BX_yz_x0,clims);
xlabel('Z axes');
ylabel('Y axes');
title('BX yz x0');
axis('square')


subplot(3,3,4)
imagesc(y_Value,x_Value,BY_xy_z0,clims);
xlabel('Y axes');
ylabel('X axes');
title('BY xy z0');
axis('square')
subplot(3,3,5)
imagesc(z_Value,x_Value,BY_xz_y0,clims);
xlabel('Z axes');
ylabel('X axes');
title('BY xz y0');
axis('square')
subplot(3,3,6)
imagesc(z_Value,y_Value,BY_yz_x0,clims);
xlabel('Z axes');
ylabel('Y axes');
title('BY yz x0');
axis('square')

subplot(3,3,7)
imagesc(y_Value,x_Value,BZ_xy_z0,clims);
xlabel('Y axes');
ylabel('X axes');
title('BZ xy z0');
axis('square')
subplot(3,3,8)
imagesc(z_Value,x_Value,BZ_xz_y0,clims);
xlabel('Z axes');
ylabel('x axes');
title('BZ xz y0');
axis('square')
subplot(3,3,9)
imagesc(z_Value,y_Value,BZ_yz_x0,clims);
xlabel('Z axes');
ylabel('Y axes');
title('BZ yz x0');
axis('square')

colormap('gray')

%% Please, be cautious with the order of the axis !

if size(z_Value,2)~= 1
    [X,Y,Z] = meshgrid(y_Value,x_Value,z_Value);
    figure
    subplot(1,3,1)
    slice(X,Y,Z,Bx,y_Value((sizeY-1)/2 +1),x_Value((sizeX-1)/2 +1),z_Value((sizeZ-1)/2 +1))
    xlabel('Y axes');
    ylabel('X axes');
    zlabel('Z axes');
    title('Bx');
    axis('square')
    colormap('gray')

    subplot(1,3,2)
    slice(X,Y,Z,By,y_Value((sizeY-1)/2 +1),x_Value((sizeX-1)/2 +1),z_Value((sizeZ-1)/2 +1))
    xlabel('Y axes');
    ylabel('X axes');
    zlabel('Z axes');
    title('By');
    axis('square')
    colormap('gray')

    subplot(1,3,3)
    slice(X,Y,Z,Bz,y_Value((sizeY-1)/2 +1),x_Value((sizeX-1)/2 +1),z_Value((sizeZ-1)/2 +1))
    xlabel('Y axes');
    ylabel('X axes');
    zlabel('Z axes');
    title('Bz');
    axis('square')
    colormap('gray')
end
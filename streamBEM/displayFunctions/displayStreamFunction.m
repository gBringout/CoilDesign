function [] = displayStreamFunction(triangle,s,node)

figure

%First approximation based on the mean on each triangle
TriAmpl = zeros(size(triangle,1),1);
for i=1:size(triangle,1)
    TriAmpl(i) = (s(triangle(i,1))+s(triangle(i,2))+s(triangle(i,3)))/3;
end


trisurf(triangle,node(:,1),node(:,2),node(:,3),TriAmpl)
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');
axis equal;
%caxis auto
%caxis([-2e5 2e5])

v = caxis;
caxis([-max(abs(v)) max(abs(v))]);
cmap = colormap;
newColorMap = gaelColormap(cmap);
colormap(newColorMap);
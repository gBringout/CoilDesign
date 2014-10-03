function [aspectRatio] = checkMeshAspectRatio(node, triangle,meshName)
% based on an idea from 
%An efficient method for evaluating BEM singular integrals on curved elements
%with application in acoustic analysis
%Junjie Rong, LihuaWen?, Jinyou Xiao

% Check aspect ratio, wehich is define as the longer length divided by the
% smallest one in the triangle

for i=1:size(triangle,2)
    a = node(triangle(i).node(1)).coord;
    b = node(triangle(i).node(2)).coord;
    c = node(triangle(i).node(3)).coord;
    sideLength(1) = norm(a-b);
    sideLength(2) = norm(a-c);
    sideLength(3) = norm(c-b);
    aspectRatio(i) = max(sideLength)/min(sideLength);
end

% figure
% plot(aspectRatio) % Study 1/A(theta to see that)
fprintf('Mesh: %s. \n Mean AR: %f. Std: %f\n',meshName,mean(aspectRatio),std(aspectRatio))
if max(aspectRatio>1.5)
    fprintf(2,'The mesh: %s is probably not good to evaluate singular and hypersingular function.\n Aspect ratio superior to 1.5 detected\n',meshName);
end

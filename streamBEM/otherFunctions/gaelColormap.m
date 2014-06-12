function [cmap2] = gaelColormap(colormap)
% This colormap  aim to have :
%   a red color for the max value
%   a blue color for the min value
%   a white color for the middle value

%
if mod(size(colormap,1),2) == 0 % if the number of color is even, we have 2 values for the white
    sizeCmap = size(colormap,1);
    for k = 1:sizeCmap
        if k<sizeCmap/2
            %add the red color
            % The first point fo through 
            x1=1;% to have 
            y1=0;
            % The second point go through 
            x2=sizeCmap/2;% to have 
            y2=1;
            %We know calculate the line going through those 2 point
            a = (y1-y2)/(x1-x2);
            b = y2-a*x2;
            % and use it to fill the color space
            cmap2(k,1) = 1;
            cmap2(k,2) = a*k + b;
            cmap2(k,3) = a*k + b;

        elseif k==sizeCmap/2 || k==sizeCmap/2+1
            cmap2(k,1) = 1;
            cmap2(k,2) = 1;
            cmap2(k,3) = 1;
        else
            %add the blue color
            % The first point fo through 
            x1=sizeCmap;% to have 
            y1=0;
            % The second point go through 
            x2=sizeCmap/2+1;% to have 
            y2=1;
            %We know calculate the line going through those 2 point
            a = (y1-y2)/(x1-x2);
            b = y2-a*x2;
            cmap2(k,1) = a*k + b;
            cmap2(k,2) = a*k + b;
            cmap2(k,3) = 1;
        end
    end
else 
    if sizeCmap == 1
        sizeCmap = 1;
    else
        sizeCmap = size(colormap,1)-1;
    end
    for k = 1:sizeCmap
        if k<=sizeCmap/2
            %add the red color
            % The first point fo through 
            x1=1;% to have 
            y1=0;
            % The second point go through 
            x2=sizeCmap/2+1;% to have 
            y2=1;
            %We know calculate the line going through those 2 point
            a = (y1-y2)/(x1-x2);
            b = y2-a*x2;
            % and use it to fill the color space
            cmap2(k,1) = 1;
            cmap2(k,2) = a*k + b;
            cmap2(k,3) = a*k + b;

        elseif k==sizeCmap/2+1
            cmap2(k,1) = 1;
            cmap2(k,2) = 1;
            cmap2(k,3) = 1;
        else
            %add the blue color
            % The first point fo through 
            x1=sizeCmap+1;% to have 
            y1=0;
            % The second point go through 
            x2=sizeCmap/2+1;% to have 
            y2=1;
            %We know calculate the line going through those 2 point
            a = (y1-y2)/(x1-x2);
            b = y2-a*x2;
            cmap2(k,1) = a*k + b;
            cmap2(k,2) = a*k + b;
            cmap2(k,3) = 1;
        end
    end
end

cmap2 = flipud(cmap2); % have to turn to to get the red as +1 and the blut at -1;
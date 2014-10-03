function [node2,triangle2,subBoundary] = NodeSorting2(node,triangle)
% Type = 1: we try to sort the node
% Type = 2: we do not try to sort the node
if nargin<3
    type = 1;
end

if type == 1
    %% First we search how much border we have
    % we should perhaps go for the vector.
    % A vector which belongs to a single triangle is a border

    % calculate the basis
    % triangle = coil.listTriangle;
    % node = coil.listNode;
    % triangle = shield.listTriangle;
    % node = shield.listNode;

    % for each vector
    % nbr triangle containing this vector
    l = 1;
    isBoundary = [];
    for i=1:size(node,1)
        for j=1:size(node,1)
            %now find all the triangle with both node A and B
            % first we look at the list with one
            [row1,~] = find(triangle == i);
            % then look for the second
            [row2,~] = find(triangle(row1,:) == j);
           if size(row2,1) == 1 %here we found one border
               % and we put the smaller index on the left
                if i<j
                    isBoundary(l,:) = [i j];
                else
                    isBoundary(l,:) = [j i];
                end
                l = l+1;
           end
        end
    end
    %remove the duplicate rows
    isBoundary = unique(isBoundary, 'rows');
    % %%
    isBoundarySave = isBoundary;
    %% Look for the sub-boundaries
    %isBoundary = isBoundarySave;
    nodeOnBoundary = unique(isBoundary(:));
    if not(isempty(nodeOnBoundary))
        % Inizialize
        subBoundary(1,1).node(1,1) = nodeOnBoundary(1);
        [row, col] = find(isBoundary == nodeOnBoundary(1));
        nodeOnBoundary(1) = [];
        i=1; % we start on the first sub-boundrie
        j=2; % we have to first fill the second point
        k=3; % next points in the target list. we always start with 1 node and
        % find 2 on the first search. 1+2 = 3

        while not(isempty(isBoundary)) %until we have no node on the borders left
            while not(isempty(row))%as long as we find nodes connected to the actual sub-boundary

                for l=1:size(row,1)% we add the new point
                    %only the first time we should find 2 row
                    % then it should be only one
                    subBoundary(i,1).node(j,1) = isBoundary(row(l),mod(col(l),2)+1);
                    [row2,~] = find(subBoundary(i,1).node(j,1) == nodeOnBoundary);
                    nodeOnBoundary(row2) = [];
                    j = j+1;
                end
                isBoundary(row,:) = []; %we remove the use row
                [row, col] = find(isBoundary == subBoundary(i,1).node(k,1));
                k=k+1;
            end

            subBoundary(i,1).node = unique(subBoundary(i,1).node);
            if not(isempty(isBoundary))
                i= i+1; % number of sub boundaries
                [row, col] = find(isBoundary == nodeOnBoundary(1));
                subBoundary(i,1).node(1,1) = nodeOnBoundary(1);
                nodeOnBoundary(1) = [];
                j=2;
                k=3;
            end
        end


        %% Then we create a new node matrix with the border node on top of it
        isBoundary = isBoundarySave;

        node2 = zeros(size(node,1),3);
        index = 1;
        nodeListe = 1:size(node,1);
        nodeConversionTable = 0; %conversion from the old list Node to the new one
        % newNode = nodeConversionTable(oldNode);

        % we first put the border node
        for j=1:size(subBoundary,1)
            for k=1:size(subBoundary(j).node,1)
                node2(index,:) = node(subBoundary(j).node(k),:);
                nodeConversionTable(subBoundary(j).node(k)) = index;
                row = find(subBoundary(j).node(k) == nodeListe);
                nodeListe(row) = [];
                index = index + 1;
            end
        end
        index2 = 1;
        for i=index:size(node,1)
            node2(i,:) = node(nodeListe(index2),:);
            nodeConversionTable(nodeListe(index2)) = i;
            index2 = index2+1;
        end
        % to check if it is working correctly
        % t = unique(nodeConversionTable);

        %% Then we create a new triangle matrix, with the converted node

        for i=1:size(triangle,1)
            triangle2(i,1) = nodeConversionTable(triangle(i,1));
            triangle2(i,2) = nodeConversionTable(triangle(i,2));
            triangle2(i,3) = nodeConversionTable(triangle(i,3));
        end
    else
        node2 = node;
        triangle2 = triangle;
        subBoundary = [];
    end
end
    
function [wire,nbrwire,minimalWireSpacing,res] = NbrWireOptimisation3(node,triangle,s,spaceBetweenWire,startingWireNumber,distanceBetween2Wire,rateIncreasingWire)

% spaceBetweenWire is the maximal minimal distance between wire in meter.
res(1).minimalWireSpacing.value = 1;
res(1).nbrwire = startingWireNumber-rateIncreasingWire;

index = 2;
while res(index-1).minimalWireSpacing.value>spaceBetweenWire
    tic
    res(index).nbrwire = res(index-1).nbrwire+rateIncreasingWire;
    res(index).wire = exctracteWire5(node,triangle,s,res(index).nbrwire,distanceBetween2Wire);
    res(index).wire = cleanWire(res(index).wire);
    res(index).minimalWireSpacing = distanceBetweenWire(res(index).wire);
    fprintf('Nb Level: %i, Minimal Spacing: %2.1i, Time spent: %5.0f.\n',res(index).nbrwire,res(index).minimalWireSpacing.value, toc);
    index = index+1;
end

% As the criteria has been overpass, we recalculate the previous one
nbrwire = res(index-2).nbrwire;
wire = res(index-2).wire;
minimalWireSpacing = res(index-2).minimalWireSpacing;
if index>3
    minimalWireSpacing.oldValue = res(index-3).minimalWireSpacing;
else
    minimalWireSpacing.oldValue = 0;
end
%displayWire(wire);
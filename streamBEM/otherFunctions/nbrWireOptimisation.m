function [Wire,nbrwire,minimalWireSpacing] = nbrWireOptimisation(Node,Triangle,I,spaceBetweenWire)

% spaceBetweenWire is the maximal minimal distance between wire in meter.
minimalWireSpacing.value = 1;
nbrwire = 10;
while minimalWireSpacing.value>spaceBetweenWire
    nbrwire = nbrwire+2;
    Wire = exctracteWire2(Node,Triangle,I,nbrwire);
    minimalWireSpacing = DistanceBetweenWire(Wire);
end
%displayWire(Wire);
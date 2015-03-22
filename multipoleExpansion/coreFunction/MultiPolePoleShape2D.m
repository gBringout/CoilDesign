function [r] = MultiPolePoleShape2D(r0,n,phi)
% r0 : free inner radius
% n order of the pole shape
% phi:phase in the symetry of the pole

% Equation of the pole shape according to
% page 15 (pdf 23) from cern-2010-004
% r^n = r0^n/sin(n*theta-phi)

theta = 0:2*pi/360:2*pi;
r = zeros(size(theta,2),1);

for i=1:size(theta,2)
    r(i)= nthroot(abs(r0^n/sin(n*theta(i)-phi)),n);
    if r(i)>4*r0
        r(i)=4*r0;
    end
end
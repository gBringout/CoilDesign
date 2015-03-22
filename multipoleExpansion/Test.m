addpath(genpath(fullfile('..')))

% creation of the magnet's shape of a pure multipole field

n = 2; % for a quadrupole
phi = 0;
r0 = 1; % with inner radius 1 meter
theta = 0:2*pi/360:2*pi;

r = MultiPolePoleShape2D(r0,n,phi);

figure
polar(theta,r')
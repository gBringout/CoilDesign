function [u,v,ck] = triGaussPoints(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function TriGaussPoints provides the Gaussian points and weights %
% for the Gaussian quadrature of order n for the standard triangles. %
% based on the Paper of Rathod and al. "Gauss Legendre quadrature over a triangle"
% Input: n - the order of the Gaussian quadrature (n<=12) %
% %
% Output: xw - a n by 3 matrix: %
% 1st column gives the x-coordinates of points %
% 2nd column gives the y-coordinates of points %
% 3rd column gives the weights %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% u = zeros(n^2,1);
% v = zeros(n^2,1);
% ck = zeros(n^2,1);
[eta,w] = gauss(n);

k=1;
for i=1:n
    for j=1:n
        u(k,1) = (1+eta(i))/2;
        v(k,1) = (1-eta(i))*(1+eta(j))/4;
        ck(k,1) = ((1-eta(i))/8)*w(i)*w(j);
        k=k+1;
    end
end

%% from http://www.abcm.org.br/pt/wp-content/anais/cobem/2007/pdf/COBEM2007-1614.pdf
% % N = 3
% 
% u(1,1) = 0.1063508;
% u(2,1) = 0.4718246;
% u(3,1) = 0.8372983;
% u(4,1) = 0.08452624;
% u(5,1) = 0.3750000;
% u(6,1) = 0.6654738;
% u(7,1) = 0.06270166;
% u(8,1) = 0.2781754;
% u(9,1) = 0.4936492;
% 
% v(1,1) = 0.1063508;
% v(2,1) = 0.08452624;
% v(3,1) = 0.06270166;
% v(4,1) = 0.4718246;
% v(5,1) = 0.3750000;
% v(6,1) = 0.2781754;
% v(7,1) = 0.8372983;
% v(8,1) = 0.6654738;
% v(9,1) = 0.4936492;
%  
% ck(1,1) = 0.06846439;
% ck(2,1) = 0.08563571;
% ck(3,1) = 0.03858025;
% ck(4,1) = 0.08563571;
% ck(5,1) = 0.09876543;
% ck(6,1) = 0.03782109;
% ck(7,1) = 0.03858025;
% ck(8,1) = 0.03782109;
% ck(9,1) = 0.008696116;


%% Generalized from
% http://www.abcm.org.br/pt/wp-content/anais/cobem/2007/pdf/COBEM2007-1614.pdf
% u = zeros(n^2,1);
% v = zeros(n^2,1);
% ck = zeros(n^2,1);
% [eta,w] = gauss(n);
% 
% k=1;
% for i=1:n
%     for j=1:n
%         u(k) = (1-(1+eta(j))/4)*((1+eta(i))/2);
%         v(k) = (1-(1+eta(i))/4)*((1+eta(j))/2);
%         ck(k) = 1/4*(1-((2+eta(i)+eta(j))/(4)))*w(i)*w(j);
%         k=k+1;
%     end
% end

%% Generalized form from the paper from Rathod "Gauss Legendre quadrature
% over a triangle"


%% Based on Poole thesis N = 2
% It seams to come from the paper from Rathod "Gauss Legendre quadrature
% over a triangle"
% u(1,1) = 0.211324865;
% u(2,1) = 0.211324865;
% u(3,1) = 0.788675134;
% u(4,1) = 0.788675134;
% 
% v(1,1) = 0.166666667;
% v(2,1) = 0.622008467;
% v(3,1) = 0.044658198;
% v(4,1) = 0.166666667;
% 
% ck(1,1) = 0.197168783;
% ck(2,1) = 0.197168783;
% ck(3,1) = 0.052831216;
% ck(4,1) = 0.052831216;

%% From Poole thesis N = 3
% It seams to come from the paper from Rathod "Gauss Legendre quadrature
% over a triangle"
% u(1,1) = 0.112701665; % Gauss-Legendre Integration quadrature coordinate
% u(2,1) = 0.112701665;
% u(3,1) = 0.112701665;
% u(4,1) = 0.500000000;
% u(5,1) = 0.500000000;
% u(6,1) = 0.500000000;
% u(7,1) = 0.887298334;
% u(8,1) = 0.887298334;
% u(9,1) = 0.887298334;
% 
% v(1,1) = 0.100000000;
% v(2,1) = 0.443649167;
% v(3,1) = 0.787298334;
% v(4,1) = 0.056350832;
% v(5,1) = 0.250000000;
% v(6,1) = 0.443649167;
% v(7,1) = 0.012701665;
% v(8,1) = 0.056350832;
% v(9,1) = 0.100000000;
% 
% ck(1,1) = 0.068464377; % Gauss-Legendre Integration quadrature weightings factor
% ck(2,1) = 0.109543004;
% ck(3,1) = 0.068464377;
% ck(4,1) = 0.061728395;
% ck(5,1) = 0.098765432;
% ck(6,1) = 0.061728395;
% ck(7,1) = 0.008696116;
% ck(8,1) = 0.013913785;
% ck(9,1) = 0.008696116;

%fprintf('\n hey it works! size(u) = [%i %i]\n',size(u,1),size(u,2))



end
function [u,v,ck] = triGaussPoints(n)
% This file aim at calculating the coordinates needed for a Gauss Legendre
% integration method on a triangle and the associated weighting factors. From:
% H. T. Rathod et al., Gauss Legendre quadrature over a triangle, J. Indian
% Inst. Sci, 2004, vol 84, pages 183-188
%
% n : order of the foreseen integration
% u,v : coordinate of the point
% ck : weighting coefficients

[eta,w] = gauss(n);

k=1;
for i=1:size(eta,1)
    for j=1:size(eta,1)
        u(k,1) = (1+eta(i))/2;
        v(k,1) = (1-eta(i))*(1+eta(j))/4;
        ck(k,1) = ((1-eta(i))/8)*w(i)*w(j);
        k=k+1;
    end
end

end
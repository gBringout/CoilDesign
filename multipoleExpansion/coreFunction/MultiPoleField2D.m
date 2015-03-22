function [Bx,By,Bz] = MultiPoleField2D(x_Value,y_Value,C,n)

Bx = zeros(size(x_Value,2),size(y_Value,2));
By = zeros(size(x_Value,2),size(y_Value,2));
Bz = zeros(size(x_Value,2),size(y_Value,2));


for k=1:size(x_Value,2)
    for j=1:size(y_Value,2)
        Bcomplex = C*(x_Value(k)+1i*y_Value(j))^(n-1);
        Bx(k,j) = imag(Bcomplex);
        By(k,j) = real(Bcomplex);
        Bz(k,j) = 0;
    end
end
        
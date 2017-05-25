function [ Y ] = resultant( X )
%Get Resultant. Size of input vector must be greater than 3.
X1 = X.^2;
Y1 = (X1(:,1)+X1(:,2)+X1(:,3));
Y = Y1.^(1/2);
end


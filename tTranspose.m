function Tran = tTranspose(X)
sX = size(X);
A = X(:,:,1);
B = X(:,:,2);
C = X(:,:,3);
Tran = zeros(sX(2),sX(1),sX(3));
Tran(:,:,1) = A';
Tran(:,:,2)=C';
Tran(:,:,3)=B';
end

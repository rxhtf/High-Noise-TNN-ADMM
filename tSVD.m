function  [U,S,V] = tSVD(X)
X = fft(X,[],3);
sZ = size(X);

U=zeros(sZ(1),sZ(1),sZ(3));
V = zeros(sZ(2),sZ(2),sZ(3));
S=zeros(sZ(1),sZ(2),sZ(3));
for i  = 1:3
   [U(:,:,i),S(:,:,i),V(:,:,i)] = svd(X(:,:,i));
end
U = ifft(U,[],3);
V=ifft(V,[],3);
S=ifft(S,[],3);

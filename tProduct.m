function C = tProduct(A,B)
A = fft(A,[],3);
B = fft(B,[],3);
sA = size(A);
sB = size(B);
C= zeros(sA(1),sB(2),sA(3));
for i = 1:3
    C(:,:,i)=A(:,:,i)*B(:,:,i);
end
C = ifft(C,[],3);
end

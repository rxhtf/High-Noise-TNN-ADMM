function [n,sum1,sum2] = TNNmap2(X,ind)
Xh = X;
s = size(X);
p = length(s);
%N = prod(s(3:p));

if ~exist('ind','var')
    [~,ind] = buildmap2(s);
end

for i = 3:p
    Xh = fft(Xh,[],i);
end

m1 = min(s(1:2));

S1 = zeros([m1 1 length(ind{1})]);
S2 = zeros([m1 1 length(ind{2})]);

for k = 1:length(ind{1}) %calculate SVD for these slices
    S1(:,1,k) = svd(Xh(:,:,ind{1}(k)));
end

for k = 1:length(ind{2}) %calculate SVD for these slices
    S2(:,1,k) = svd(Xh(:,:,ind{2}(k)));
end

sum1 = sum(S1(:));
sum2 = sum(S2(:));
n = sum1+2*sum2;

%[~,S,~] = tSVDFmap2(X,ind,'diag');
%n = sum(S(:));

end


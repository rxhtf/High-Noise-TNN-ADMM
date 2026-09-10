%computes the Fro norm of a p-order tensor A
function y = tNormF(A)
% [~,S,~] = tSVD(A);
% s = size(S);
% p = length(s);
% N = prod(s(3:p));
% y = 0;
% for i = 1:N
%     y = y + sum(diag(S(:,:,i)).^2); %sum(diag(A(:,:,i)).^2);
% end
% y = y^0.5;

N = numel(A);
y = 0;
for i = 1:N
    y = y + abs(A(i))^2;
end
y = y^0.5;


end


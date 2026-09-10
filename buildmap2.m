%builds map of conjugate slices for a given set of dimensions
%ind is a partition of indicies:
%ind{1} is list of indicies of self conjugate slices
%ind{2} is list of indicies of first slices in pairs
%ind{3} is list of indicies of second slices in pairs
%s is size/dimensions of tensor
function [map,ind] = buildmap2(s)
p = length(s);
N = prod(s(3:end));
map = ones(N,1);
%indi = zeros(1,N); %0 if self-conjugate, 1 is first slice in pair, -1 if second slice in pair

%determine number of even dimensions
eCounter = 0;
for i = 3:p
    if mod(s(i),2) == 0
        eCounter = eCounter+1;
    end
end

%prepare index partition
temp = 2^eCounter;
ind = cell(3,1);
ind{1} = ones(1,temp);
ind{2} = zeros(1,(N-temp)/2);
ind{3} = zeros(1,(N-temp)/2);
iCounter1 = 2; %already know first slice is self conj
iCounter = 1;

%loop over each frontal slice
%exclude first slice since always own conjugate
for k = 2:N
    if map(k) == 1 %has not yet been identified as a conjugate slice
        map(k) = k;
        
        %v = zeros(p-2,1); %used to store the 'coordinates' of the slice
        ks = 0; %var for conjugate slice number
        c = k; %temp var for slice number
        n = N; %temp var for product of dimensions
        for i = p:-1:4
            n = n/s(i); %update product
            %v(i-2) = ceil(c/n); %calculate coordinate
            %ks = ks + mod(s(i)-(v(i-2)-1),s(i))*n; %iterate conjugate slice number
            %c = c - (v(i-2)-1)*n; %update var corresponding with slice number
            
            coor = ceil(c/n); %calculate coordinate
            c = c - (coor-1)*n; %update var corresponding with slice number
            ks = ks + mod(s(i)-(coor-1),s(i))*n; %update conjugate slice number
        end
        %v(3) = c;
        ks = ks + mod(s(3)-(c-1),s(3))+1; %conjugate slice number
        
        map(ks) = k; %identify conjugate slice
        
        if ks == k
            ind{1}(iCounter1) = k;
            iCounter1 = iCounter1+1;
        else
            ind{2}(iCounter) = k;
            ind{3}(iCounter) = ks;
            iCounter = iCounter+1;
        end
    end
end

%temp = map-(1:N)';
%ind1 = find(~temp); %indicies of each first slice
%ind2 = find(temp); %indicies of conjugate slices

%ind0 = find(~indi);
%ind1 = find(indi > 0);
%ind2 = find(indi < 0);

%ind = cell(3,1);
%ind{1} = find(~indi); % == 0
%ind{2} = find(indi > 0);
%ind{3} = find(indi < 0);

end
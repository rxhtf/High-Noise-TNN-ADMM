%Y is data
%A is sampling operator (as a function!)
%T is true, passed to compute residual
%iM is max iteration
%Trying to recover X from observed Y and sampling operator A.
%res is residual vector at each iteration, obj is objective function, GTE
%ground truth error, totalTime is runtime of alg

%v7 uses Hadamard Inverse for denoising (entrywise division)
%v6 groups vars into structs. Allows for rGTE calc and allows for epsilon to be specified.
%v5 reorders the operations to make more sense, i.e. X update last! Z
%unneeded for input and output!
%v4 reduces operations in the standard domain by operating in the Fourier
%domain (use buildmap2/ind)
%v3 assumes no underlying data is known
%v2 returns B,Z for continuation and allows them as inputs
%stores Xold and computes convergence ratio=||X-Xold||/||Xold||
%Uses map (of conjugates) for improved runtime
%Uses econ SVD!
function [contD,AV] = tADMMv7(initD,iterP,contD)%(Y,A,iM,mode,rho,ind,B,X)
%Additional Variables: AV = {res,rGTE(opt),obj,convR,totalTime}

%Initial Data: initD = {Y,A,T(opt),ind}

sY = size(initD.Y);
sm = min(sY(1:2));

p = length(sY);
%N = prod(sY(3:p));

if ~isfield(initD,'ind')
    [~,initD.ind] = buildmap2(sY);
end


%Iteration Parameters: iterP = {iM,mode,rho,epsilon}
if ~isfield(iterP,'mode')
    iterP.mode = 'constrained';
end

if ~isfield(iterP,'rho')
    iterP.rho = 0.0001; %1, 0.1, ..., 0.0001 %rho needs to be > 1 for unconstrained problem!
end
disp(['rho = ',num2str(iterP.rho)])

if strcmp(iterP.mode,'constrained')
    lambda = 1;
elseif strcmp(iterP.mode,'unconstrained')
    %lambda = 0.00001*2; %can vary .001,10,...

    if ~isfield(iterP,'lambda')
        lambda = 1e-4*iterP.rho;
    else
        lambda = iterP.lambda;
    end
    disp(['lambda = ',num2str(lambda)])
    
    unconCon = iterP.rho + initD.A(ones(sY)); %unconstrained constant
    %disp(['unconCon(1) = ',num2str(unconCon(1))]) %for debugging
    %figure
    %imagesc(unconCon(:,:,1))
    
    %check = (iterP.rho + initD.A(ones(sY)))./unconCon - ones(S);
    
end

if ~isfield(iterP,'epsilon')
    iterP.epsilon = 1e-3; %1e-4; %1e-5; %for convergence rate
end

if ~isfield(iterP,'iM')
    iterP.iM = 10;
end
%iterMax = iterP.iM; %10,30,50,100
iter = 1;


%Continuation Data: contD = {B,X}
if ~exist('contD','var')
    contD = struct;
end

if ~isfield(contD,'B')%~exist('B','var')
    contD.B = zeros(sY);
end
if ~isfield(contD,'X')%~exist('X','var')
    contD.X = initD.Y; %zeros(sY);
end


tol = 10^-10; %tolerence for truncating sing values in economic SVD

%additional parameters only needed if updating rho
% mu = 10;
% tinc = 1; %1,2
% tdec = 1; %1,2


AV = struct;
AV.res = zeros(iterP.iM+1,1);
AV.obj = zeros(iterP.iM+1,1);
AV.convR = zeros(iterP.iM+1,1);

if isfield(initD,'T')
    Tn = tNormF(initD.T);
    AV.rGTE = zeros(iterP.iM+1,1);
    AV.rGTE(1) = tNormF(initD.T-contD.X)/Tn;
end

%Initial values
AV.res(1) = tNormF(initD.A(contD.X)-initD.Y); %res(1) = tNormF(A(Y)-Y) = 0;
if strcmp(iterP.mode,'constrained')
    AV.obj(1) = TNNmap2(contD.X,initD.ind);%TNN(X);
elseif strcmp(iterP.mode,'unconstrained')
    AV.objT1 = zeros(iterP.iM+1,1);
    AV.objT1(1) = TNNmap2(contD.X,initD.ind);
    AV.obj(1) = lambda/2*AV.objT1(1)+AV.res(1)^2/2;
end



%display initial iteration info, 'iter 0'
dispStr = ['Iter ', num2str(iter), ': res ', num2str(AV.res(iter)), ', obj ', num2str(AV.obj(iter))];
if isfield(AV,'rGTE')
    dispStr = [dispStr, ', rGTE ', num2str(AV.rGTE(iter))];
end
dispStr = [dispStr, ', Conv ', num2str(AV.convR(iter))];
disp(dispStr)

AV.totalTime = 0;


%set econ sizes of U,S,V
sU = [sY(1) sm];
sV = [sY(2) sm];
sS = [sm sm];

indVec = [initD.ind{2} initD.ind{1}];
lind = length(indVec);

%tic
while iter <= iterP.iM %GTE(iter) > epsilon && iter <= iterMax %problem with GTE(iter)! %for iter = 1:iterMax
    tic
    %Z update
    D = contD.X+contD.B;

    %SVD with shrinking step
    %loop to calculate fft along each dimension
    for i = 3:p
        D = fft(D,[],i);
    end
    %toc
    
    U = zeros([sU lind]);
    V = zeros([sV lind]);
    S = zeros([sS lind]);
    r = 0; %rank
    
    %loop for SVD of each frontal slice
    %Note, for >3D tensors, Matlab's (:,:,i) notation 'unfolds' into 3D, then
    %gets ith slice
    eps = lambda/iterP.rho; %1/rho;
    for i = 1:lind
        [U(:,:,i),S(:,:,i),V(:,:,i)] = svd(D(:,:,indVec(i)),'econ');
        for j = 1:sm
            S(j,j,i) = shrink(S(j,j,i),eps);
        end
        %r update
        if r < sm
            j = r;
            r = sm; %assume slice is full rank and update if not
            while j < sm
                j = j+1;
                if S(j,j,i) < tol
                    r = j-1; %keep 1:j-1, discard j:end since starting from (j,j) sing value is too small
                    j = sm;
                end
            end
        end
    end
    %toc
    
    %cut padding
    if r < sm
        U = U(:,1:r,:);
        S = S(1:r,1:r,:);
        V = V(:,1:r,:);
        %reshape not needed, since only the slices are referenced and their
        %product is stored in the appropriate dimension
    end
    %conj and ifft not needed, stay in Fourier domain for tProduct calc
    %toc
    
    %Zpre = Z; %store previous Z for later calc of rho update
    Z = zeros(sY);
    if sY(1) > sY(2)
        for i = 1:lind
            Z(:,:,indVec(i)) = U(:,:,i)*(S(:,:,i)*V(:,:,i)'); %The transpose operator ' incorporates conjugate!!!!!
        end
    else
        for i = 1:lind
            Z(:,:,indVec(i)) = (U(:,:,i)*S(:,:,i))*V(:,:,i)';
        end
    end

    for i = 1:length(initD.ind{3})
        Z(:,:,initD.ind{3}(i)) = conj(Z(:,:,initD.ind{2}(i)));
    end

    %either order of ifft works
    for i = 3:p %p:-1:3
        Z = ifft(Z,[],i);
    end
    
    %B update
    contD.B = contD.B+contD.X-Z;
    
    %X Update
    Xold = contD.X;
    if strcmp(iterP.mode,'constrained')
        %X = tProduct(teye(sA)-A, Z-B)+Y; %constrained case
        temp = Z-contD.B;
        contD.X = temp-initD.A(temp)+initD.Y; %constrained case
    elseif strcmp(iterP.mode,'unconstrained')
        %contD.X = initD.A(iterP.rho*(Z-contD.B)+initD.Y)./(iterP.rho+1); %unconstrained
        %temp = (Z-contD.B)+(initD.Y./iterP.rho);
        %contD.X = temp - initD.A(temp)./(iterP.rho+1);
        temp = iterP.rho*(Z-contD.B) + initD.Y;
        %disp(['temp(1) = ',num2str(temp(1))])
        contD.X = (temp)./unconCon;
        %disp(['X(1) = ',num2str(contD.X(1))])
    end
    %toc
    AV.totalTime = AV.totalTime + toc;
    
    iter = iter+1; %update iteration count
    
    %Iteration Vars and Values for output
    %s = ['Iter ', num2str(iter), ': Z ', num2str(tNormF(Z)), ' B ', num2str(tNormF(B))]
    AV.res(iter) = tNormF(initD.A(contD.X)-initD.Y); %tNormF(tProduct(A,X)-Y);
    if strcmp(iterP.mode,'constrained')
        AV.obj(iter) = TNNmap2(contD.X,initD.ind); %set = 1 to save time if don't care about objective
    elseif strcmp(iterP.mode,'unconstrained')
        AV.objT1(iter) = TNNmap2(contD.X,initD.ind);
        AV.obj(iter) = lambda/2*AV.objT1(iter)+AV.res(iter)^2/2; %may need to change to 2-norm...
    end
    
    if isfield(AV,'rGTE')
        AV.rGTE(iter) = tNormF(initD.T-contD.X)/Tn;
    end
    
    temp = tNormF(Xold);
    if temp ~= 0
        AV.convR(iter) = tNormF(contD.X-Xold)/temp;
    end
    if mod(iter,iterP.outF) == 0
        %disp(['Iter ', num2str(iter), ': res ', num2str(res(iter)), ', obj ', num2str(obj(iter)), ', GTE rel ', num2str(rGTE(iter)),', Conv ', num2str(convR(iter))])
        %disp(['Iter ', num2str(iter), ': res ', num2str(res(iter)), ', obj ', num2str(obj(iter)),', Conv ', num2str(convR(iter))])
        if isfield(AV,'rGTE')
            disp(['Iter ', num2str(iter), ': res ', num2str(AV.res(iter)),', obj ', num2str(AV.obj(iter)),', rGTE ', num2str(AV.rGTE(iter)),', Conv ', num2str(AV.convR(iter))]);
        else
            disp(['Iter ', num2str(iter), ': res ', num2str(AV.res(iter)),', obj ', num2str(AV.obj(iter)),', Conv ', num2str(AV.convR(iter))])
        end
    end
    
%     %rho update
%     sN = tNormF(rho*(Z-Zpre));
%     rN = tNormF(X-Z);
%     if rN > mu*sN
%         rho = tinc*rho;
%     elseif sN > mu*rN
%         rho = rho/tdec;
%     end

    %trim vectors if escape condition is met
    if AV.convR(iter-1) < iterP.epsilon && AV.convR(iter-1) > 0
        AV.res = AV.res(1:iter);
        AV.obj = AV.obj(1:iter);
        if isfield(AV,'objT1')
            AV.objT1 = AV.objT1(1:iter);
        end
        if isfield(AV,'rGTE')
            AV.rGTE = AV.rGTE(1:iter);
        end
        AV.convR = AV.convR(1:iter);
        break
    end
end

%In previous code had an extra X update here.

end
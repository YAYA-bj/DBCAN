function [seita, y, A, evs] = DBCAN(X, c, k, r, islocal)
% X: dim*num data matrix, each column is a data point
% c: number of clusters
% k: number of neighbors to determine the initial graph, and the parameter r if r<=0
% r: paremeter, which could be set bo a large enough value. If r<0, then it is determined by algorithm with k
% islocal: 
%           1: only update the similarities of the k neighbor pairs, the neighbor pairs are determined by the distances in the original space 
%           0: update all the similarities
% seita: dim*dim weighted matrix
% y: num*1 cluster indicator vector
% A: num*num learned symmetric similarity matrix
% evs: eigenvalues of learned graph Laplacian in the iteration

NITER = 100;
[dim, num] = size(X);  

if nargin < 6
    islocal = 1;
end;
if nargin < 5
    r = -1;
end;
if nargin < 4
    k = 15;
end;

% Initialize the similarity matrix A
distX = L2_distance_1(X,X);
distsqrtX = sqrt(distX);
sumdistX = sum(distX,2);   
distX = diag(sumdistX) * distsqrtX;
[distX1, idx] = sort(distX,2);
A = zeros(num);
rr = zeros(num,1);
for i = 1:num
    di = distX1(i,2:k+2);
    rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
    id = idx(i,2:k+2);
    A(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);%(LZR)initialize A
end;

% Initialize matrix seita
seita = zeros(dim);
ZJ = zeros(1,num);
W = zeros(dim);
W1 = zeros(dim);
for m =1:dim
    Xm = X(m,:);
    distXm = L2_distance_1(Xm,Xm);   
    dXm = sum(distXm);   
    ZJ = diag(dXm') * distsqrtX .* A;
    W(m,m) =sum(ZJ(:))+eps;
    W1(m,m)=1./W(m,m);    
end;  
for m =1:dim
    seita(m,m) = 1/( W(m,m) * sum(W1(:)) );
end;

if r <= 0
    r = mean(rr);
end;
lambda = r;

% Initialize the similarity matrix F
A0 = (A+A')/2;
D0 = diag(sum(A0));
L0 = D0 - A0;
[F, temp, evs]=eig1(L0, c, 0);

if sum(evs(1:c+1)) < 0.00000000001
    error('The original graph has more than %d connected component', c);
end;

for iter = 1:NITER   
    distf = L2_distance_1(F',F');
    diststx = L2_distance_1(seita*X,seita*X);
    dk = sum(diststx,2);
    distx = diag(dk) * distsqrtX;
    if iter>5
        [temp, idx] = sort(distx,2);
    end;
    A = zeros(num);
    for i=1:num
        if islocal == 1
            idxa0 = idx(i,2:k+1);
        else
            idxa0 = 1:num;
        end;
        dfi = distf(i,idxa0);
        dxi = distx(i,idxa0);
        ad = -(dxi+lambda*dfi)/(2*r);
        A(i,idxa0) = EProjSimplex_new(ad);% (LZR) Update A
    end;

    A = (A+A')/2;
    D = diag(sum(A));
    L = D-A;% (LZR) Update L
    
    for m =1:dim
        Xm = X(m,:);
        distXm = L2_distance_1(Xm,Xm);
        dXm = sum(distXm);
        ZJ = diag(dXm') * distsqrtX .* A;
        W(m,m) =sum(ZJ(:)) + eps;
        W1(m,m)=1./W(m,m); 
    end;
    for m =1:dim
        seita(m,m) = 1/( W(m,m) * sum(W1(:)) );    % (LZR) Update matrix seita
    end;
                      
    F_old = F;
    [F, temp, ev]=eig1(L, c, 0);% (LZR) Update F
    evs(:,iter+1) = ev;
    

    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    if fn1 > 0.00000000001
        lambda = 2*lambda;
    elseif fn2 < 0.00000000001
        lambda = lambda/2;  F = F_old;
    else
        break;
    end;
end;

[clusternum, y]=graphconncomp(sparse(A)); y = y';
if clusternum ~= c
    sprintf('Can not find the correct cluster number: %d', c)
end;


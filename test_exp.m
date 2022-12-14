clc;
clear all;

% sample
load('F:\Ñ§Ï°ÎÄÏ×\¾ÛÀà\CAN code\SWCAN\data\glass.data');
X = glass(:,2:10);
label = glass(:,11);
c = 6;

y = label;
[m,n] = size(X);

for i =1:n
    if max(X(:,i)) == min(X(:,i))
        X(:,i) = 0;
    else
        X(:,i) = 10 * ( X(:,i) - min( X(:,i))) / (max(X(:,i)) - min(X(:,i)));
    end;
end;


NITER = 100;          
for iter = 1:NITER
    k = iter    
    [seita,la, A, evs] = DBCAN(X', c ,k);  
    MIhat = nmi(y', la');
    [Fscore,Accuracy] = Fmeasure(y', la');
end;


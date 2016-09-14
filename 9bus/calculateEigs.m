%% Script takes data and calculates Eigs
%   More precisely, takes low rank updates and runs Arnoldi
%   with Sherman-Morrison-Woodbury calculations to the inverse

numcontigs = 6;
Adata = dlmread('data/matrix');
A = spconvert(Adata);

for i = 1:numcontigs
    udata = dlmread(strcat('data/u',int2str(i)));
    u = spconvert(udata);
    vdata = dlmread(strcat('data/v',int2str(i)));
    v = spconvert(vdata);
    n = length(A);
    [Q,H] = InvArnoldi(A,ones(n,1),n,0,u,v,sprank(u));
    display(2*pi./imag(eig(H)));
end
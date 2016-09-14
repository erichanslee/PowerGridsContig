function out = ShermanIteration(l,u,p,x,y,b,rank)
%% ~~SPARSE VERSION~~


%   ShermanIteration() solves the equation
%   (A^-1 + A^-1X(I + Y^TA^-1X)^-1U^tA^-1)b = out 
%   Where the LU factorization of A^-1 is already given
%   PA^-1 = LU
%   and rank is the size of x and y


temp = l\p*b;
invAb = u\temp;

invAx = x;
for i = 1:rank
%     size(l)
%     size(p)
%     size(x)
    temp = l\p*x(:,i);
    invAx(:,i) = u\temp;
end
W = eye(rank) + y'*invAx;
invWy = W\y';

out = invAb - invAx*invWy*invAb;

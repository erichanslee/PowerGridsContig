function [Q,H] = InvArnoldi(A,q1,m,shift,x,y,rank)

%   A = original matrix
%   q1 = initial guess
%   m = iters
%   shift = shift for inverse iteration
%   x,y = low rank update x*y'
%   rank = rank of x,y


%   [Q,H] = ARNOLDI(A,q1,M) carries out M iterations of inv Arnoldi
%   iteration with N-by-N matrix A+xy^T and starting vector q1
%   (which need not have unit 2-norm).  For M < N it produces
%   an N-by-(M+1) matrix Q with orthonormal columns and an
%   (M+1)-by-M upper Hessenberg matrix H such that
%   Ainv*Q(:,1:M) = Q(:,1:M)*H(1:M,1:M) + H(M+1,M)*Q(:,M+1)*E_M',
%   where E_M is the M'th column of the M-by-M identity matrix.
%   and where Ainv = (A+xy^T - shift*I)^-1

%   Note that we usually refer to a low rank update as u*v', but
%   as we are using the variable 'u' for lu factorization, we instead
%   use x*y' for the low rank update

eps = 10^(-8);
n = length(A);
q1 = q1/norm(q1);
Q = zeros(n,m); Q(:,1) = q1;
H = zeros(min(m+1,m),n);

[l,u,p] = lu(A-shift*eye(n));

for k=1:m
    b = Q(:,k);
    z = ShermanIteration(l,u,p,x,y,b,rank);
    for i=1:k
        H(i,k) = Q(:,i)'*z;
        z = z - H(i,k)*Q(:,i);
    end
    if k < n
       H(k+1,k) = norm(z);
       if abs(H(k+1,k)) < eps, return, end
       Q(:,k+1) = z/H(k+1,k);
   end
end
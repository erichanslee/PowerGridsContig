L1 = -.1 + .2*1i;
L2 = -.2 + .4*1i;

t = 0:.1:10; 
x1 = 2*exp(L1*t)'+ 4*exp(L2*t)'; 
x2 = exp(L1*t)'+4*exp(L2*t)';
x3 = 2*exp(L1*t)'+ 2*exp(L2*t)';
z = iddata([x1 x2 x3],zeros(length(t),1),.1);
m = n4sid(z,3,'Form','modal','DisturbanceModel','none');    

A = m.C*m.A*inv(m.C);
[x,d] = eig(A)
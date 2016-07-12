% Starting point: Discrete time system with two frequencies
L1 = -.1 + .2*1i;
L2 = -.2 + .4*1i;
xi1 = exp(L1/10);
xi2 = exp(L2/10);

% Scramble A a little to make things interesting.
[Q,~] = qr(randn(2));
A = Q*diag([xi1, xi2])*Q';

% Generate a signal
xs = zeros(2,20);
xs(:,1) = randn(2,1);
for j = 2:20
  xs(:,j) = A*xs(:,j-1);
end
% xs = real(xs);

% Now try to identify
u = zeros(1,20);
z = iddata(xs', u');
m = n4sid(z,4,'Form','modal','DisturbanceModel','none');
[mx,md] = eigs(m.A);

predvecs = normalizematrix(m.C*mx);
fprintf('Checking that frequencies match up\n');
display(log(diag(md)));
display([L1 L2]);
fprintf('Checking that eigenvalues match up\n');

disp(Q'*conj(predvecs));
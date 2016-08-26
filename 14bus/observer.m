function observerdata = observer(contignum, matrixnum, noise, PMUidx)

maxfreq = .5;
minfreq = .05;
tstep = .05;
load metadata.mat;

%%  Randomly Place PMUs and Offset data
win = place_PMU(contignum, PMUidx);
rangebus = (differential + numlines + 1):(differential + numlines + numlines);

%   Load and offset data
filename = ['data/busdata_' num2str(contignum) '.mat'];
load(filename);
offset = 50;
data = data(offset:end, win - (differential + numlines));
data = data';
% Calculate Eigenvalue and Eigenvector Predictions from State Matrix
% from the reduced state matrix
I = eye(differential);
E = zeros(algebraic + differential);
E(1:differential,1:differential) = I;
A = full(matrix_read(sprintf('data/matrixfull%d', contignum)));
[vi,di] = eig(A,E); 

% Organize data from State Matrix properly
[linvecsEntire, linvals] = filter_eigpairs(minfreq, maxfreq, diag(di), vi);
linvecs = linvecsEntire(rangebus,:); 
linearvecs = linvecsEntire;

M = form_ode(linvecs, linvals);
x0 = 1.2*ones(numlines,1);
K = -4*M;
[len,~] = size(data);
observerdata = zeros(size(data));
observerdata(:,1) = x0;
for i = 1:len
	xobs_old = observerdata(:,i);
	xdata_old = data(:,i);
	Prop = exp(M*i*tstep);
	xnew = Prop*xobs_old + Prop^4*(xdata_old - xobs_old);
	observerdata(:,i+1) = xnew;
end
hold on;
plot(observerdata(1,:));
plot(data(1,:));

end
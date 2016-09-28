% Checks correctness of id_contig.m

clear all;
load metadata.mat
maxfreq = .5;
minfreq = .05;
contignum = 2;
n = differential + algebraic;
method = 'OrthReg';
win = (differential + numlines + 1):(differential + numlines + numlines);
noise = 0;

%% Test Case: Power System (Doesn't work well)
I = eye(differential);
E = zeros(algebraic + differential);
E(1:differential,1:differential) = I;
A = full(matrix_read(sprintf('data/matrixfull%d', contignum)));
[v2,d2] = eig(A,E); 
[v2_subset, d2_subset] = filter_eigpairs(minfreq, maxfreq, diag(d2), v2);


%%   Load and offset data
filename = ['data/busdata_' num2str(contignum) '.mat'];
load(filename);
offset = 50;
data = data(offset:end, win - (differential + numlines));
[evals, evecs]  = run_n4sid(data, noise, timestep, numlines, maxfreq, minfreq);
evecs = normalizematrix(evecs);
[res, evecs_fitted] = id_contig(A, E, method, evals, evecs, win);




disp('Dot Product of Fitted and Real eigenvectors:');
disp(abs(evecs_fitted'*normalizematrix(v2_subset)));

plot_eigvecs(normalizematrix(v2_subset),evecs_fitted);

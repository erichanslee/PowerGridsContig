% isolatecase is an auxillary function for debugging purposes that checks
% runs a time-domain simulation for contingency <contignum> and then
% calculates residuals for the linearized system <matrixnum>

% That is, it checks data from contingency identification routine on a
% cases-by-case basis for debugging purposes.

% NOTE: code contains quite a bit of similarity with other files but
% is separate due to large amount of output. 

% ~~~~~~~~~INPUTS~~~~~~~~~ %

% method = method type number one would like to use
% contignum = contingency to simulate
% matrixnum = linearized system to use for residual calculation
% noise = percentage of max amplitude to add as gaussian noise

% ~~~~~~~~~OUTPUTS~~~~~~~~~ %

% linearvecs = eigenvectors from the linearized system Jacobian
% empvecsfull = eigenvectors from fitting
% empresidual = residual from fittings

function [linearvecs, empvecsfull, empresidual] = isolatecase(method, contignum, matrixnum, noise, window)

maxfreq = .5;
minfreq = .05;
load('metadata.mat')

%% Basic Pre-Run Checks
if(method > 5 || method < 1)
    error('Problems with parameter "Method". Please an integer in [1,3]')
end

if(noise > 1 || noise < 0)
    error('Problems with parameter "Noise". Please enter in a real number in the range of [0,1]')
end

if(window > 1 || window < 0)
    error('Problems with parameter "Window". Please enter in a real number in the range of [0,1]')
end

if(contignum > numcontigs)
    error('Contigency Number not found! Please enter a smaller integer')
end

%%  Randomly Place PMUs and Offset data
PMU = place_PMU(contignum, window);
rangebus = (differential + numlines + 1):(differential + numlines + numlines);

%   Load and offset data
filename = ['data/busdata_' num2str(contignum) '.mat'];
load(filename);
offset = 50;
data = data(offset:end, PMU - (differential + numlines));

[empvals, empvecs]  = run_n4sid(data, noise, timestep, numlines, maxfreq, minfreq);
empvecs = normalizematrix(empvecs);
data_dump = zeros(1,numcontigs);

%% Calculate Eigenvalue and Eigenvector Predictions from State Matrix
% from the reduced state matrix
I = eye(differential);
E = zeros(algebraic + differential);
E(1:differential,1:differential) = I;
A = full(matrix_read(sprintf('data/matrixfull%d', contignum)));
[vi,di] = eig(A,E); %solve generalized eigenvalue problem

%% Organize data from State Matrix properly
linvals = (diag(di));
rangepred = find(abs(imag(linvals)/2/pi) > minfreq & abs(imag(linvals)/2/pi) < maxfreq);
linvals = linvals(rangepred);
[~, idx1] = sort(abs(imag(linvals)));

%cut and sort eigenpairs of state matrix
linvals = linvals(idx1);
linvecs = vi(rangebus,rangepred); %where the important eigenvectors are
linvecs = linvecs(:,idx1);
linvecsEntire = vi(:,rangepred);
linvecsEntire = linvecsEntire(:,idx1);
linearvecs = linvecsEntire;

% form DAE matrix E
I = eye(differential);
E = zeros(algebraic + differential);
E(1:differential,1:differential) = I;
A = full(matrix_read(sprintf('data/matrixfull%d', matrixnum)));
format long

[empresidual, empvecsfull] = id_contig(A, E, method, empvals, empvecs, PMU);

end

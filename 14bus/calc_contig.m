%% calc_contig loads 

% ~~~~~~~~~INPUTS~~~~~~~~~ %
% method = method type number one would like to use
% data = voltage readings from PMUs
% PMU = indices of PMU placements
% rangerest = indices of everything else
% noise = boolean value indicating presence of noise

% ~~~~~~~~~OUTPUTS~~~~~~~~~ %
% predcontig = the cotingency the chosen method predicts
% confidence = the confidence levels for correctly identified contigs

function [predcontig, confidence, list] = calc_contig(method, data, PMU, noise)


maxfreq = .5;
minfreq = .05;
load('metadata.mat');
list = cell(1,numcontigs);

%use n4sid
[empvals, empvecs]  = run_n4sid(data, noise, timestep, numlines, maxfreq, minfreq);
empvecs = normalizematrix(empvecs);
data_dump = zeros(1,numcontigs);
clc;

% form DAE matrix E
I = eye(differential);
E = zeros(algebraic + differential);
E(1:differential,1:differential) = I;

for k = 1:numcontigs    
    % Read in matrix
    A = full(matrix_read(sprintf('data/matrixfull%d', k)));
    format long
    
    %% Calculate Backward Error
    [fittedres, fittedvecs] = id_contig(A, E, method, empvals, empvecs, PMU);
    data_dump(k) = norm(fittedres);
    list{k} = fittedvecs;
end

% calculate contingency and confidence measure
[val1, idx1] = min(data_dump);
data_dump(idx1) = [];
[val2, ~] = min(data_dump);
predcontig = idx1;
confidence = abs(val2 - val1)/abs(val1);

end

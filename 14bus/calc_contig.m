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

function [predcontig, confidence] = calc_contig(method, data, PMU, rangerest, noise)

maxfreq = .5;
minfreq = .05;
load('metadata.mat');

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
    Ifull = eye(differential + algebraic);
    order = [PMU, rangerest];
    P = Ifull(order,:);
    out = zeros(length(empvals),1);
    for j = 1:length(empvals)

        % form variables to pass into calc_residual
        lambda = empvals(j);
        Ashift = A-lambda*E;
        xfull = zeros(differential + algebraic,1);
        x1 = empvecs(:,j);
        res = calc_residual(method, Ashift, x1, PMU, rangerest, xfull, P);
        out(j) = norm(res);
    end
    data_dump(k) = mean(out);
end

% calculate contingency and confidence measure
[val1, idx1] = min(data_dump);
data_dump(idx1) = [];
[val2, ~] = min(data_dump);
predcontig = idx1;
confidence = abs(val2 - val1)/abs(val1);

end

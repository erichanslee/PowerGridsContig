%% NOTE: tesbed is not used anymore to generate confusion matrix

% ~~~~~~~~~INPUTS~~~~~~~~~ %

% method = method type number one would like to use
% noise = percentage of max amplitude to add as gaussian noise
% percentage = percentage of PMUs visible

% ~~~~~~~~~OUTPUTS~~~~~~~~~ %

% predcontig = the cotingency the chosen method predicts
% actualcontig = the contingency that was actually simulated
% confidence = the confidence levels for correctly identified contigs

function [predcontig, actualcontig, confidence] = testbed(method, noise, percentage)

load('metadata.mat')

%% Basic Pre-Run Checks
if(method > 5 || method < 1)
    error('Problems with parameter "Noise". Please an integer in [1,3]')
end

if(noise > 1 || noise < 0)
    error('Problems with parameter "Noise". Please enter in a real number in the range of [0,1]')
end

if(contignum > numcontigs)
    error('Contigency Number not found!')
end

%   randomly pick contig
contignum = ceil(numcontigs*rand);
actualcontig = contignum;


%   Load contig data
filename = ['data/busdata_' num2str(contignum) '.mat'];
load(filename);
offset = 50;

%%  Randomly Place PMUs and Offset data
PMU = gen_PMUidx(percentage, numbuses)
PMU = place_PMU(contignum, PMU);
data = data(offset:end, PMU - (differential + numlines));

% predict contingency
[predcontig, confidence, ~, ~] = calc_contig(method, data, PMU, noise);

end

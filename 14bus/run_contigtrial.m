
function [scores, ranking, vecs, res] = run_contigtrial(fitting_method, score_method, contignum, PMUidx)

load metadata.mat
noise = 0;

%   Load data
filename = ['data/busdata_' num2str(contignum) '.mat'];
load(filename);
offset = 50;
win = place_PMU(contignum, PMUidx);
data = data(offset:end, win - (differential + numlines));

% calc fit data
[vecs, res] = calc_contig(fitting_method, data, win, noise);

% predict contingency
[ranking, scores] = calc_scores(score_method, res);

end

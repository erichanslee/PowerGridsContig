%   place_PMU simulates virtual placement of PMUs by generating
%   a set of indices called "PMU" in which voltages can be read

%   ~~~TODO: Make more realistic taking into account network topology~~~

function [PMU, rangerest] = place_PMU(rangebus,window)
load('metadata.mat');

PMUnum = round(window*length(rangebus));
PMU = randsample(rangebus,PMUnum);
PMU = sort(PMU);
outrange = setdiff(rangebus, PMU);
rangerest = [1:(differential + numlines), (differential + numlines + numlines + 1): (differential + algebraic), outrange];

end
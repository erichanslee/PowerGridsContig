%   place_PMU simulates virtual placement of PMUs by generating
%   a set of indices called "PMU" in which voltages can be read

%   ~~~TODO: Make more realistic taking into account network topology~~~

function [PMU] = place_PMU(rangebus,window)

PMUnum = round(window*length(rangebus));
PMU = randsample(rangebus,PMUnum);
PMU = sort(PMU);


end
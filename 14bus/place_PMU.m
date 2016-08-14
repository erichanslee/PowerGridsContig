%   place_PMU simulates virtual placement of PMUs by generating
%   a set of indices called "PMU" in which voltages can be read

% ~~~~~~~~~INPUTS~~~~~~~~~ %

% rangebus
% contignum = contig number for knowledge of line faliure
% window = percentage PMU placement

% ~~~~~~~~~OUTPUTS~~~~~~~~~ %

% PMU = indices for PMU placement

function PMU = place_PMU(rangebus, contignum, window)

len = length(rangebus);
PMUnum = round(window*len);
PMUidx = randsample(1:len,PMUnum);
run(sprintf('contig%d.m',contignum));
Lines = Line.con(:,1:2);
Lines(contignum,:) = [];

Agg_Neighbor = [];
for i = 1:length(PMUidx)
	busnum = PMUidx(i);
	[rowidx, colidx] = find(Lines == busnum);
	neighbors = Lines(rowidx,:);
	neighbors = neighbors(:)';
	neighbors = unique(neighbors);
	Agg_Neighbor = [Agg_Neighbor, neighbors]; 
end

PMUidx = sort(unique(Agg_Neighbor));
PMU = rangebus(PMUidx);

end
	
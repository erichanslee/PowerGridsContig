% Run time-domain simulations for each contingency in the 14 bus system
%
function psat_runtrials()
  load metadata.mat
  psat_runtrial(basefilename, 'data/sim14');
  for k = 1:numcontigs
    psat_runtrial(sprintf('contig%d', k), ...
                  sprintf('data/sim14_%d.mat', k))
  end
end

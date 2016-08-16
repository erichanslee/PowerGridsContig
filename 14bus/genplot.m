
% ~~~~~~~~~INPUTS~~~~~~~~~ %

% method = method type number one would like to use
% numtrials = number of trials desired
% noise = percentage of max amplitude to add as gaussian noise
% window = percentage of PMUs visible

% ~~~~~~~~~OUTPUTS~~~~~~~~~ %

% A png file with the desired image
% out = The Confusion Matrix
% confidence = array of confidence measures from each
%               trial instance


% TODO: fix calc_contig, this as right now M is printing out the norms
% of the eigenvectors instead.
function [M] = genplot(method, noise, PMUidx)

  load metadata.mat

  %   temp variables
  M = zeros(numcontigs);

  for i = 1:numcontigs
      %   Load data
      filename = ['data/busdata_' num2str(i) '.mat'];
      load(filename);
      offset = 50;
      win = place_PMU(i, PMUidx);
      % Offset and Isolate data
      data = data(offset:end, win - (differential + numlines));

      % predict contingency
      [predcontig, confidence, ~, list] = calc_contig(method, data, win, noise);

      % Populate Confusion Matrix M 
      for j = 1:numcontigs
        M(i,j) = norm(list{j});
      end
  end


end

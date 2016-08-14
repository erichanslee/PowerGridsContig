
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
function [M] = genplot(method, noise, PMU)

  load metadata.mat

  %   temp variables
  M = zeros(numcontigs);

  for i = 1:numcontigs
      %   Load data
      filename = ['data/busdata_' num2str(i) '.mat'];
      load(filename);
      offset = 50;

      %%  Randomly Place PMUs and Offset data
      data = data(offset:end, PMU - (differential + numlines));

      % predict contingency
      [predcontig, confidence, list] = calc_contig(method, data, PMU, noise);

      % Populate Confusion Matrix M 
      for j = 1:numcontigs
        M(i,j) = norm(list{j});
      end
  end


end

function out = genplot(numtrials, noise, window)
load metadata.mat

%%  Case with full PMU
M = zeros(numcontigs);
for n = 1:numtrials
    [i,j] = testbed(noise, window);
    M(i,j) = M(i,j) + 1;
end


out = M;

figure 
imagesc(M);
colormap(bone);
set(gca,'XTick',[]); % Remove the ticks in the x axis!
set(gca,'YTick',[]); % Remove the ticks in the y axis
set(gca,'Position',[0 0 1 1]); % Make the axes occupy the hole figure
saveas(gcf,'results14bus','png');


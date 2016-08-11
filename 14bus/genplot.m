function [out,confidence] = genplot(numtrials, noise, window)
load metadata.mat

%%  Case with full PMU
M = zeros(numcontigs);
C = zeros(1,numtrials);
counter = 1;
for n = 1:numtrials
    [i,j,con] = testbed3(noise, window);
    M(i,j) = M(i,j) + 1;
    
    if(i == j)
        C(counter) = con;
        counter = counter + 1;
    end
end


out = M;
confidence = C;
figure 
imagesc(M);
colormap(bone);
set(gca,'XTick',[]); % Remove the ticks in the x axis!
set(gca,'YTick',[]); % Remove the ticks in the y axis
set(gca,'Position',[0 0 1 1]); % Make the axes occupy the hole figure
saveas(gcf,'results14bus','png');


function out = genplot()
load metadata.mat

%%  Case with full PMU
M = zeros(numcontigs);
noise = 0;
window = 1;
for i = 1:numcontigs
    for j = 1:numcontigs
        out = testbed(i,j,0,window);
        M(i,j) = mean(out);
    end
end

out = M;

figure 
imagesc(M);
colormap(bone);
set(gca,'XTick',[]); % Remove the ticks in the x axis!
set(gca,'YTick',[]); % Remove the ticks in the y axis
set(gca,'Position',[0 0 1 1]); % Make the axes occupy the hole figure
saveas(gcf,'results14bus_100percent','png');

clear all;
load metadata.mat
%%  Case with most PMUs
M = zeros(numcontigs);
noise = 0;
window = .9;
for i = 1:numcontigs
    for j = 1:numcontigs
        out = testbed(i,j,0,window);
        M(i,j) = mean(out);
    end
end

out = M;

figure 
imagesc(M);
colormap(bone);
set(gca,'XTick',[]); % Remove the ticks in the x axis!
set(gca,'YTick',[]); % Remove the ticks in the y axis
set(gca,'Position',[0 0 1 1]); % Make the axes occupy the hole figure
saveas(gcf,'results14bus_90percent','png');

clear all;
load metadata.mat
%%  Case with half PMUs
M = zeros(numcontigs);
noise = 0;
window = .5;
for i = 1:numcontigs
    for j = 1:numcontigs
        out = testbed(i,j,0,window);
        M(i,j) = mean(out);
    end
end

out = M;

figure 
imagesc(M);
colormap(bone);
set(gca,'XTick',[]); % Remove the ticks in the x axis!
set(gca,'YTick',[]); % Remove the ticks in the y axis
set(gca,'Position',[0 0 1 1]); % Make the axes occupy the hole figure
saveas(gcf,'results14bus_50percent','png');

clear all;


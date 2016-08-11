%% RUNME (note: everything will be slow due to simulation times)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~README~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%   SET UP: given (M-lambda*E) and a subset of a potential eigenvector x1,
%   how best fit x_\2 to minimize || (M - lambda*E)[x1;x2] ||?

%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%       Description of (Contingency Identification) Methods:
%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%   Method 1 = least squares fit subject to || x1 || = 1

%   Method 2 = least squares fit subject to || x1 || = 1 with a
%   normalization after

%   Method 3 = least squares fit subject to || [alpha*x1; x2] || = 1

%   Method 4 = least squares fit subject to || [alpha*x1; x2] || = 1
%   and then an reverse-normalization via scaling by 1/alpha to bring
%   || x1 || = 1

%   Method 5 = dummy method simply calculating "true" eigenvectors of
%   linearized system

%       ~~~~~~~~~~~~~~~~~~~~~~~~
%       Description of Functions
%       ~~~~~~~~~~~~~~~~~~~~~~~~

%   testbed.m: randomly picks contingency and then runs contingency
%   identification

%   isolatecontingency.m: picks contingency and linearized system subject
%   to your choice, runs contingency identification

%   genplot.m: Generates confusion matrix using n trials of testbed.m

%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%       Other Equally Important Things
%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%   **There are a number of sections in this runme, so you can run a section
%   at a time if you'd like or just run everything at once (although that
%   will take some time)

%   **Please don't set noise to anything but 0 right now it's not correctly
%   implemented :) 

%   **The numbering of the Methods is the same as the report with the
%   exception that Method 4 is the modified Method 3 and Method 5 is just a
%   dummy. 

%   **Figures/plots will be saved in the "Figures" folder.

%   ****RUN genscript14 FIRST BEFORE DOING ANYTHING ELSE AS IT GENERATES NECESSARY FILES****

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~END README~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Just to make things clearer
Method1 = 1; Method2 = 2; Method3 = 3; Method4 = 4; Method5 = 5;

%% RUN genscript14 FIRST... SUPER IMPORTANT
genscript14;


%% Generates Confusion Matrix of your choice, and confidence plot
load metadata.mat;
%   Confusion matrix will be saved in confusion.png
%   Confidence plot will be saved in conplot.png

window = 1;
noise = 0;
numtrials = 15;

%   Run genplot (genplot itself prints the confusion matrix)
%   Compares Methods 1 and 3 by default (same plot as report)
[confusion1, confidence1] = genplot(Method1,numtrials,noise,window);
[confusion3, confidence3] = genplot(Method3,numtrials,noise,window);

%   Plot confusion matrices and confidence graph
f = figure();
hold on;
plot(confidence1,'-*','LineWidth',1.5); plot(confidence3,'-rs','LineWidth',1.5);
legend('Method1', 'Method3');
hold off
saveas(f,'figures/conplot.png');
close

figure
imagesc((-confusion1));
colormap(bone);
set(gca,'XTick',[]); % Remove the ticks in the x axis!
set(gca,'YTick',[]); % Remove the ticks in the y axis
set(gca,'Position',[0 0 1 1]); % Make the axes occupy the hole figure
saveas(gcf,'figures/confusion1.png','png');
close

figure
imagesc((-confusion3));
colormap(bone);
set(gca,'XTick',[]); % Remove the ticks in the x axis!
set(gca,'YTick',[]); % Remove the ticks in the y axis
set(gca,'Position',[0 0 1 1]); % Make the axes occupy the hole figure
saveas(gcf,'figures/confusion3.png','png');
close

%%  Generates Graph from "Sanity Checks" Section
%   Tests to see the difference in fits between methods. Plots the
%   size of x1 relative to the entire vector i.e.
load metadata.mat;

%   Parameters: Modify as needed. 
window = 1;
noise = 0;

% (prop stands for Proportion)
prop1 = [];
prop2 = [];
prop3 = [];

% If there is a discrepancy between number of empirical pairs and
% theoretical pairs, isolatecase will return a vector of zeros
for i = 1:numcontigs
    prop = isolatecase(Method1,i,i,noise,window);
    prop1 = [prop1 prop];
    prop = isolatecase(Method3,i,i,noise,window);
    prop2 = [prop2 prop];
    prop = isolatecase(Method5,i,i,noise,window);
    prop3 = [prop3 prop];
end

% Plot and save 
f = figure();
hold on;
plot(prop1,'-*','LineWidth',1.5); hold on; plot(prop2,'-rs','LineWidth',1.5);
plot(prop3,'-go','LineWidth',1.5); legend('Method 1/2', 'Method 3', 'Linearized System')
hold off;
saveas(f,'figures/fitplot.png');
close

%% Basic Skeletons of things you might want to run

% contignum = i;
% matrixnum = j;
% methodnum = k
% noise = 0;
% window = w;
% numtrials = n;
% 
% 
% % Using isolatecase to test different cases
% prop = isolatecase(methodnum,contignum,matrixnum,noise,window);
% 
% % Running instances of testbed
% [predcontig,actualcontig,confidence] = testbed(methodnum,noise,window);
% 
% % Complete Comparison of All Methods
% [confusion1, confidence1] = genplot(Method1,numtrials,noise,window);
% [confusion2, confidence2] = genplot(Method2,numtrials,noise,window);
% [confusion3, confidence3] = genplot(Method3,numtrials,noise,window);
% [confusion4, confidence4] = genplot(Method4,numtrials,noise,window);

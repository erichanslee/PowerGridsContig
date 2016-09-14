%% RUNME

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


%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%        Important Things
%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%   **Please don't set noise to anything but 0 right now it's not correctly
%   implemented :) 

%   **The numbering of the Methods is the same as the report with the
%   exception that Method 4 is the modified Method 3.

%   **Figures/plots will be saved in the "Figures" folder.

%   ****RUN genscript14 FIRST BEFORE DOING ANYTHING ELSE AS IT GENERATES NECESSARY FILES****

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~END README~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Just to make things clearer
Method1 = 1; Method2 = 2; Method3 = 3; Method4 = 4; 
load metadata.mat;

% PMUs at all buses
percentage = 1;
noise = 0;
idx = gen_PMUidx(percentage, numbuses);


% Contig 3 is a case where misidentification occurs (Contig 4 identified instead)
% NOTE: PMU data is on rows 64:77 of linearvecs/empvecs/empresidual

rangebus = (differential + numlines + 1):(differential + numlines + numlines);

% The output from the correctly-fitted data (that is "worse" by the averaging metric compared to
% the false positive)
method = Method3;
contignum = 3; % Contingency to Simulate
matrixnum = 3; % Jacobian in which to fit eigenvectors
[linearvecs_true, empvecs_true, empresidual_true] = isolatecase(method, contignum, matrixnum, noise, idx);

% The false positive 
method = Method3;
contignum = 3; % Contingency to Simulate
matrixnum = 4; % Jacobian in which to fit eigenvectors
[linearvecs_false, empvecs_false, empresidual_false] = isolatecase(method, contignum, matrixnum, noise, idx);
clc;

% For example, compare portions of false positive eigenvectors seen by buses: 
pmode1 = 'entire';
pmode2 = 'PMU';
j = 2;
printvecs(linearvecs_false, empvecs_false, j, pmode1);


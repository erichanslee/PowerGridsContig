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
idx = gen_PMUidx(percentage, numbuses);


% Contig 3 is a case where misidentification occurs (Contig 4 identified instead)
% NOTE: PMU data is on rows 64:77 of linearvecs/empvecs/empresidual
method = Method3;
contignum = 3;
matrixnum = 3;
noise = 0;
[linearvecs_true, empvecs_true, empresidual_true] = isolatecase(method, contignum, matrixnum, noise, idx);
 
method = Method3;
contignum = 3;
matrixnum = 4;
noise = 0;
[linearvecs_false, empvecs_false, empresidual_false] = isolatecase(method, contignum, matrixnum, noise, idx);
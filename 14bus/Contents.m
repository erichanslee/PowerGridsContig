% Analyze contingencies in 14-bus test case by analyzing empirical modes.
%
% Data files:
%   bus14.m      - PSAT file to specify 14 bus system
%   contig*.m    - PSAT files for the contingency systems
%   metadata.mat - Specify params for 14 bus case
%   results.dat  - Results??  Not documented anywhere.  Just a 10-by-10
%
% Static analysis:
%   genscript   - Generate data/ submatrix with power flow Jacobians
%   genscript14 - Call genscript with inputs for the 14-bus system
%
% Helper functions:
%   normalizematrix - Scale matrix columns to unit norm
%   matrix_write    - Write a matrix to file (dense or sparse)
%   matrix_read     - Read a matrix from file (dense or sparse)

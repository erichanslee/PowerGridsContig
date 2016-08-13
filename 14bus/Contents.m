% Analyze contingencies in 14-bus test case by analyzing empirical modes.
%
% Data files:
%   bus14.m      - PSAT file to specify 14 bus system
%   contig*.m    - PSAT files for the contingency systems
%   metadata.mat - Specify params for 14 bus case
%
% PSAT analysis:
%   genscript      - Generate data/ submatrix with power flow Jacobians
%   gendata        - Save simulation data to busdata*.mat
%   genscript14    - Call genscript with inputs for the 14-bus system
%   psat_runtrial  - Run PSAT sim and save results to sim14*.mat file
%   psat_runtrials - Run PSAT sims for every contingency
%	testbed		   	- Run contingency identification routine
%	isolatecase		- Run contingency identification routine for fixed linearized system	
%
% Helper functions:
%   normalizematrix - Scale matrix columns to unit norm
%   matrix_write    - Write a matrix to file (dense or sparse)
%   matrix_read     - Read a matrix from file (dense or sparse)
%   run_n4sid       - Runs n4sid
%   place_PMU       - "Places" PMUs on data
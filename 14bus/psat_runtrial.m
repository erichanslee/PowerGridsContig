% psat_runtrial(basename, matfname)
%
% Run a time domain simulation.
%
% Inputs:
%   psatname: Name of PSAT input file
%   matfname: Output mat file name
%
function psat_runtrial(psatname, matfname)
  initpsat;
  Settings.freq = 60;     % System frequency is 60 Hz
  Settings.fixt = 1;      % Simulate one second
  Settings.tstep = 0.05;  % Time step is 50 ms
  runpsat(psatname, 'data');
  runpsat('td');
  save(matfname);
end

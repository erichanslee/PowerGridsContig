genscript9;
initpsat;
Settings.freq = 60;
runpsat('CONTIG_ONE','data');
runpsat('td');
for i = 1:numbuses
string = strcat('simulation/sim9bus' , num2str(i), '.txt');
fid = fopen(string, 'w');
fprintf(fid, '%f\n', Varout.vars(10:end,33+i));
end
%A = dlmread('matrix1'); A = spconvert(A); A = full(A);

data = Varout.vars(10:end,34:42);
[length,num] = size(data);
z = iddata(data,zeros(length,1));
m = n4sid(z,9,'Form','modal','DisturbanceModel','none');
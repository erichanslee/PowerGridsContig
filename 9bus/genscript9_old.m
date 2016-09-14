%% script takes in a test system and then outputs the contigency data
% the format being:
%   -The contingencies, separated into low rank changes and put into 
%   csv filies 'ui' and 'vi' respectively, with 'i' being the 
%   contigency number
%
%   -'ui' and 'vi' are calculated via an SVD on the difference between
%   the original and modified (contig) Jacobian, with the singular values
%   multiplied into 'ui'

%   the filename containing the original system information is 'bus9.m'
%   and is the IEEE 9-bus system, with 6 total possible contingencies for
%   line failure. 
basefilename = 'bus9';
filename = strcat(basefilename,'_');
numcontigs = 6;

%   Generates Scratch Space Files (not really important in terms of 
%   the end product)
for i = 1:numcontigs
     temp = strcat(filename,num2str(i));
     fid = fopen(strcat(temp,'.m'),'w');
     fid = fclose(fid);
     copyfile(strcat(basefilename,'.m'),strcat(temp,'.m'));
     cssm(strcat(temp,'.m'),i+1);
end

%   Use PSAT to calculate original Jacobian
initpsat
runpsat(basefilename,'data');
runpsat('pf');
Aoriginal = DAE.Fx - DAE.Fy*(DAE.Gy\DAE.Gx)- 1e-6*speye(DAE.n);
GY = DAE.Gy; GX = DAE.Gx;
Aoriginal = full(Aoriginal);

%   For each contigency i, form the state matrix As and store
%   inside the csv named 'matrixi'. Then calculate the svd of the
%   difference matrix and store as 'ui' and 'vi'. Matrices stored in
%   non-sparse format
filename = strcat('contig_',filename);
for i = 1:numcontigs
    temp = strcat(filename,num2str(i));
    runpsat(temp,'data');
    runpsat('pf');
    As = DAE.Fx - DAE.Fy*(DAE.Gy\DAE.Gx) - 1e-6*speye(DAE.n); 
%   Print ranks and stuff    
%   rank(full(GY - DAE.Gy))
%   rank(full(GX - DAE.Gx)) 
    diff = Aoriginal - full(As);
    [u,d,v] = svd(diff);
    display(d)
    dlmwrite(strcat('u',num2str(i)), u*d);
    dlmwrite(strcat('v',num2str(i)), v);
    dlmwrite(strcat('matrix',num2str(i)),full(As));
    
    %   Print ranks and stuff    
    %   rank(full(GY - DAE.Gy))
    %   rank(full(GX - DAE.Gx)) 
    %   rank(diff)
end

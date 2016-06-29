%% script takes in a test system and then outputs the contigency data 
%   ~~~SPARSE VERSION~~~
%   the format being:
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

clear all;
load('metadata.mat');
filename = strcat(basefilename,'_');


%%   Generates Scratch Space Files (not really important in terms of 
%   the end product)
for i = 1:numcontigs
     temp = strcat(filename,num2str(i));
     fid = fopen(strcat(temp,'.m'),'w');
     fid = fclose(fid);
     copyfile(strcat(basefilename,'.m'),strcat(temp,'.m'));
     cssm(strcat(temp,'.m'),i+1);
end

%%   Use PSAT to calculate original Jacobian
initpsat
runpsat(basefilename,'data');
runpsat('pf');
Aoriginal = DAE.Fx - DAE.Fy*(DAE.Gy\DAE.Gx)- 1e-6*speye(DAE.n);
Aoriginal(abs(Aoriginal) < 10^-10) = 0;
[col,row,val] = find(Aoriginal);
data_dump = [col,row,val];
dlmwrite('data/matrix',data_dump,'precision',16);

%%   For each contigency i, form the state matrix As and store
%   inside the csv named 'matrixi'. Then calculate the svd of the
%   difference matrix and store as 'ui' and 'vi'. Matrices stored in
%   sparse format
filename = strcat('contig_',filename);
format long;
for i = 1:numcontigs
    temp = strcat(filename,num2str(i));
    runpsat(temp,'data');
    runpsat('pf');
    As = DAE.Fx - DAE.Fy*(DAE.Gy\DAE.Gx) - 1e-6*speye(DAE.n); 
    diff = Aoriginal - As;
    
    %   Write out As in a sparse format
    [col,row,val] = find(As);
    data_dump = [col,row,val];
    dlmwrite(strcat('data/matrix',num2str(i)),data_dump,'precision',16);
    [u,d,v] = svds(diff);
    
    %   Write out ui = u*d in a sparse format
    %   making sure it's the correct size by
    %   padding with extra zeros
    temp = u*d;
    [n,m] = size(temp);
    temp(abs(temp) < 10^-10) = 0;
    temp(n,m) = 10^-15;
    [col,row,val] = find(temp);
    data_dump = [col,row,val];
    dlmwrite(strcat('data/u',num2str(i)),data_dump,'precision',16);
    
    %   Write out vi = v in a sparse format
    [n,m] = size(v);
    v(abs(v) < 10^-10) = 0;
    v(n,m) = 10^-15;
    [col,row,val] = find(v);
    data_dump = [col,row,val];
    dlmwrite(strcat('data/v',num2str(i)),data_dump,'precision',16);
    
    %   Print ranks and stuff    
    %   rank(full(GY - DAE.Gy))
    %   rank(full(GX - DAE.Gx)) 
    %   rank(diff)
end

%%   Delete workspace files

for i = 1:numcontigs
    name = strcat(basefilename,'_',num2str(i),'.m');
    delete(name);
    name = strcat('contig_',basefilename,'_',num2str(i),'.m');
    delete(name);
    
end

basefilename = 'bus14';
filename = strcat(basefilename,'_');
numcontigs = 20;



for i = 1:numcontigs
     temp = strcat(filename,num2str(i));
     fid = fopen(strcat(temp,'.m'),'w');
     fid = fclose(fid);
     copyfile(strcat(basefilename,'.m'),strcat(temp,'.m'));
     cssm(strcat(temp,'.m'),i+1);
end

initpsat
runpsat(basefilename,'data');
runpsat('pf');
Aoriginal = DAE.Fx - DAE.Fy*(DAE.Gy\DAE.Gx) - 1e-6*speye(DAE.n);
GY = DAE.Gy; GX = DAE.Gx;
Aoriginal = full(Aoriginal);

filename = strcat('contig_',filename);
for i = 1:numcontigs
    temp = strcat(filename,num2str(i));
    runpsat(temp,'data');
    runpsat('pf');
    As = DAE.Fx - DAE.Fy*(DAE.Gy\DAE.Gx) - 1e-6*speye(DAE.n);
    rank(full(GY - DAE.Gy))
    rank(full(GX - DAE.Gx))
    diff = Aoriginal - full(As);
    rank(diff)
    [u,d,v] = svd(diff);
    dlmwrite(strcat('u',num2str(i)), u*d);
    dlmwrite(strcat('v',num2str(i)), v);
end


for i = 1:numcontigs;
    
end
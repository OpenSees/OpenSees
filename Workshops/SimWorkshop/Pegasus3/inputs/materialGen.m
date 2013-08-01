# open the input file
fid=fopen("materialProperties.in","r");

# read steel properties
numOut=str2num(fgetl(fid))
numEle=str2num(fgetl(fid))
ksi=str2num(fgetl(fid))
EMean=str2num(fgetl(fid))
EStdDev=str2num(fgetl(fid))
Etype=fgetl(fid)
fyMean=str2num(fgetl(fid))
fyStdDev=str2num(fgetl(fid))
fyType=fgetl(fid)
bMean=str2num(fgetl(fid))
bStdDev=str2num(fgetl(fid))
bType=fgetl(fid)


# create array of random numbers with mean and std dev
if (strcmp(Etype,'n'))
  a=normrnd(fyMean, fyStdDev, numOut*numEle,1);
 else
   m=fyMean;
   v= fyStdDev*fyStdDev;
   mu=log(m^2/sqrt(v+m^2));
   sigma=sqrt(log((v/m^2) + 1));
   a=lognrnd(mu, sigma, numOut*numEle,1);
end

if (strcmp(fyType,'n'))
  b=normrnd(bMean, bStdDev, numOut*numEle,1);
 else
   m=bMean;
   v= bStdDev*bStdDev;
   mu=log(m^2/sqrt(v+m^2));
   sigma=sqrt(log((v/m^2) + 1));
   b=lognrnd(mu, sigma, numOut*numEle,1);
end

if (strcmp(fyType,'n'))
  c=normrnd(EMean, EStdDev, numOut*numEle,1);
 else
   m=EMean;
   v= EStdDev*EStdDev;
   mu=log(m^2/sqrt(v+m^2));
   sigma=sqrt(log((v/m^2) + 1));
   c=lognrnd(mu, sigma, numOut*numEle,1);
end

# create output files
count = 1
for i=1:numOut
	outFileName=sprintf('matProperties.tcl.%d',i);
        outFile=fopen(outFileName,"w");
        if i == 1
	  for j=1:numEle
		  fprintf(outFile,'uniaxialMaterial Steel01 %d %f %f %f\n',j, fyMean, EMean, bMean);
          end
	else
	  for j=1:numEle
		  fprintf(outFile,'uniaxialMaterial Steel01 %d %f %f %f\n',j,a(count),c(count),b(count));
             count=count+1;
          end
        end
        fclose(outFile);
end

# close the file
fclose(fid);

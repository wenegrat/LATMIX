% Reads in statistics outputted by diablo
clear turbo;
clear n;

turbo=dlmread('../../latmix_front_tracer1/tke_les.txt');

% Set the domain size
%NY=201;
%NYM=NY-1;

% Enter the viscosity
%NU=1/60000;

% Set the starting time in code units for start of averaging
%tstart=0;

% Determine the number of records in the file based on its length
nk=floor(length(turbo(:,1))/(NY))

row=1;
for k=1:nk
  tii(k)=turbo(row,2);
  dt(k)=turbo(row,3);
  ubulk(k)=turbo(row,4);
  row=row+1;
  for j=2:NY
    gy_tke(j)=turbo(row,2);
    tke_8(j,k)=turbo(row,3);
    tke_9(j,k)=turbo(row,4);
    
    row=row+1;
  end
end


% Get the time index based on start time
kstart=0;
for k=1:nk
  if (tii(k) <= tstart)
     kstart=k;
  end
end
if (kstart == 0)
  kstart=1;
end
'Start of time average: ',tii(kstart)

for j=1:NY
  tke_8_mean(j)=mean(tke_8(j,kstart:nk));
  tke_9_mean(j)=mean(tke_9(j,kstart:nk));
end
  
   
%for j=2:NY
%  gy(j)=(gyf(j)+gyf(j-1))/2;    
%end
%gy(1)=(gyf(1)+0)/2;
%gy(NY)=(1+gyf(NYM))/2;
%for j=1:NYM
%  dyf(j)=(gy(j+1)-gy(j));
%end



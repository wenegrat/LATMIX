% Reads in statistics outputted by diablo
% Run after readmean.m

clear turbo;
clear n;

turbo=dlmread('../front/run1/tke.txt');

n=1;
row=1;
for k=1:nk
  tii(k)=turbo(row,2);
  dt(k)=turbo(row,3);
  ubulk(k)=turbo(row,4);
  row=row+1;
  for j=1:NY
    gy_tke(j)=turbo(row,2);
    tke_1(j,k)=turbo(row,3);
    tke_2(j,k)=turbo(row,4);
    tke_3(j,k)=turbo(row,5);
    tke_4(j,k)=turbo(row,6);
    for n=1:N_TH
      tke_5(j,k,n)=turbo(row,6+n);
    end
    if (N_TH==0) 
      n=0;
    end
    tke_6_1(j,k)=turbo(row,6+n+1);
    tke_6_1_1(j,k)=turbo(row,6+n+2);
    tke_6_1_2(j,k)=turbo(row,6+n+3);
    tke_6_1_3(j,k)=turbo(row,6+n+4);
    tke_6_2(j,k)=turbo(row,6+n+5);
    tke_6_2_1(j,k)=turbo(row,6+n+6);
    tke_6_2_2(j,k)=turbo(row,6+n+7);
    tke_6_2_3(j,k)=turbo(row,6+n+8);
    tke_7(j,k)=turbo(row,6+n+9);
    tke_3_1(j,k)=turbo(row,6+n+10);
    tke_3_2(j,k)=turbo(row,6+n+11);
    tke_3_3(j,k)=turbo(row,6+n+12);
    tke_7_1(j,k)=turbo(row,6+n+13);
    tke_7_2(j,k)=turbo(row,6+n+14);
    tke_7_3(j,k)=turbo(row,6+n+15);
    tke_8_1(j,k)=turbo(row,6+n+16);
    tke_8_2(j,k)=turbo(row,6+n+17);
    tke_8_3(j,k)=turbo(row,6+n+18);
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

eta(1:NY,1:nk)=0;
for j=1:NY
for k=1:nk
  eta(j,k)=abs(tke_7(j,k))^(-1/4)*(sqrt(NU));
end
end

%for n=1:N_TH
%for j=1:NY
%for k=1:nk
%  if (brunt(j,k,n)~=0)
%  loz(j,k)=sqrt(-(tke_7(j,k)+nu_t_sgs(j,k)*tke_7(j,k)/NU)/brunt(j,k,n)^3);
%  end
%end
%end
%end

for j=1:NY
  tke_1_mean(j)=mean(tke_1(j,kstart:kend));
  tke_2_mean(j)=mean(tke_2(j,kstart:kend));
  tke_3_mean(j)=mean(tke_3(j,kstart:kend));
  tke_3_1_mean(j)=mean(tke_3_1(j,kstart:kend));
  tke_3_2_mean(j)=mean(tke_3_2(j,kstart:kend));
  tke_3_3_mean(j)=mean(tke_3_3(j,kstart:kend));
  tke_4_mean(j)=mean(tke_4(j,kstart:kend));
  for n=1:N_TH
    tke_5_mean(j,n)=mean(tke_5(j,kstart:kend,n));
  end
  if (N_TH==0)
    tke_5_mean(j,1)=0;
  end
  tke_6_1_mean(j)=mean(tke_6_1(j,kstart:kend));
  tke_6_1_1_mean(j)=mean(tke_6_1_1(j,kstart:kend));
  tke_6_1_2_mean(j)=mean(tke_6_1_2(j,kstart:kend));
  tke_6_1_3_mean(j)=mean(tke_6_1_3(j,kstart:kend));
  tke_6_2_mean(j)=mean(tke_6_2(j,kstart:kend));
  tke_6_2_1_mean(j)=mean(tke_6_2_1(j,kstart:kend));
  tke_6_2_2_mean(j)=mean(tke_6_2_2(j,kstart:kend));
  tke_6_2_3_mean(j)=mean(tke_6_2_3(j,kstart:kend));
  tke_7_mean(j)=mean(tke_7(j,kstart:kend));
  tke_7_1_mean(j)=mean(tke_7_1(j,kstart:kend));
  tke_7_2_mean(j)=mean(tke_7_2(j,kstart:kend));
  tke_7_3_mean(j)=mean(tke_7_3(j,kstart:kend));
  tke_8_1_mean(j)=mean(tke_8_1(j,kstart:kend));
  tke_8_2_mean(j)=mean(tke_8_2(j,kstart:kend));
  tke_8_3_mean(j)=mean(tke_8_3(j,kstart:kend));
  eta_mean(j)=mean(eta(j,kstart:kend));
end
 
   
%for j=2:NY
%  gy(j)=(gyf(j)+gyf(j-1))/2;    
%end
%gy(1)=(gyf(1)+0)/2;
%gy(NY)=(1+gyf(NYM))/2;
%for j=1:NYM
%  dyf(j)=(gy(j+1)-gy(j));
%end

% Now, compute the potential energy in the system

pe(1:NY,1:nk)=0;
for n=1:N_TH
for k=1:nk
  pe(j,k)=0;
  pe_total(k)=0;
  for j=2:NY-2
    pe(j,k)=pe(j,k)-RI(1)*thme(j,k,n)*gyf(j)*(gy(j+1)-gy(j));
    pe_total(k)=pe_total(k)-RI(1)*thme(j,k,n)*gyf(j)*(gy(j+1)-gy(j));
  end
end
for k=2:nk-2
  dpedt(k)=(pe_total(k)-pe_total(k-1))/(tii(k)-tii(k-1));
end

for k=1:nk
  pe_drift(k)=pe_total(1)-0.5*(gy(NY-1)-gy(2))^2*RI(n)*(tii(k)-tii(1))/(5/NU);
end
end


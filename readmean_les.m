% Reads in statistics outputted by diablo from the LES subroutine
% Run after readmean.m

turbo=dlmread('../front/run1/mean_les.txt');

row=1;
for k=1:nk
  tii(k)=turbo(row,2);
  dt(k)=turbo(row,3);
  row=row+1;
  for j=1:NY
    gyf(j)=turbo(row,2);
    tau11(j,k)=turbo(row,3); 
    tau22(j,k)=turbo(row,4); 
    tau33(j,k)=turbo(row,5); 
    tau12(j,k)=turbo(row,6); 
    tau13(j,k)=turbo(row,7); 
    tau23(j,k)=turbo(row,8); 
    nu_t_sgs(j,k)=turbo(row,9);
    nu_u1(j,k)=turbo(row,10);
    nu_u3(j,k)=turbo(row,11);
    cdyn(j,k)=turbo(row,12);
    cdyn_h(j,k)=turbo(row,13);
    cdyn_v(j,k)=turbo(row,14);
    row=row+1;
  end
end

% Add the eddy viscosity part to tau_ij
% This part was treated implicitly and was not saved in tau_ij
for k=1:nk
  for j=1:NY
    tau12(j,k)=tau12(j,k)-nu_t_sgs(j,k)*dudy(j,k);
    tau23(j,k)=tau23(j,k)-nu_t_sgs(j,k)*dwdy(j,k);
  end
end

for k=1:nk
  for j=1:NY
    urms_sgs(j,k)=sqrt(abs(tau11(j,k)));
    vrms_sgs(j,k)=sqrt(abs(tau22(j,k)));
    wrms_sgs(j,k)=sqrt(abs(tau33(j,k)));
    uv_sgs(j,k)=tau12(j,k);
    uw_sgs(j,k)=tau13(j,k);
    wv_sgs(j,k)=tau23(j,k);
%    if (dudy(j,k)~=0) 
%      nu_t_sgs(j,k)=-uv_sgs(j,k)/dudy(j,k);  %readmean.m must be done first
%    else
%      nu_t_sgs(j,k)=0;
%    end
    tke_sgs(j,k)=0.5*(abs(tau11(j,k))+abs(tau22(j,k))+abs(tau33(j,k)));
  end
end

for k=1:nk
  for j=2:NY-1
    for n=1:N_TH
        if ((RI(1)*(thme(j,k,n)-thme(j-1,k,n))/(gyf(j)-gyf(j-1))*eta_rms(j,k))~=0)
        frv_sgs(j,k,n)=vrms_sgs(j,k)/(RI(1)*(thme(j,k,n)-thme(j-1,k,n))/(gyf(j)-gyf(j-1))*eta_rms(j,k));
        else
        frv_sgs(j,k,n)=0;
        end
    end
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
  urms_sgs_mean(j)=mean(urms_sgs(j,kstart:kend));
  vrms_sgs_mean(j)=mean(vrms_sgs(j,kstart:kend));
  wrms_sgs_mean(j)=mean(wrms_sgs(j,kstart:kend));
  uv_sgs_mean(j)=mean(uv_sgs(j,kstart:kend));
  uw_sgs_mean(j)=mean(uw_sgs(j,kstart:kend));
  wv_sgs_mean(j)=mean(wv_sgs(j,kstart:kend));
  tke_sgs_mean(j)=mean(tke_sgs(j,kstart:kend));
%  if (dudy_mean(j)~=0)
%    nu_t_sgs_mean(j)=sqrt(uv_sgs_mean(j)^2+wv_sgs_mean(j)^2)/sqrt(dudy_mean(j)^2+dwdy_mean(j)^2);
%    nu_t_total_mean(j)=sqrt((uv_mean(j)+uv_sgs_mean(j))^2+(wv_mean(j)+wv_sgs_mean(j))^2)/sqrt(dudy_mean(j)^2+dwdy_mean(j)^2);
%  end
  nu_t_sgs_mean(j)=mean(nu_t_sgs(j,kstart:kend));
  cdyn_mean(j)=mean(cdyn(j,kstart:kend));
end

for n=1:N_TH
for j=2:NY-1
  frv_sgs_mean(j,n)=mean(frv_sgs(j,kstart:kend,n));
end
end

beta=3;
LX=0.1454;
NX=128;
LZ=0.1454;
NZ=128;
for j=2:NY-1
  delta_yf(j)=sqrt((beta*LX/NX)^2+(2*dyf(j))^2+(beta*LZ/NZ)^2);
  cdyn_only(j)=sqrt(cdyn_mean(j)/(delta_yf(j)^2));
end


%% Reads in statistics outputted by diablo 
% turbo=dlmread('../front/run2/mean.txt');
% turbo_th=dlmread('../front/run2/mean_th.txt');
disp(path)
turbo=dlmread([path, 'mean.txt']);
turbo_th=dlmread([path, 'mean_th.txt']);
% Set the domain size
NY=100;
LY=1.0;
NYM=NY-1;

% Enter the viscosity
NU=1e-6;
% Prandtl number
PR=1;
% For latmix
f=9.3e-5;
zeta=0.8*f;
f_eff=sqrt(f*(f+zeta));

% Enter the number of scalars
N_TH=3;

% Enter the richardson number for each scalar
RI(1)=1.0;
RI(2)=0.0;
RI(3)=0.0;

% Set the background lateral density gradient
% For latmix
dbdx=-5e-7;
dbdz=0;

% Thermal wind balance
dwgdy=RI(1)*dbdx/f;
dugdy=-RI(1)*dbdz/f;

% Determine the number of records in the file based on its length
nk=ceil(length(turbo(:,1))/(NY+3))

row=1;
for k=1:nk
  tii(k)=turbo(row,2);
  dt(k)=turbo(row,3);
  row=row+1;
  ubulk(k)=turbo(row,1);
  row=row+1;
  for j=1:NY
    gyf(j)=turbo(row,2);
    ume(j,k)=turbo(row,3);
    vme(j,k)=turbo(row,4);
    wme(j,k)=turbo(row,5);
    urms(j,k)=turbo(row,6);
    vrms(j,k)=turbo(row,7);
    wrms(j,k)=turbo(row,8);
    uv(j,k)=turbo(row,9);
    duv(j,k)=turbo(row,10);
    uw(j,k)=turbo(row,11);
    wv(j,k)=turbo(row,12);
    dwv(j,k)=turbo(row,13);
    dudy(j,k)=turbo(row,14);
    dwdy(j,k)=turbo(row,15);
    cp(j,k)=turbo(row,16);
    shear(j,k)=turbo(row,17);
    omega_x(j,k)=turbo(row,18);
    omega_y(j,k)=turbo(row,19);
    omega_z(j,k)=turbo(row,20);
    prms(j,k)=turbo(row,21);  
    pv(j,k)=turbo(row,22); 
    row=row+1;
  end
end

for j=1:NY
  for k=1:nk
    vel(j,k)=sqrt(ume(j,k)^2+(wme(j,k))^2);
  end
end

% Add background geostrophic shear to shear:
for j=1:NY
	for k=1:nk
		shear(j,k)=sqrt((dudy(j,k)+dugdy)^2+(dwdy(j,k)+dwgdy)^2);
	end
end


% Now read in the scalar statistics
if (N_TH>0)
% Determine the number of records in the file based on its length
%  nk_th=ceil(length(turbo_th(:,1))/(N_TH*(NY+3)))
  nk_th=nk;
else
  nk_th=0;
end


row=1;
for k=1:nk_th
  tii(k)=turbo_th(row,2);
  dt(k)=turbo_th(row,3);
  row=row+1;
  ubulk(k)=turbo_th(row,1);
  row=row+1;
  for n=1:N_TH
  for j=1:NY
    thme(j,k,n)=turbo_th(row,3);
    dthdy(j,k,n)=turbo_th(row,4); % Add one here if a background was subtracted
    thrms(j,k,n)=turbo_th(row,5);
    thv(j,k,n)=turbo_th(row,6); 
    dthv(j,k,n)=turbo_th(row,7);
% For new versions only
    if (RI(n) ~= 0) 
      pe_diss(j,k,n)=turbo_th(row,8)*(NU/PR)/RI(n);
    else
      pe_diss(j,k,n)=0;
    end
    qme(j,k,n)=turbo_th(row,9);
    qv(j,k,n)=turbo_th(row,10); 
    thu(j,k,n)=turbo_th(row,11);
    thw(j,k,n)=turbo_th(row,12);
    qme_vert(j,k,n)=turbo_th(row,13);
    qme_bc(j,k,n)=turbo_th(row,14);
    thth(j,k,n)=turbo_th(row,15);
    row=row+1;
  end
  end
end

% Compute secondary quantities
for k=1:nk
  for j=1:NY
    tke(j,k)=0.5*(urms(j,k)^2.+vrms(j,k)^2.+wrms(j,k)^2.);
    if (dudy(j,k)~=0)
      nu_t(j,k)=-(uv(j,k)+wv(j,k))/sqrt(dudy(j,k)^2+(dwdy(j,k)+dwgdy)^2);
    else
      nu_t(j,k)=0;
    end
    if ((urms(j,k)~=0)&&(vrms(j,k)~=0))
      uv_corr(j,k)=uv(j,k)/(urms(j,k)*vrms(j,k));
    else
      uv_corr(j,k)=0;
    end
    if ((wrms(j,k)~=0)&&(vrms(j,k)~=0))
      wv_corr(j,k)=wv(j,k)/(wrms(j,k)*vrms(j,k));
    else
      wv_corr(j,k)=0;
    end
    for n=1:N_TH
      if ((thrms(j,k,n)~=0)&&(vrms(j,k)~=0))
        thv_corr(j,k,n)=thv(j,k,n)/(thrms(j,k,n)*vrms(j,k));
      else
        thv_corr(j,k,n)=0;
      end
      if (dthdy(j,k,n)~=0) 
        eta_rms(j,k,n)=thrms(j,k)/dthdy(j,k,n); 
      else
        eta_rms(j,k,n)=0;
      end
    end

% Calculate the vertical taylor scale
    if (shear(j,k)~=0) 
      taylor(j,k)=sqrt((ume(j,k)^2.+wme(j,k)^2.+urms(j,k)^2.+wrms(j,k)^2.)/shear(j,k));
    else
      taylor(j,k)=0;
    end

    if (N_TH > 0)
      for n=1:N_TH
        brunt(j,k,n)=sqrt(RI(n)*(dthdy(j,k,n))); 
        if (shear(j,k)~=0) 
          grarich(j,k,n)=brunt(j,k,n)^2./shear(j,k)^2; 
        else
          grarich(j,k,n)=0;
        end
        if (dthdy(j,k,n)~=0)
          kappa_t(j,k,n)=-thv(j,k,n)/dthdy(j,k,n);
        else
          kappa_t(j,k,n)=0;
        end 

      end

    end
  end
end

for k=1:nk
  for j=2:NY-1
    for n=1:N_TH
        if ((RI(n)*(thme(j,k,n)-thme(j-1,k,n))/(gyf(j)-gyf(j-1))*eta_rms(j,k,n))~=0)
        frv(j,k,n)=vrms(j,k)/(RI(n)*(thme(j,k,n)-thme(j-1,k,n))/(gyf(j)-gyf(j-1))*eta_rms(j,k,n));
        else
        frv(j,k,n)=0;
        end
    end
  end
end

for j=1:NY
for k=1:nk
ri_bulk(j,k)=brunt(j,k)^2*f^2/((RI(1)*dbdz)^2+(RI(1)*dbdx)^2);
end
end

tstart=0;
tend=(67-64.5)*3600*24;

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

kend=0;
for k=1:nk
  if (tii(k) <= tend)
     kend=k;
  end
end
if (kend == 0)
  kend=1;
end


'Start of time average: ',tii(kstart)
'End of time average: ',tii(kend)


for j=1:NY
  ume_mean(j)=mean(ume(j,kstart:kend));
  vme_mean(j)=mean(vme(j,kstart:kend));
  wme_mean(j)=mean(wme(j,kstart:kend));
  urms_mean(j)=mean(urms(j,kstart:kend));
  vrms_mean(j)=mean(vrms(j,kstart:kend));
  wrms_mean(j)=mean(wrms(j,kstart:kend));
  dudy_mean(j)=mean(dudy(j,kstart:kend));
  dwdy_mean(j)=mean(dwdy(j,kstart:kend));
  tke_mean(j)=mean(tke(j,kstart:kend));
  uv_mean(j)=mean(uv(j,kstart:kend));
  wv_mean(j)=mean(wv(j,kstart:kend));
  qme_mean(j)=mean(qme(j,kstart:kend));
  uv_corr_mean(j)=mean(uv_corr(j,kstart:kend));
  wv_corr_mean(j)=mean(wv_corr(j,kstart:kend));
  vel_mean(j)=sqrt(ume_mean(j)^2+wme_mean(j)^2);
  rib_mean(j)=mean(ri_bulk(j,kstart:kend));
  prms_mean(j)=mean(prms(j,kstart:nk));
  pv_mean(j)=mean(pv(j,kstart:nk));

  if (dudy_mean(j)~=0) 
    nu_t_mean(j)=sqrt(uv_mean(j)^2+wv_mean(j)^2)/sqrt(dudy_mean(j)^2+(dwdy_mean(j)+dwgdy)^2);
  else
    nu_t_mean(j)=0;
  end
  for n=1:N_TH
    thv_mean(j,n)=mean(thv(j,kstart:kend,n));
    thu_mean(j,n)=mean(thu(j,kstart:kend,n));
    thw_mean(j,n)=mean(thw(j,kstart:kend,n));
    thv_corr_mean(j,n)=mean(thv_corr(j,kstart:kend,n));
    dthdy_mean(j,n)=mean(dthdy(j,kstart:kend,n));
    thrms_mean(j,n)=mean(thrms(j,kstart:kend,n));
    thbar(j,n)=mean(thme(j,kstart:kend,n));
    pe_diss_mean(j,n)=mean(pe_diss(j,kstart:kend,n));
    if (dthdy_mean(j,n)~=0) 
      kappa_t_mean(j,n)=-thv_mean(j,n)/dthdy_mean(j,n);
    else
      kappa_t_mean(j,n)=0;
    end 
    eta_rms_mean(j,n)=mean(eta_rms(j,kstart:kend,n));
  end
  shear_mean(j)=mean(shear(j,kstart:kend));
  omega_x_mean(j)=mean(omega_x(j,kstart:kend));
  omega_y_mean(j)=mean(omega_y(j,kstart:kend));
  omega_z_mean(j)=mean(omega_z(j,kstart:kend));
end
  


for j=2:NY-1
  for n=1:N_TH
    frv_mean(j,n)=mean(frv(j,kstart:kend,n));
  end
end
 
for j=2:NY
  gy(j)=(gyf(j)+gyf(j-1))/2;    
end
for j=2:NY-1
  dyf(j)=(gy(j+1)-gy(j));
end
dyf(NY)=dyf(NY-1);
for j=2:NY
  dy(j)=gyf(j)-gyf(j-1);
end

for j=2:NY-1
  ry(j)=(gyf(j+1)-gyf(j))/(gyf(j)-gyf(j-1));
end

for n=1:N_TH
  for j=1:NY-1
    for k=1:nk
      dthv2(j,k,n)=(thv(j+1,k,n)-thv(j,k,n))/dyf(j);
    end
  end
end
for j=1:NY-1
  for k=1:nk
    dwv(j,k)=(wv(j+1,k)-wv(j,k))/dyf(j);
    duv(j,k)=(uv(j+1,k)-uv(j,k))/dyf(j);
  end
end

for n=1:N_TH
  for j=2:NY-1
    for k=1:nk
      d2thv(j,k,n)=(dthv(j,k,n)-dthv(j-1,k,n))/(gyf(j)-gyf(j-1));
    end
  end
end
for n=1:N_TH
  for j=2:NY-1
    for k=1:nk
      thv_ratio2(j,k,n)=abs(d2thv(j,k,n))/(abs(dbdx*(ume(j,k)-ume(j-1,k))/(gyf(j)-gyf(j-1)))+abs(d2thv(j,k,n)));
    end
  end
end

for k=1:nk
surface_th(k)=thme(NY-1,k)+(thme(NY,k)-thme(NY-1,k))*(gy(NY)-gyf(NY-1))/(gyf(NY)-gyf(NY-1));
end


k1=1;
k2=nk;
for j=1:NY
  thv_int(j)=0;
  dudy_int(j)=0;
  for k=k1:k2-1
    thv_int(j)=thv_int(j)+0.5*(thv(j,k)+thv(j,k+1))*(tii(k+1)-tii(k));
    dudy_int(j)=dudy_int(j)+0.5*(dudy(j,k)+dudy(j,k+1))*(tii(k+1)-tii(k));
  end
end

k1=1;
k2=nk;
for j=1:NY
  wme_int(j)=0;
  uv_int(j)=0;
  for k=k1:k2-1
    wme_int(j)=wme_int(j)+0.5*(wme(j,k)+wme(j,k+1))*(tii(k+1)-tii(k));
    uv_int(j)=uv_int(j)+0.5*(uv(j,k)+uv(j,k+1))*(tii(k+1)-tii(k));
  end
end

for j=2:NY-1
  duv_int(j)=(uv_int(j+1)-uv_int(j))/dyf(j);
end

% dthv_int defined at GYF points
for j=1:NY-1
  dthv_int(j)=(thv_int(j+1)-thv_int(j))/dyf(j);
end

%d2thv_int defined at GY points
for j=2:NY-1
  d2thv_int(j)=(dthv_int(j)-dthv_int(j-1))/dy(j);
end

for j=1:NY
for k=1:nk
ellison(j,k,:)=thrms(j,k,:)./dthdy(j,k,:);
end
end

for j=1:NY
elli_mean(j,:)=thrms_mean(j,:)./dthdy_mean(j,:);
end

for k=1:nk-1
  dthmedt(:,k,:)=(thme(:,k+1,:)-thme(:,k,:))/(tii(k+1)-tii(k));
end

for j=1:NY
for n=1:N_TH
thme_mean(j,n)=mean(thme(j,kstart:kend,n));
end
end

GSP_int(1)=0;
ASP_int(1)=0;
for k=2:nk
GSP_int(k)=trapz(tii(1:k),mean(-wv(:,1:k)*dwgdy,1));
ASP_int(k)=trapz(tii(1:k),mean(-uv(:,1:k).*dudy(:,1:k)-wv(:,1:k).*dwdy(:,1:k),1));
end

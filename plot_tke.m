

for j=2:NY-1
  tke_1_sum(j)=0.5*(tke_1_mean(j)+tke_1_mean(j+1));
end
plot(-1*tke_1_sum(2:NY-1),gy_tke(2:NY-1),'k-');
hold on
%plot(tke_2_mean(2:NY),gy_tke(2:NY),'k:');
for j=2:NY-1
% Here, the SGS geostrophic production is subtracted since it wasn't accounted for in tke_8
%  tke_3_sum(j)=tke_3_1_mean(j)+0.5*(tke_3_2_mean(j)+tke_3_2_mean(j+1))+tke_3_3_mean(j)-wv_mean(j)*dvgdz-nu_t_sgs_mean(j)*dwdy_mean(j)*dvgdz;
  tke_3_sum(j)=tke_3_1_mean(j)+0.5*(tke_3_2_mean(j)+tke_3_2_mean(j+1))+tke_3_3_mean(j)-wv_mean(j)*dvgdz;
end
plot(tke_3_sum(2:NY-1),gy_tke(2:NY-1),'b.-');
%plot(nu_t_sgs_mean(2:NY-1).*dwdy_mean(2:NY-1)*dvgdz,gy_tke(2:NY-1),'c*-');
%plot(nu_t_sgs_mean(2:NY-1).*(dudy_mean(2:NY-1).^2+dwdy_mean(2:NY-1).^2),gy_tke(2:NY-1),'b--');
plot(-(wv_mean(2:NY-1))*dvgdz,gy_tke(2:NY-1),'c.-');
%plot(-wv_mean(2:NY-1)*dvgdz+nu_t_sgs_mean(2:NY-1).*dwdy_mean(2:NY-1)*dvgdz+tke_5_mean(2:NY-1)',gy_tke(2:NY-1),'g-');
for j=2:NY
%  tke_5_sum(j)=tke_5_mean(j)-1*RI(1)*nu_t_sgs_mean(j)*dthdy_mean(j);
  tke_5_sum(j)=tke_5_mean(j);
end
plot(tke_5_sum(2:NY)',gy_tke(2:NY),'rs-');
%plot(-1*RI(1)*nu_t_sgs_mean(2:NY).*dthdy_mean(2:NY)',gy_tke(2:NY),'r--');
for j=2:NY-1
  tke_6_1_sum(j)=tke_6_1_1_mean(j)+0.5*(tke_6_1_2_mean(j)+tke_6_1_2_mean(j+1))+tke_6_1_3_mean(j);
end
%plot(tke_6_1_sum(2:NY-1),gyf(2:NY-1),'g^-');
for j=2:NY-1
  tke_6_2_sum(j)=tke_6_2_1_mean(j)+0.5*(tke_6_2_2_mean(j)+tke_6_2_2_mean(j+1))+tke_6_2_3(j,5);
end
%plot(tke_6_1_mean(2:NY-1),gyf(2:NY-1),'gv-');
%plot(tke_6_2_mean(2:NY-1),gyf(2:NY-1),'g^-');
for j=2:NY
dtke(j,:)=(tke(j,:)-tke(j-1,:))/(gyf(j)-gyf(j-1));
end
for j=2:NY-1
d2tke(j,:)=(dtke(j+1,:)-dtke(j,:))/(gy(j+1)-gy(j));
end
for j=2:NY-1
d2tke_mean(j)=mean(d2tke(j,kstart:kend));
end
for j=2:NY-1
  tke_6_3_mean(j)=d2tke_mean(j)*0.5*(nu_t_sgs_mean(j)+nu_t_sgs_mean(j+1));
end
for j=2:NY-1
  tke_7_sum(j)=tke_7_1_mean(j)+0.5*(tke_7_2_mean(j)+tke_7_2_mean(j+1))+tke_7_3_mean(j);
end
%plot(tke_7_sum(2:NY-1),gyf(2:NY-1),'r*-');
for j=2:NY-1
  tke_8_sum(j)=tke_8_1_mean(j)+0.5*(tke_8_2_mean(j)+tke_8_2_mean(j+1))+tke_8_3_mean(j);
end

for j=2:NY-2
  tke_sum(j)=-1*tke_1_sum(j)+0.5*(tke_2_mean(j)+tke_2_mean(j+1))+0.5*(tke_3_sum(j)+tke_3_sum(j+1))+0.5*(tke_5_sum(j)+tke_5_sum(j+1))+tke_6_1_mean(j)+tke_6_2_mean(j)+tke_6_3_mean(j)+tke_7_sum(j)+tke_8_sum(j);
end
j=NY-1
tke_sum(j)=-1*tke_1_sum(j)+0.5*(tke_2_mean(j)+tke_2_mean(j+1))+0.5*(tke_3_sum(j)+0)+0.5*(tke_5_sum(j)+tke_5_sum(j+1))+tke_6_1_sum(j)+tke_6_2_sum(j)+tke_7_sum(j)+tke_8_sum(j);

plot(tke_8_sum(2:NY-1),gyf(2:NY-1),'m*-');
plot(-tke_sum(2:NY-1),gyf(2:NY-1),'k*-');
plot(tke_6_1_mean(2:NY-1)+tke_6_2_mean(2:NY-1)+tke_6_3_mean(2:NY-1),gyf(2:NY-1),'gs-');

nk=length(tke_8_1(1,:));
for j=2:NY-1
  for k=1:nk
%    epsilon_si(j,k)=tke_8_1(j,k)+0.5*(tke_8_2(j,k)+tke_8_2(j+1,k))+tke_8_3(j,k)+tke_3_1(j,k)+0.5*(tke_3_2(j,k)+tke_3_2(j+1,k))+tke_3_3(j,k)+tke_6_1_1(j,k)+0.5*(tke_6_1_2(j,k)+tke_6_1_2(j+1,k))+tke_6_1_3(j,k)+tke_6_2_1(j,k)+0.5*(tke_6_2_2(j,k)+tke_6_2_2(j+1,k))+tke_6_2_3(j,k);
    epsilon_si(j,k)=tke_8_1(j,k)+0.5*(tke_8_2(j,k)+tke_8_2(j+1,k))+tke_8_3(j,k)+tke_3_1(j,k)+0.5*(tke_3_2(j,k)+tke_3_2(j+1,k))+tke_3_3(j,k);
    epsilon(j,k)=tke_8_1(j,k)+0.5*(tke_8_2(j,k)+tke_8_2(j+1,k))+tke_8_3(j,k) ...
        +tke_7_1(j,k)+0.5*(tke_7_2(j,k)+tke_7_2(j+1,k))+tke_7_3(j,k);
    P_AG(j,k)=tke_3_1(j,k)+0.5*(tke_3_2(j,k)+tke_3_2(j+1,k))+tke_3_3(j,k);
  end
end

disp('paused');
pause

for j =2:NY-1
  epsilon_si_mean(j)=mean(epsilon_si(j,kstart:kend));
  epsilon_si_std(j)=std(epsilon_si(j,kstart:kend));
  epsilon_mean(j)=mean(epsilon(j,kstart:kend));
  epsilon_std(j)=std(epsilon(j,kstart:kend));
end
for j=2:NY-1
  tke_3_std(j)=std(tke_3_1(j,kstart:kend)+0.5*(tke_3_2(j,kstart:kend)+tke_3_2(j+1,kstart:kend))+tke_3_3(j,kstart:kend));
  tke_5_std(j)=std(tke_5(j,kstart:kend));
  tke_1_std(j)=std(tke_1(j,kstart:kend));
  wv_std(j)=std(wv(j,kstart:kend));
end
figure
plot(tke_5_sum(2:NY-1),gy_tke(2:NY-1),'r-');
hold on
plot(tke_1_sum(2:NY-1),gy_tke(2:NY-1),'k-');
plot(tke_3_sum(2:NY-1),gy_tke(2:NY-1),'c-');
%plot(epsilon_si_mean(2:NY-1),gyf(2:NY-1),'g-');
plot(epsilon_mean(2:NY-1),gyf(2:NY-1),'m-');
plot(-wv_mean(2:NY-1)*dvgdz,gy_tke(2:NY-1),'b-');
plot(tke_5_sum(2:NY-1)-tke_5_std(2:NY-1),gy_tke(2:NY-1),'r:');
plot(tke_5_sum(2:NY-1)+tke_5_std(2:NY-1),gy_tke(2:NY-1),'r:');
plot(tke_1_sum(2:NY-1)-tke_1_std(2:NY-1),gy_tke(2:NY-1),'k:');
plot(tke_1_sum(2:NY-1)+tke_1_std(2:NY-1),gy_tke(2:NY-1),'k:');
plot(tke_3_sum(2:NY-1)-tke_3_std(2:NY-1),gy_tke(2:NY-1),'c:');
plot(tke_3_sum(2:NY-1)+tke_3_std(2:NY-1),gy_tke(2:NY-1),'c:');

%plot(epsilon_si_mean(2:NY-1)+epsilon_si_std(2:NY-1),gyf(2:NY-1),'g:');
%plot(epsilon_si_mean(2:NY-1)-epsilon_si_std(2:NY-1),gyf(2:NY-1),'g:');
plot(epsilon_mean(2:NY-1)+epsilon_std(2:NY-1),gyf(2:NY-1),'m:');
plot(epsilon_mean(2:NY-1)-epsilon_std(2:NY-1),gyf(2:NY-1),'m:');
plot(-wv_mean(2:NY-1)*dvgdz-wv_std(2:NY-1)*dvgdz,gy_tke(2:NY-1),'b:');
plot(-wv_mean(2:NY-1)*dvgdz+wv_std(2:NY-1)*dvgdz,gy_tke(2:NY-1),'b:');


ume_int(1:NY,1:nk)=0;
thme_int(1:NY,1:nk)=0;
wme_int(1:NY,1:nk)=0;

for k=1:nk
  for j=1:NY-2
    ume_int(j,k)=trapz(gyf(j:NY-1),ume(j:NY-1,k));
    wme_int(j,k)=trapz(gyf(j:NY-1),wme(j:NY-1,k));
    thme_int(j,k)=trapz(gyf(j:NY-1),thme(j:NY-1,k));
  end
end

for k=2:nk
  dthme_int(:,k)=(thme_int(:,k)-thme_int(:,k-1))/(tii(k)-tii(k-1));
  dwme_int(:,k)=(wme_int(:,k)-wme_int(:,k-1))/(tii(k)-tii(k-1));
end

jstart=10;
jend=NY-1;

for k=1:nk
  epsilon_si_int(k)=trapz(gyf(jstart:jend),epsilon_si(jstart:jend,k));
  epsilon_int(k)=trapz(gyf(jstart:jend),epsilon(jstart:jend,k));
  P_AG_int(k)=trapz(gyf(jstart:jend),P_AG(jstart:jend,k));
  PG_res_int(k)=trapz(gyf(jstart:jend),-wv(jstart:jend,k)*dvgdz);
  PG_sgs_int(k)=trapz(gy(jstart:jend+1),nu_t_sgs(jstart:jend+1,k).*dwdy(jstart:jend+1,k)*dvgdz);
  u_int(k)=trapz(gyf(jstart:jend),ume_int(jstart:jend,k));
  w_int(k)=trapz(gyf(jstart:jend),wme_int(jstart:jend,k));
  thme_int(k)=trapz(gyf(jstart:jend),thme_int(jstart:jend,k));
  bw_int(k)=RI(1)*trapz(gy(jstart-1:jend+1),thv(jstart-1:jend+1,k));
  bw_sgs_int(k)=RI(1)*nu_t_sgs(jend+1,k)*(thme(jend+1,k)-thme(jend,k))/(gyf(jend+1)-gyf(jend))-RI(1)*nu_t_sgs(jstart-1,k)*(thme(jstart-1,k)-thme(jstart-2,k))/(gyf(jstart-1)-gyf(jstart-2));
  tke_int(k)=trapz(gyf(jstart-1:jend),tke(jstart-1:jend,k));
end

for k=2:nk
  du_int(k)=RI(1)*drhodx/f*(u_int(k)-u_int(k-1))/(tii(k)-tii(k-1));
  dw_int(k)=RI(1)*drhodx/f*(w_int(k)-w_int(k-1))/(tii(k)-tii(k-1));
  db_int(k)=RI(1)*(thme_int(k)-thme_int(k-1))/(tii(k)-tii(k-1));
  dtke_int(k)=(tke_int(k)-tke_int(k-1))/(tii(k)-tii(k-1));
end

for k=2:nk
  dbdt(:,k)=RI(1)*(thme(:,k)-thme(:,k-1))/(tii(k)-tii(k-1));
  dwdt(:,k)=(wme(:,k)-wme(:,k-1))/(tii(k)-tii(k-1));
end

figure
plot((tii(2:nk)-tii(1))/(2*pi/f),epsilon_int(2:nk),'b-')     
hold on
plot((tii(2:nk)-tii(1))/(2*pi/f),epsilon_si_int(2:nk),'c-')     
plot((tii(2:nk)-tii(1))/(2*pi/f),PG_res_int(2:nk)+PG_sgs_int(2:nk),'g-')     
plot((tii(2:nk)-tii(1))/(2*pi/f),P_AG_int(2:nk),'m-')     
plot((tii(2:nk)-tii(1))/(2*pi/f),bw_int(2:nk),'r-');
plot((tii(2:nk)-tii(1))/(2*pi/f),dtke_int(2:nk),'k-');
plot((tii(2:nk)-tii(1))/(2*pi/f),epsilon_int(2:nk)+PG_res_int(2:nk)+PG_sgs_int(2:nk)+P_AG_int(2:nk)+bw_int(2:nk)-dtke_int(2:nk),'k.-');
legend('\int \epsilon dz','\int \epsilon_{SI} dz','\int GSP dz','\int AGSP dz','\int <b''w''> dz','\int -d<k>/dt dz','SUM');
line([0 2.5],[0 0]);   
set(gca,'FontName','Times');
set(gca,'FontSize',14);
xlabel('t/\tau_{inertial}');
title('Depth-integrated energy budget [m^3/s^3]');

figure
plot((tii(2:nk)-tii(1))/(2*pi/f),PG_res_int(2:nk)+PG_sgs_int(2:nk),'g-');
hold on
plot((tii(2:nk)-tii(1))/(2*pi/f),bw_int(2:nk)+bw_sgs_int(2:nk),'r-');
plot((tii(2:nk)-tii(1))/(2*pi/f),-1*ones(nk-1,1)*4.39453e-7*(gyf(end)-gyf(2)),'c-');
plot((tii(2:nk)-tii(1))/(2*pi/f),+dw_int(2:nk),'b-');
plot((tii(2:nk)-tii(1))/(2*pi/f),-db_int(2:nk),'g-');
plot((tii(2:nk)-tii(1))/(2*pi/f),PG_res_int(2:nk)+PG_sgs_int(2:nk)+bw_int(2:nk)+bw_sgs_int(2:nk)-4.39453e-7*(gyf(end)-gyf(2))+dw_int(2:nk)-db_int(2:nk),'k.-');
set(gca,'FontName','Times');
set(gca,'FontSize',14);
 xlabel('t/\tau_{inertial}');
title('Depth-integrated momentum,buoyancy balance [m^3/s^3]');
line([0 2.5],[0 0]);
legend('\int GSP dz','\int <b''w''> dz','H*WDBF','-\int d<b>/dt dz','M^2/f \int d<v>/dt dz','SUM');


figure
plot((tii(2:nk)-tii(1))/(2*pi/f),db_int(2:nk),'r-');
hold on
plot((tii(2:nk)-tii(1))/(2*pi/f),u_int(2:nk)*RI(1)*drhodx,'g-');
plot((tii(2:nk)-tii(1))/(2*pi/f),RI(1)*thv_int,'g-');
plot((tii(2:nk)-tii(1))/(2*pi/f),0.5*(-RI(1)*(nu_t_sgs(jstart+1,2:nk)+NU/PR).*dthdy(jstart+1,2:nk)-RI(1)*(nu_t_sgs(jstart,2:nk)+NU/PR).*dthdy(jstart,2:nk)),'c-');
plot((tii(2:nk)-tii(1))/(2*pi/f),0.5*(-RI(1)*(nu_t_sgs(jend+1,2:nk)+NU/PR).*dthdy(jend+1,2:nk)-RI(1)*(nu_t_sgs(jend,2:nk)+NU/PR).*dthdy(jend,2:nk)),'m-');




% LES_inertial_analysis.m

% load les_front

% to=64.5;
% te=yearday(end);
% clear yearday

% 
% load /home/jacob/LATMIX/matlab/LES_COARE35_forcing.mat
% 
% Ig=find(yearday>=to & yearday<=te);
% 
% t=(yearday(Ig)-yearday(Ig(1)))*86400;




f=9.3e-5;
rho_o=1024;
%f=1.21e-4;
Ti=2*pi/f;
dt=Ti/60;
ti=t(1):dt:t(end);

%interpolate stress to uniform time grid

tau_xi=interp1(t,tx.*rho_o,ti);
tau_yi=interp1(t,ty.*rho_o,ti);


%calculate centered difference derivative of tau_x
Ip=3:length(ti);
Im=1:length(ti)-2;
I=2:length(ti)-1;

tmp=(tau_xi(Ip)-tau_xi(Im))/(2*dt);
dtauxdt=zeros(size(tau_xi));
dtauxdt(I)=tmp;
dtauxdt(1)=dtauxdt(2);
dtauxdt(end)=dtauxdt(end-1);

%calculate depth integrate x-velocity

G=dtauxdt/rho_o+f*tau_yi/rho_o;
U=zeros(size(ti));
for ii=2:length(ti)
    tt=ti(ii);
    F=sin(f*(tt-ti))/f;
   U(ii)=trapz(ti(1:ii),G(1:ii).*F(1:ii)); 
end
U=U+tau_xi(1)/rho_o/f*sin(f*ti);

tmp=(U(Ip)-U(Im))/(2*dt);
V=zeros(size(U));
V(I)=tmp/f;
V(1)=V(2);
V(end)=V(end-1);
V=V-tau_xi/rho_o/f;

U = interp1(ti, U, t)./25;
% DUDZ = U./90;
%%
plot(tm, cumtrapz(t, U), 'linewidth', 2);
hold on;
plot(tm, xposcomunwrap(1:length(tm)), 'linewidth', 2);
hold off
legend('$\int_t U_{SLAB}$', '$X_{COM}$');
xlabel('t'); ylabel('x');
%%

ydi=interp1(t,yearday(Ig),ti);
H=90;

%load matlab figure


hgload('ageo_vel_timeseries.fig');
subplot(2,1,1)
hold on
plot(ydi,U/H,'r','linewidth',2)
set(gca,'fontsize',14);
ylabel('u (m s^{-1})','fontsize',14);
grid on
subplot(2,1,2)
hold on
plot(ydi,V/H,'r','linewidth',2)
set(gca,'fontsize',14);
ylabel('v (m s^{-1})','fontsize',14);
xlabel('yearday','fontsize',14);
grid on

%print -depsc LES_inertial_analysis_COARE.eps
% Calculate diffusivities from second moments
dyevar = dye1;
thetai = x./max(x).*2*pi; % Convert to angle
xi = cos(thetai);
zi = sin(thetai);
ubar = trapz(z(1:end-1), u(:, 1:end-1, :), 2)./abs(z(1));
for i =1:1:length(dye1) % Calculate at each time step.
    
   zeta(i,:) = squeeze(trapz(z(1:end-1), dyevar(:,1:end-1,i), 2)).';
      TT(i) = trapz(x, zeta(i,:)); % Total Tracer
    zxbar(i,:) = squeeze(trapz(x, dyevar(:,:,i),1));
    
%    ZCM(i) = trapz(x, x.*zeta(i,:))./TT(i); % Zonal Center Of mass
   VCM(i) = trapz(z, z.*zxbar(i,:))./TT(i); % Vertical Center Of Mass
%    bpi(i,:) = trapz(z(1:end-1), bp(:,1:end-1,i), 2)./abs(z(1)); 
   
   % Circular center of mass
   xii = xi.*zeta(i,:);
   zii = zi.*zeta(i,:);
   tbar = atan2(-mean(zii), -mean(xii)) + pi;
   ZCM(i) = max(x).*tbar./(2*pi);
   
   % Find index of zonal center of mass
   ind = find(x>=ZCM(i), 1, 'first');
   indshift = ind -length(x)./2;
   if ~isempty(ind)
   xposcom(i) = x(ind);
   end
   indvert = find(z>=VCM(i), 1, 'first'); %Index of vertical COM
   
   % Calculate total b at xposcom
   bpcom(i) = bp(ind, indvert, i);
   
%    ind = 0;
   zetas(i,:) = circshift(zeta(i,:), -indshift, 2); % This properly centers the dye (note not applicable in buoyancy space).
   xs(i,:) = circshift(x, -indshift, 2);
   ubars(:,i) = circshift(ubar(:,:,i), -indshift, 1);
   ddxdye1s(:,:,i) = circshift(ddxdye1(:,:,i), -indshift, 1);
   M20 = trapz(x, x.^2.*zetas(i,:)); % X Pos around Center of Mass
%    M20 = trapz(x, (x-ZCM(i)).^2.*zeta(i,:));
  
%    ZCM(i) = x(length(x)./2);
%    XVar(i) = (M20 - 0*ZCM(i).^2)./TT(i);
    
    XVar(i) = M20./TT(i) - ZCM(i).^2;
   XVar(i) = M20./TT(i) - x(length(x)./2); % Variance around center of x (which is appropriate for centered zetas).
   

end
xposcomunwrap = unwrap(xposcom.*2*pi./5000 - pi).*5000/(2*pi); % Unwrap to a constant track.
%%
ucom = gradient(xposcomunwrap(1:length(t)), t);
for i=1:length(t)
%    GT(i) = trapz(x, (x-ZCM(i)).*zeta(i,:).*(ubar(:,i)-ucom(i)).')./TT(i);
%    GT(i) = trapz(x, (x-ZCM(i)).^2.*squeeze(trapz(z(1:end-1), (u(:,1:end-1,i)-repmat(ucom(i), [nx nz-1 1])).*ddxdye1(:,1:end-1,i), 2)).')./TT(i);
    GT(i) = trapz(xs(i,:), xs(i,:).^2.*squeeze(ubars(:,i)-ucom(i)).'.*squeeze(trapz(z(1:end-1), ddxdye1s(:,1:end-1,i), 2)).')./TT(i);
end

%Total buoyancy field following the center of mass
bcom = bpcom + xposcomunwrap.*dbdx;

% Assuming no strain (ad hoc)
kappa = 1/2.*gradient(smooth(XVar(1:length(t)), 1), t);
%%
[dudz, ~, ~] = gradient(u, z, 1e-11, 1e-11);
alpha = squeeze(nanmean(nanmean(dudz(:,1:end-2,:))));
factor = 1/2.*(alpha./9.2e-5).^2;
%%
% figure
% subplot(3,1,1)
% plot(t, xposcom(1:length(t)));
% subplot(3,1,2)
% plot(xposcomunwrap);
% subplot(3,1,3)
% plot(t./86400,bcom(1:length(t)))
% hold on
% % plot(t./86400,bpcom(1:length(t)))
% hold off
% %%
% figure
% pcolor(tm, x, zetas(1:length(t), :).'); shading interp
%%
% gamma = mean(dudx);
% kappas = 1/2.*gradient(XVar(1:length(t)), t) - gamma(1:length(t)).*XVar(1:length(t));
% 
% kappas = 1/2.*gradient(XVar(1:length(t)), t) +  1/2.*GT(1:length(t)); - 1./2.*XVar(1:length(t)).*TT(1:length(t)).*gradient(1./TT(1:length(t)), t);

% %%
% figure
% subplot(4,1,1)
% plot(gradient(1./TT(1:length(t)), t));
% subplot(4,1,2)
% plot(t, ZCM(1:length(t)));
% subplot(4,1,3)
% plot(t, XVar(1:length(t)))
% subplot(4,1,4)
% % semilogy(t, kappa, 'o')
% semilogy(t, kappas, 'o');
% hold on; 
% semilogy(t, kappa, 'x');
% %  plot(t,-GT(1:length(t)), '--'); hold off
% semilogy(t, kxt, 'd')
% % kxt = kxa; kxt(abs(kxt)>1e4) = NaN;
% grid on
% % hold on; plot(t, kxt, 'r'); hold off
%%
tm = t./86400;
gap = [0.03 0.01]; margh = 0.1; margw = 0.1;
figure
subtightplot(3,1,1, gap, margh, margw)
plot(tm, abs(tx+1i.*ty), 'LineWidth', 2);
% plot(tm, alpha(end-length(t)+1:end))
% plot(tm, cumtrapz(t, squeeze(nanmean(u(:,1:10:end,1:length(t)))), 2))
% hold on; plot(tm, ucom); hold off
ylabel('$|\tau/\rho|$');
grid on
set(gca, 'FontSize', 18)
set(gca, 'xlim', [tm(1) 2]);
set(gca, 'xticklabels', []);

subtightplot(3,1,1, gap, margh, margw)
plot(tm, bcom(1:length(t))./bcom(1), 'LineWidth', 2);
ylabel('$b_{COM}/b_o$');
hold on
plot(tmnf, bcomnf(1:length(tmnf))./bcomnf(1), 'LineWidth', 2)
plot(tmntw, bcomntw(1:length(tmntw))./bcomntw(1), 'LineWidth', 2)
hold off
grid on
set(gca, 'FontSize', 18)
set(gca, 'xlim', [tm(1) 2]);
set(gca, 'xticklabels', []);
legend('Front', 'No Front', 'No SI', 'Location', 'NorthEast')

subtightplot(3,1,2, gap, margh, margw)
plot(tm, XVar(1:length(t))./XVar(1), 'LineWidth', 2);
ylabel('$\sigma^2_{x}/\sigma^2_o$')
grid on
set(gca, 'FontSize', 18)
hold on
plot(tmnf, XVarnf(1:length(tmnf))./XVarnf(1), 'LineWidth', 2)
plot(tmntw, XVarntw(1:length(tmntw))./XVarntw(1), 'LineWidth', 2)
% tm2 = tm.^3.*4e6; tm2 = tm2 - tm2(1) + XVar(1);
% plot(tm, tm2)
hold off
set(gca, 'xlim', [tm(1) 2]);
set(gca, 'xticklabels', []);

subtightplot(3,1,3, gap, margh, margw)
plot(tm, kappa, 'o')
hold on
plot(tmnf, kappanf, 'x');
plot(tmntw, kappantw, '+')
hold off
set(gca, 'ylim', [0 1].*300);
% hold on; plot(tm, kx); hold off
hold on; plot(tm, factor(1:length(tm)).*0.06); hold off;
xlabel('t  (days)');
ylabel('$\kappa$ $(m^2 s^{-1})$');
set(gca, 'FontSize', 18)
grid on
set(gca, 'xlim', [tm(1) 2]);

set(gcf, 'Color','w', 'Position', [    670   255   998   716]);
%%
% kv = [kappa.' kxa];
% kv(abs(kv)>1e4) = NaN;
% hist(kv, 100)
% %%
% scatter(kappa, kappas)
% set(gca, 'xlim', [-1 1].*1e3, 'ylim', [-1 1].*1e3)
% %
% plot(squeeze(nanmean(kx(20:end-1,:)))./1000);
% 
% plot(kxa, 'x');
% hold on
% plot(kappas, 'o');
% hold off
% set(gca, 'ylim', [-1 1].*1000)
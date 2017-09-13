for i=1:length(latmix.t)
    subplot(1,2,1)
pcolor(latmix.x, latmix.y, squeeze(latmix.dye1(:,:,i)).'); shading interp
hold on
contour(latmix.x, latmix.y, squeeze(latmix.b(:,:,i)).', linspace(-.01, .01, 10),'k');
hold off
colorbar;
set(gca, 'clim', [0 1]);
title(num2str(latmix.t(i)./86400))
    subplot(1,2,2)
pcolor(latmix.x, latmix.y, squeeze(latmix.dye2(:,:,i)).'); shading interp
hold on
contour(latmix.x, latmix.y, squeeze(latmix.b(:,:,i)).', linspace(-.01, .01, 10),'k');
hold off
colorbar;
set(gca, 'clim', [0 1]);
title(num2str(latmix.t(i)./86400))

drawnow
% pause();
    
end

%%
ntp = 16;
for i=1:ntp
    
    ti = 1+ floor(length(latmix.t)./(2*ntp))*(i-1);
subplot(4,4,i)
pcolor(latmix.x, latmix.y, real(squeeze(log10(latmix.dye2(:,:,ti)))).'); shading interp
hold on
contour(latmix.x, latmix.y, squeeze(latmix.b(:,:,ti)).', linspace(-.01, .01, 10),'k');
hold off
colorbar;
set(gca, 'clim', [-4 -2]);
title(num2str(latmix.t(ti)./86400+64.5))

end

%%
xvec = latmix.x;
bback = -dbdx.*xvec;
bback = repmat(bback.', [1 100 1420]);
btotals = latmix.b; + bback;

gap =[0.05 0.05]; margh = .1; margw=.2;

subtightplot(3,1,1, gap, margh, margw)
cl = [-2 -0];
ti = 1;
tc = 'w';
nc = 25;
pcolor(latmix.x, latmix.y, real(squeeze(log10(latmix.dye2(:,:,ti)))).'); shading interp
hold on
contour(latmix.x, latmix.y, squeeze(btotals(:,:,ti)).', linspace(-.01, .01, nc),tc);
hold off
cb = colorbar;
set(gca, 'clim', cl);
title('LES of LATMIX Dye Release')
set(gca, 'FontSize', 18);
set(get(cb, 'Ylabel'), 'String', '$\mathrm{log_{10}(Dye)}$', 'Interpreter', 'Latex');
ylabel({['\makebox[1in][c]{Day: ', num2str(latmix.t(ti)./86400+64.5,3),'}'], '\makebox[1in][c]{z (m)}'}, 'Rotation', 0);
% xlabel('km');
set(gca, 'XTickLabel', '')

subtightplot(3,1,2, gap, margh, margw)
cl = [-3 -1];
ti = 90;
pcolor(latmix.x, latmix.y, real(squeeze(log10(latmix.dye2(:,:,ti)))).'); shading interp
hold on
contour(latmix.x, latmix.y, squeeze(btotals(:,:,ti)).', linspace(-.01, .01, nc),tc);
hold off
cb = colorbar;
set(gca, 'clim', cl);
% title(num2str(latmix.t(ti)./86400+64.5))
set(gca, 'FontSize', 18);
set(get(cb, 'Ylabel'), 'String', '$\mathrm{log_{10}(Dye)}$', 'Interpreter', 'Latex');
ylabel({['\makebox[1in][c]{Day: ', num2str(latmix.t(ti)./86400+64.5,3),'}'], '\makebox[1in][c]{z (m)}'}, 'Rotation', 0);
set(gca, 'XTickLabel', '')
subtightplot(3,1,3, gap, margh, margw)
cl = [-5 -2];
ti = 450;
pcolor(latmix.x, latmix.y, real(squeeze(log10(latmix.dye2(:,:,ti)))).'); shading interp
hold on
contour(latmix.x, latmix.y, squeeze(btotals(:,:,ti)).', linspace(-.01, .01, nc),tc);
hold off
cb = colorbar;
set(gca, 'clim', cl);
% title(num2str(latmix.t(ti)./86400+64.5))
set(gca, 'FontSize', 18);
set(get(cb, 'Ylabel'), 'String', '$\mathrm{log_{10}(Dye)}$', 'Interpreter', 'Latex');
ylabel({['\makebox[1in][c]{Day: ', num2str(latmix.t(ti)./86400+64.5,3),'}'], '\makebox[1in][c]{z (m)}'}, 'Rotation', 0);
xlabel('km');
set(gca, 'XTickLabel', '')
set(gca, 'XTickLabel', 1:5)

colormap(cptcmap('Deep_Sea.cpt'))

set(gcf, 'Color', 'w');
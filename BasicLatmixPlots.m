counter = 1;
for i=1:50:length(t)
    subplot(2,2,1:2)
    plot(t, tx);
    hold on
    plot(t, ty);
    yt = get(gca, 'YTick');
    plot(t(i).*ones(size(yt)), yt);
    hold off
    
    subplot(2,2,3)
pcolor(x, z, squeeze(dye1(:,:,i)).'); shading interp
hold on
contour(x, z, squeeze(b(:,:,i)).', linspace(-.01, .01, 25), 'k');
hold off
colorbar;
% set(gca, 'clim', [0 1]);
title(num2str(t(i)./86400))
    subplot(2,2,4)
pcolor(x, z, squeeze(dye2(:,:,i)).'); shading interp
hold on
contour(x, z, squeeze(b(:,:,i)).', linspace(-.01, .01, 25),'k');
hold off
colorbar;
% set(gca, 'clim', [0 1]);
title(num2str(t(i)./86400))
set(gcf, 'Color', 'w', 'Position',[   675   241   968   733]);
drawnow
% pause();
    if counter<21 & false
        saveas(gcf, ['FR_', num2str(i), '.png']);
        counter = counter + 1;
    end
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

%% BASIC PLOT OF  FORCING
gap = [0.05 0.02]; margh=0.2; margw=0.1;

figure
subtightplot(2,1,1, gap, margh, margw)
plot(tm, abs(tx+1i.*ty).*1030, 'LineWidth', 2);
ylabel('$|\tau|$ $(N m^{-2})$');
% xlabel('Days');
set(gca, 'xlim', [tm(1) 2]);
set(gca, 'FontSize', 18)
grid on
set(gca, 'XTickLabels', []);
subtightplot(2,1,2, gap, margh, margw)
plot(tm, Bo, 'LineWidth', 2);
ylabel('$B_o$ $(m s^{-3})$');
xlabel('Days');
set(gca, 'xlim', [tm(1) 2]);
set(gca, 'FontSize', 18)
grid on
set(gcf, 'Color', 'w', 'Position', [ 670   374   970   597]);


%% DYE DISPERSAL OVERVIEW PLOTS
gap = [0.02 0.01]; margh=0.1; margw=0.15;
nc = 25;
figure
conts = linspace(0, 1, 20);

for i=1:4
    if i==1
            tind=1;
    else
            tind=find(tm>=0.4+(0.2).*(i-1), 1, 'first');
%             disp(tind)
%         case i==4
    end
subtightplot(2, 4, (1:2) + 2*(i-1), gap, margh, margw)
dmax = squeeze(max(max(dye1(:,:,tind))));
[c, h]=contourf(x./1000, z, dye1(:,:,tind).'./dmax);shading interp
set(h, 'edgecolor', 'none')
set(gca, 'xTick', 0:1:5);
% colorbar
hold on
grid on
contour(x./1000, z, squeeze(b(:,:,tind)).', linspace(-.01, .01, nc),'w');
hold off
if (i==1 | i==3)
    ylabel('z (m)');
else
    set(gca, 'YTickLabels', []);
end
if i<3
    set(gca, 'XTickLabels', []);
else
    xlabel('km');
end
set(gca, 'clim', [conts(1) conts(end)]);
set(gca, 'FontSize', 16)
end
cb = colorbar;
set(cb, 'TickLabelInterpreter', 'latex')
set(get(cb, 'ylabel'), 'string', {'Normalized concentration', '$\;\;\;\;\;\;$     C/max[C(t)]'}, 'Interpreter', 'Latex','FontSize', 16)
set(gcf, 'Color','w', 'Position',[         322         246        1200         670]);
set(cb, 'Position', [         0.8648    0.2799    0.0195    0.4448]);

%% DYE DISPERSAL OVERVIEW PLOTS LNG
gap = [0.02 0.01]; margh=0.2; margw=0.1;
nc = 25;
figure
conts = linspace(0, 1, 20);

for i=1:4
    if i==1
            tind=1;
    else
            tind=find(tm>=0.4+(0.2).*(i-1), 1, 'first');
%             disp(tind)
%         case i==4
    end
subtightplot(1, 8, (1:2) + 2*(i-1), gap, margh, margw)
dmax = squeeze(max(max(dye1(:,:,tind))));
[c, h]=contourf(x./1000, z, dye1(:,:,tind).'./dmax);shading interp
set(h, 'edgecolor', 'none')
set(gca, 'xTick', 0:1:5);
% colorbar
hold on
grid on
contour(x./1000, z, squeeze(b(:,:,tind)).', linspace(-.01, .01, nc),'w');
hold off
if (i==1)
    ylabel('z (m)');
else
    set(gca, 'YTickLabels', []);
end
% if i<3
%     set(gca, 'XTickLabels', []);
% else
    xlabel('km');
% end
set(gca, 'clim', [conts(1) conts(end)]);
set(gca, 'FontSize', 16)
end
cb = colorbar;
set(cb, 'TickLabelInterpreter', 'latex')
set(get(cb, 'ylabel'), 'string', {'Normalized concentration', '$\;\;\;\;\;\;$     C/max[C(t)]'}, 'Interpreter', 'Latex','FontSize', 16)
set(gcf, 'Color','w', 'Position',[           74         426        1847         468]);
set(cb, 'Position', [      0.9057    0.1991    0.0086    0.6019]);
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

%% 
% Density class histograms
edges = linspace(-1e-2, 1e-2, 50);
AvgD1 = NaN(length(edges), length(t));
for i=1:length(t)
    bvec = reshape(squeeze(b(:,:,i)), nx.*nz, 1);
    d1vec = reshape(squeeze(dye2(:,:,i).*repmat(dz, [nx 1])), nx.*nz, 1);

[N,edges,bin] = histcounts(bvec, edges);

m = (bin > 0);
inds = ismember(1:length(edges), bin(m));
temp = accumarray(bin(m),d1vec(m),[],@nansum);

AvgD1(unique(bin(m)),i) = temp(unique(bin(m)));
end

dz = diff(out1.y);
dz = [dz(1) dz];
df = squeeze(nansum(dye2(:,1:end-1,:).*repmat(dz(1:end-1), [nx 1 nt]), 2));

%%
figure
subplot(3,1,1)
 plot(t, abs(tx+1i.*ty));
 xlabel('t');
 title('|\tau|')
% pcolor(t, z, uflux1); shading interp
% plot(t, nanmean(ua(1:end-1,:)))
% hold on
% plot(t, nanmean(wa(1:end-1,:)));
% hold off
% plot(t, nansum(df(:,1:length(t))))
set(gca, 'xlim', [t(1) t(end)])
subplot(3,1,2)
pcolor(t, x, df(:,1:length(t))); 
shading interp
hold on
contour(t, x, squeeze(b(:,end-2,1:length(t))),'k');
hold off
 xlabel('t');
 ylabel('x');
 title('Depth averaged Dye')
 
subplot(3,1,3)
pcolor(t, edges, AvgD1)
shading interp
 xlabel('t');
title('dye by buoyancy class')

 ylabel('b');
 set(gcf, 'Color', 'w', 'Position', [    670   173   815   798]);
 
 %%
 figure
 subplot(2,1,1)
 plot(t, abs(tx+1i.*ty));
 set(gca, 'xlim', [t(1) t(end)]);
 
 subplot(2,1,2)
 plot(t, ua(1:end-1,:).');
 set(gca, 'xlim', [t(1) t(end)]);
 
 %%
 figure
 subplot(2,1,1)
 semilogy(t, kza);
 subplot(2,1,2)
 semilogy(t, kxa);
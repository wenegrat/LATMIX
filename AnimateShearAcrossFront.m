for i=1:length(t)
   pcolor(x, z(1:end-2), squeeze(dudz(:,1:end-2,i)).');shading interp
   cl = get(gca, 'clim');
   cl = [-1 1].*1e-2;
   hold on
   contour(x, z(1:end-2), squeeze(b(:,1:end-2,i)).', 30, 'k');
   contour(x, z(1:end-2), squeeze(dye1(:,1:end-2,i)).', linspace(0.01,1,20),'w');
   set(gca, 'clim', cl);
   hold off
   title(num2str(tm(i)));
   colorbar;
   xlabel('x');ylabel('z');
   drawnow
%    pause();
    
end
% Set the name of the file to read
moviefile=dlmread('../front/run1/movie_x3.txt');

NX=1280; 
NY=100;
N_TH=0;
LX=5000;
LY=150;

% Get x-vector and y-vector with uniform spacing
for i=1:NX
  xvec(i)=LX*i/NX;
end
for j=1:NY-1
  yvec(j)=gyf(j);
end
yvec(NY)=yvec(NY-1)+yvec(NY-1)-yvec(NY-2);

null(1:NX,1:NY)=0;

colormat=get(gca,'ColorOrder');


% Get the number of timesteps
nk=length(moviefile(:,1))/(NX*NY);

count=0;

% zero matrix mat
mat(1:NX,1:NY,1:nk)=0;

for k=1:nk
k
for i=1:NX
for j=1:NY
count=count+1;
  mat(i,j,k)=moviefile(count);%-thme(j,k)-drhodx*xvec(i);
end
end 
end

clear moviefile;

f1=gcf;
clear M

for k=1:nk
surf(xvec,yvec,zeros(NX,NY,1)'-90,mat(:,:,k)','EdgeColor','none');
view(0,90)
caxis([0 0.01]);
axis([0 LX -LY 0]);
shading interp;
colorbar;
xlabel('Cross-front distance (m)');
ylabel('Depth (m)');
title(['Dye concentration, yearday ' num2str(tii(k)/3600/24+64.5,'%1.1f')]);
M(k)=getframe(gcf);
clf
end




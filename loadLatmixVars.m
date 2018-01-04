%% FOR 1 PARTICULAR MODEL RUN
% Code sticks together the 2 run portions.
nofront = false;
path = '/home/jacob/LATMIX/front/run1/';
out1 = loadLESVars(path);

%%
%%%%%%%%%%%%%%%% RUN 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = '/home/jacob/LATMIX/front/run2/';
out2 = loadLESVars(path);

%%
% nofront = true;
if ~nofront
dbdx = -5e-7;
else
    dbdx = 0;
end
x = out1.x;
z = out1.y;
tx = [out1.tx out2.tx];
ty = [out1.ty out2.ty];
dye1 = cat(3, out1.dye1, out2.dye1);
dye2 = cat(3, out1.dye2, out2.dye2);
b = cat(3, out1.b, out2.b); % XX- doesn't include background gradient
[nx nz nt] = size(b);
bp = b - repmat((dbdx.*(x)).', [1 nz nt]);
t = [out1.ttot out2.ttot];
ua = cat(2, out1.ua, out2.ua);
wa = cat(2, out1.wa, out2.wa);
u = cat(3, out1.u, out2.u);
v = cat(3, out1.v, out2.v);
B2temp = out2.Bo; B2temp(1) = out1.Bo(end);
Bo = [out1.Bo B2temp];


uflux1 = cat(2, squeeze(out1.uflux_a(:,:, 2)), squeeze(out2.uflux_a(:,:,2)));
vflux1 = cat(2, squeeze(out1.wflux_a(:,:, 2)), squeeze(out2.wflux_a(:,:,2)));
[ddzdye1, ddxdye1, ~] = gradient(dye1,z, x, 1e11);
[bz, bx, ~] = gradient(b, z, x, 1e11);

% kz = - vflux1(:,1:length(t))./squeeze(nanmean(ddzdye1(:,:,1:length(t))));
% kz = - squeeze(nanmean(vflux1(1:end-2,:))).'./squeeze(nanmean(nanmean((ddzdye1(:,1:end-2,1:1407)))));

% kx = - uflux1(:,1:length(t))./squeeze(nanmean(ddxdye1(:,:,1:length(t))));
% % kx = -squeeze(nanmean(uflux1(:,1:length(t))))./squeeze(nanmean(nanmean(ddxdye1(:,:,1:length(t)))));
% dia = kx.*squeeze(nanmean(bx(:,:,1:length(t)))) + kz.*squeeze(nanmean(bz(:,:,1:length(t))));

% yl = 2:97;
% ns = 10;
% kza = - smooth(nanmean(vflux1(yl,:)).', ns)./smooth(squeeze(nanmean(nanmean(ddzdye1(:,yl,1:1407)))), ns);
% kxa = - smooth(trapz(z(yl), uflux1(yl,:))./abs(z(1)), ns)./smooth(squeeze(trapz(z(yl), nanmean(ddxdye1(:,yl,1:1407)), 2)./abs(z(1))), ns);

%%
SecondMomentCalcs

%%

% No Front
% tmnf = tm;
% bcomnf = bcom;
% XVarnf = XVar;
% kappanf = kappa;

% No TW
% tmntw = tm;
% bcomntw = bcom;
% XVarntw = XVar;
% kappantw = kappa;

%%
% path = '/data/thomas/thomas1/jacob13/LATMIX/front/run2/';
% readmean
% %%
% btot = cat(3, b, mat);
% ttot = [ttot tii];
% 
% %%
% dye1tot = cat(3, dye1, mat);
% 
% %%
% dye2tot = cat(3, dye2, mat);
% %%
% latmix.b = btot;
% latmix.t = ttot;
% latmix.x = xvec;
% latmix.y = yvec;
% latmix.dye1 = dye1tot;
% latmix.dye2 = dye2tot;
% 
% %%
% save('latmixvars.mat', 'latmix', '-v7.3')
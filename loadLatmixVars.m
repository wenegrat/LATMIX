%% FOR 1 PARTICULAR MODEL RUN
% Code sticks together the 2 run portions.

path = '/data/thomas/thomas1/jacob13/LATMIX/front/run1/';
out1 = loadLESVars(path);

%%
%%%%%%%%%%%%%%%% RUN 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = '/data/thomas/thomas1/jacob13/LATMIX/front/run2/';
out2 = loadLESVars(path);


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
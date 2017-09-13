function out = loadLESVars(path)

readmean;

% 1) readmean.m
% Load Variables

%%%%%%%%%%%%%%%% RUN 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
moviefile=dlmread([path, 'movie_x1.txt']); % buoyancy
plot_movie_xy;

%% Buoyancy is movie_x1;
out.x = xvec;
out.y = yvec; % remember y = z!
out.ttot = tii;

out.b = mat;
out.tx = NU./PR.*dudy(NY,:);
out.ty = NU./PR.*dwdy(NY,:);
out.Bo = squeeze(NU./PR.*dthdy(NY,:,1));
out.tracer_a = thme; % Horizontally averaged Tracer Concentrations
out.wflux_a = thv; % Vertical Tracer Fluxes
out.uflux_a = thu; % Across front tracer fluxes
out.vflux_a = thw; % Along front tracer fluxes
out.uw_a = uv; % Vertical Reynolds Stresses
out.uv_a = uw; % Across front Reynolds stresses
out.vw_a = wv; % Along Front Reynolds Stresses

%% Dye 1
moviefile=dlmread([path, 'movie_x2.txt']); % buoyancy
plot_movie_xy;
out.dye1 = mat;

%% Dye 2
moviefile=dlmread([path, 'movie_x3.txt']); % buoyancy
plot_movie_xy;
out.dye2 = mat;


end
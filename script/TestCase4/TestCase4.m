%% TestCase 4 - MAGNETOSTATIC
clear all;
close all;
clc;

%% PRE-PROCESSING
% Elements along axes
nn = 200; % elements

% Air gap
air_gap = 2e-3; % [m]

% Armature
arm_Z_start = 25e-3; % [m]
arm_Z_end = 33e-3; % [m]
arm_R_start = 0; % [m]
arm_R_end = 30e-3; % [m]

% Coil
coil_Z_start = 8e-3; % [m]
coil_Z_end = 23e-3; % [m]
coil_R_start = 15e-3; % [m]
coil_R_end = 25e-3; % [m]

% Stator
stat_Z_start = 0; % [m]
stat_Z_end = 23e-3; % [m]
stat_R_start = 0; % [m]
stat_R_end = 30e-3; % [m]

% Bounding-box
bound_Z_start = -8.5e-3; % [m]
bound_Z_end = 41.5e-3; % [m]
bound_R_start = 0; % [m]
bound_R_end = 50e-3; % [m]

bound_R = linspace(bound_R_start, bound_R_end, nn);
bound_Z = linspace(bound_Z_start, bound_Z_end, nn);

% Meshgrid combination x and y
[R,Z] = meshgrid(bound_R, bound_Z);

% Coordinate matrix
coo = [R(:), Z(:)];
% Connectivity matrix
conn = delaunay(coo(:,1), coo(:,2));

% Number of nodes
Nn = size(coo,1);
% Number of elements
Ne = size(conn,1);

% Bound nodes
bound_nodes = convhull(coo(:,1), coo(:,2))';

% Triangle object
Triangle = triangulation(conn, coo);
% Compute center coordinate
center = incenter(Triangle);

% Find elements of internal medium using the baricenter points
% Armature
armature_elems = find(center(:,1) >= arm_R_start & center(:,1) <= arm_R_end & ...
    center(:,2) >= arm_Z_start & center(:,2) <= arm_Z_end);

% Coil
coil_elems = find(center(:,1) >= coil_R_start & center(:,1) <= coil_R_end & ...
    center(:,2) >= coil_Z_start & center(:,2) <= coil_Z_end);

% Stator
stator_elems = find(center(:,1) >= stat_R_start & center(:,1) <= stat_R_end & ...
    center(:,2) >= stat_Z_start & center(:,2) <= stat_Z_end);
stator_elems = setdiff(stator_elems, coil_elems);

% Bounding-box
bound_elems = setdiff(1:Ne, [armature_elems; coil_elems; stator_elems])';

% Permeability vacuum
mu_0 = pi*4e-7; % [H/m]
% Armature and stator permeability
mu_core = 2000*mu_0;
% Assign different permeabilty to each sub-domain
mu = ones(Ne,1) .* mu_0;
% Stator
mu(stator_elems) = mu_core;
% Armature
mu(armature_elems) = mu_core;

% Source
% Number of turns
N_coil = 2000; % [turns]
% Current of coil
I_coil = 1; % [A]

% Reluctance path
% Armature
arm1_l = 4e-3; % [m]
arm1_w = 30e-3; % [m]
arm2_l = 20e-3; % [m]
arm2_w = 8e-3; % [m]
arm3_l = arm1_l; % [m]
arm3_w = 30e-3; % [m]

% Air-gap
g_l = air_gap; % [m]

% Stator
stat1_l = 19e-3; % [m]
stat2_l = 20e-3; % [m]
stat2_w = 8e-3; % [m]
stat3_l = stat1_l; % [m]

% Cross section
stat_rad_start = 15e-3; % [m]
coil_rad = 25e-3; % [m]
stat_rad_end = 30e-3; % [m]

% Inner (z-axis) air-gap cross section
g1_s = pi * stat_rad_start^2;
% Farther (z-axis) air-gap cross section
g2_s = pi*stat_rad_end^2 - pi*coil_rad^2;
% Stator and armature horizontal cross section
horiz_s = pi*(stat_rad_start + coil_rad);

%% PROCESSING
% FEM SECTION
% Scaling parameters to r factor (distance of element centroid from the z-axis)
r = center(:,1);
% Vector of material's different reluctivity
vi = (1./mu)./r;

% Source
% Compute value of Jz
Jz = (N_coil*I_coil) / ((coil_Z_end - coil_Z_start) * (coil_R_end - coil_R_start));
forcing_func = zeros(Ne,1);
forcing_func(coil_elems) = Jz;

% Compute the LHS and RHS global matices with fast method
tic;
[LHS_global, RHS_global] = f_FEM_global_fast(conn, coo, vi, forcing_func);
fprintf('Computed time of global matrices with fast method: %.3f sec\n', toc);
fprintf('\n');

% Dirichlet condition exploiting
free_vars = setdiff(1:Nn, bound_nodes);
u = zeros(Nn,1);
% u(bound_nodes) = 0;
% RHS_global = RHS_global - LHS_global*u;
% Solution vector
u(free_vars) = LHS_global(free_vars,free_vars) \ RHS_global(free_vars);

% RELUCTANCE METHOD
% Armature reluctance
armature_rel = (1/mu_core*arm1_l) / (pi*arm1_w^2) + ...
    (1/mu_core*arm2_l) / (horiz_s*arm2_w) + ...
    (1/mu_core*arm3_l) / (pi*arm3_w^2);

% Stator reluctance
stator_rel = ((1/mu_core)*stat1_l) / g2_s + ...
    ((1/mu_core)*stat2_l) / (horiz_s*stat2_w) + ...
    ((1/mu_core)*stat3_l) / g1_s;

% Air-gap reluctance
airgap_rel = (1/mu_0*g_l)/g1_s + (1/mu_0*g_l)/g2_s;

% Flux through the tube flux of reluctance path [Wb]
flux = (N_coil*I_coil) / (armature_rel + stator_rel + airgap_rel);

F_1 = -0.5*(1/mu_0)*(flux/g1_s)^2*g1_s;
F_2 = -0.5*(1/mu_0)*(flux/g2_s)^2*g2_s;

fprintf('Force computed with reluctance method:\n');
fprintf('F inner Z-axis: %.3f N\n', F_1);
fprintf('F outer Z-axis: %.3f N\n', F_2);
fprintf('\n');

%% POST-PROCESSING
% Plots section
% Mesh plot
figure('Name', 'Mesh');
triplot(conn, R, Z, 'c', 'DisplayName', 'Mesh triangle');
title('Mesh figure');
hold on;
% Armature
armature = [[arm_R_start arm_R_end arm_R_end arm_R_start arm_R_start]', ...
   [arm_Z_start arm_Z_start arm_Z_end arm_Z_end arm_Z_start]'];
plot(armature(:,1), armature(:,2), 'g', 'LineWidth', 3, 'DisplayName', 'Armature');
%scatter(center(armature_elems,1), center(armature_elems,2), '*', 'g', 'HandleVisibility','off');
% Stator
stator = [[stat_R_start stat_R_end stat_R_end stat_R_start stat_R_start]', ...
   [stat_Z_start stat_Z_start stat_Z_end stat_Z_end stat_Z_start]'];
plot(stator(:,1), stator(:,2), 'k', 'LineWidth', 3, 'DisplayName', 'Stator');
%scatter(center(stator_elems,1), center(stator_elems,2), '*', 'k', 'HandleVisibility','off');
% Coil
coil = [[coil_R_start coil_R_end coil_R_end coil_R_start coil_R_start]', ...
   [coil_Z_start coil_Z_start coil_Z_end coil_Z_end coil_Z_start]'];
plot(coil(:,1), coil(:,2), 'r', 'LineWidth', 3, 'DisplayName', 'Coil');
%scatter(center(coil_elems,1), center(coil_elems,2), '*', 'r', 'HandleVisibility','off');
% Bounding-box
plot(coo(bound_nodes,1), coo(bound_nodes,2), 'm', 'LineWidth', 3, 'DisplayName', 'Bounds');
legend;
hold off;

% Magnetic vector potential A'_theta plot
figure('Name', "Magnetic vector potential A'-theta [Wb]");
hold on;
% Geometry
plot(armature(:,1), armature(:,2), 'g', 'LineWidth', 3, 'DisplayName', 'Armature');
plot(stator(:,1), stator(:,2), 'k', 'LineWidth', 3, 'DisplayName', 'Stator');
plot(coil(:,1), coil(:,2), 'r', 'LineWidth', 3, 'DisplayName', 'Coil');
plot(coo(bound_nodes,1), coo(bound_nodes,2), 'm', 'LineWidth', 3, 'DisplayName', 'Bounds');
% Contour potential plot
contour(R, Z, reshape(u,[nn,nn]), 'DisplayName', "A'-theta [Wb]");
title("Magnetic vector potential A'-theta [Wb]");
xlabel('R-axis [m]');
ylabel('Z-axis [m]');
% legend('Potential [V]', 'Top plate', 'Bottom plate', 'Internal medium');
axis equal;
colorbar;
legend;
hold off;

% Magnetic flux density field B plot
% Magnetic flux density components in element's nodes
[B_r, B_z] = f_Bfield_axisymm_interp_T3(coo, coo, conn, u);
B_norm = sqrt(B_r.^2 + B_z.^2);

% Magnetic flux density B contour plot
figure('Name', 'Magnetic flux density B [T]');
hold on;
% Magnetic flux density B field magnitude plot
imagesc(bound_R, bound_Z, reshape(B_norm, nn, nn), 'CDataMapping', 'scaled', 'Interpolation', 'bilinear');
% Geometry
plot(armature(:,1), armature(:,2), 'g', 'LineWidth', 3, 'DisplayName', 'Armature');
plot(stator(:,1), stator(:,2), 'k', 'LineWidth', 3, 'DisplayName', 'Stator');
plot(coil(:,1), coil(:,2), 'r', 'LineWidth', 3, 'DisplayName', 'Coil');
plot(coo(bound_nodes,1), coo(bound_nodes,2), 'm', 'LineWidth', 3, 'DisplayName', 'Bounds');
% Contour potential plot with A'_theta
contour(R, Z, reshape(u,[nn,nn]), 'w', 'DisplayName', "A'-theta");
title('Magnetic flux density B');
xlabel('R-axis [m]');
ylabel('Z-axis [m]');
legend;
axis equal;
colorbar;
hold off;

% Magnetic flux density B on air-gap line
% Line in the middle of air-gap
middle_ag = stat_Z_end + air_gap/2;
ag_line = [linspace(0, 40e-3, Nn)' ones(Nn,1)*middle_ag];
% Magflux A'_theta
magflux_ag_line = f_pot_interp_T3(ag_line, coo, conn, u);
% Magnetic flux density components
[B_r, B_z] = f_Bfield_axisymm_interp_T3(ag_line, coo, conn, u);

% Import COMSOL output
comsol_A = readtable("./comsol_A'_airgap.txt");
comsol_A = table2array(comsol_A);
comsol_B_norm = readtable('./comsol_Bnorm_airgap.txt');
comsol_B_norm = table2array(comsol_B_norm);
comsol_Br = readtable('./comsol_Br_airgap.txt');
comsol_Br = table2array(comsol_Br);
comsol_Bz = readtable('./comsol_Bz_airgap.txt');
comsol_Bz = table2array(comsol_Bz);

% Magnetic vector potential A'-theta on air-gap line
figure('Name', "Magnetic vector potential A'-theta on air-gap line");
hold on;
plot(ag_line(:,1), magflux_ag_line, 'b', 'LineWidth', 1, 'DisplayName', "A'-theta [Wb]");
plot(comsol_A(:,1), comsol_A(:,2), '--m', 'DisplayName', 'Comsol data');
hold off;
title("Magnetic vector potential A'-theta [Wb] on air-gap line");
xlabel('Air-gap line [m]');
ylabel("A'-theta [Wb]");
legend;

figure('Name', 'Magnetic flux density B on air-gap line');
subplot(3,1,1);
% Norm of B
B_norm = sqrt(B_r.^2 + B_z.^2);
hold on;
plot(ag_line(:,1), B_norm, 'b', 'LineWidth', 1, 'DisplayName', 'Norm of B [T]');
plot(comsol_B_norm(:,1), comsol_B_norm(:,2), '--m', 'DisplayName', 'Comsol data');
hold off;
title('Magnetic flux density norm B [T] on air-gap line');
xlabel('Air-gap line [m]');
ylabel('norm B [T]');
legend;
% Component along r of B
subplot(3,1,2);
hold on;
plot(ag_line(:,1), B_r, 'g', 'DisplayName', 'Component Br [T]');
plot(comsol_Br(:,1), comsol_Br(:,2), '--m', 'DisplayName', 'Comsol data');
hold off;
title('Magnetic flux density Br component [T] on air-gap line');
xlabel('Air-gap line [m]');
ylabel('Br [T]');
legend;
% Component along z of B
subplot(3,1,3);
hold on;
plot(ag_line(:,1), B_z, 'k', 'DisplayName', 'Component Bz [T]');
plot(comsol_Bz(:,1), comsol_Bz(:,2), '--m', 'DisplayName', 'Comsol data');
hold off;
title('Magnetic flux density Bz component [T] on air-gap line');
xlabel('Air-gap line [m]');
ylabel('Bz [T]');
legend;
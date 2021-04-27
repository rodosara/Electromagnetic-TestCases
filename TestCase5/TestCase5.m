%% TestCase 5 - MAGNETODYNAMICS
clear all;
close all;
clc;

%% PRE-PROCESSING
% Elements along axes
nn = 100; % elements

% Core
core_R_start = 0; % [m]
core_R_end = 5e-2; % [m]
core_Z_start = -5e-2; % [m]
core_Z_end = 5e-2; % [m]

% Coil
coil_R_start = 6e-2; % [m]
coil_R_end = 8e-2; % [m]
coil_Z_start = -2e-2; % [m]
coil_Z_end = 2e-2; % [m]

% Bounding-box
bound_R_start = 0; % [m]
bound_R_end = 20e-2; % [m]
bound_Z_start = -15e-2; % [m]
bound_Z_end = 15e-2; % [m]

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
% Core
core_elems = find(center(:,1) >= core_R_start & center(:,1) <= core_R_end & ...
    center(:,2) >= core_Z_start & center(:,2) <= core_Z_end);

% Coil
coil_elems = find(center(:,1) >= coil_R_start & center(:,1) <= coil_R_end & ...
    center(:,2) >= coil_Z_start & center(:,2) <= coil_Z_end);

% Bounding-box
bound_elems = setdiff(1:Ne, [core_elems; coil_elems])';

% Materials parameters
% Permeability vacuum
mu_0 = pi*4e-7; % [H/m]
% Assign all sub-domain mu_0
mu = ones(Ne,1) .* mu_0;
% Conducibility
sigma_core = 37.7e6; % [S/m]
sigma = zeros(Ne,1);
% Core
sigma(core_elems) = sigma_core;

% Frequency
freq = 100; % [Hz]

% Source
I_coil_RMS = 4000; % [A]

%% PROCESSING

% FEM DYNAMIC SECTION
% Scaling parameters to r factor (distance of element centroid from the z-axis)
r = center(:,1);
% Vector of material's different reluctivity rescaled
vi = (1./mu)./r;

% Source
section_coil = (coil_R_end - coil_R_start) * (coil_Z_end - coil_Z_start);
% Current density
Jz = I_coil_RMS / section_coil; % [A/m^2]
forcing_func = zeros(Ne,1);
forcing_func(coil_elems) = Jz;

% Compute the LHS and RHS global matices with fast method
tic;
[LHS_global, RHS_global] = f_FEM_dynamic_fast(conn, coo, vi, forcing_func);
Mass_matrix = f_Mass_matrix(conn, coo, sigma./r);
LHS_global = LHS_global + (1i*2*pi*freq).*Mass_matrix;
fprintf('Computed time of global matrices with fast method: %.3f sec\n', toc);
fprintf('\n');

% Dirichlet condition exploiting
free_vars = setdiff(1:Nn, bound_nodes);
u = zeros(Nn,1);
% Solution vector
u(free_vars) = LHS_global(free_vars,free_vars) \ RHS_global(free_vars);

%% POST-PROCESSING
% Mesh plot
figure('Name', 'Mesh');
triplot(conn, R, Z, 'c', 'DisplayName', 'Triangle mesh');
title('Mesh figure');
hold on;
% Core
core = [[core_R_start core_R_end core_R_end core_R_start core_R_start]', ...
   [core_Z_start core_Z_start core_Z_end core_Z_end core_Z_start]'];
plot(core(:,1), core(:,2), 'g', 'LineWidth', 3, 'DisplayName', 'Core');
scatter(center(core_elems,1), center(core_elems,2), '*', 'g', 'HandleVisibility','off');
% Coil
coil = [[coil_R_start coil_R_end coil_R_end coil_R_start coil_R_start]', ...
   [coil_Z_start coil_Z_start coil_Z_end coil_Z_end coil_Z_start]'];
plot(coil(:,1), coil(:,2), 'r', 'LineWidth', 3, 'DisplayName', 'Coil');
scatter(center(coil_elems,1), center(coil_elems,2), '*', 'r', 'HandleVisibility','off');
% Bounding-box
plot(coo(bound_nodes,1), coo(bound_nodes,2), 'm', 'LineWidth', 3, 'DisplayName', 'Boundary');
xlabel('R-axis [m]');
ylabel('Z-axis [m]');
legend;
hold off;

% Joule losses in core function
joule_losses = f_joule_losses(coo, conn, sigma, freq, u, core_elems);
fprintf('Joule losses in core: %.3f [W]\n', joule_losses);

% Magnetic vector potential A'_theta contour plot
figure('Name', "Magnetic vector potential A'-theta");
% Absolute values
subplot(1,3,1);
hold on;
% Geometry
plot(core(:,1), core(:,2), 'g', 'LineWidth', 3);
plot(coil(:,1), coil(:,2), 'r', 'LineWidth', 3);
plot(coo(bound_nodes,1), coo(bound_nodes,2), 'm', 'LineWidth', 3);
% Contour potential plot
contour(R, Z, reshape(abs(u),[nn,nn]));
title("Magnetic vector potential A'-theta norm [Wb]");
xlabel('R-axis [m]');
ylabel('Z-axis [m]');
legend('Core', 'Coil', 'Boundary', "A'-theta norm [Wb]");
axis equal;
colorbar;
hold off;

% Real part values
subplot(1,3,2);
hold on;
% Geometry
plot(core(:,1), core(:,2), 'g', 'LineWidth', 3);
plot(coil(:,1), coil(:,2), 'r', 'LineWidth', 3);
plot(coo(bound_nodes,1), coo(bound_nodes,2), 'm', 'LineWidth', 3);
% Contour potential plot
contour(R, Z, reshape(real(u),[nn,nn]));
title("Magnetic vector potential A'-theta real values [Wb]");
xlabel('R-axis [m]');
ylabel('Z-axis [m]');
legend('Core', 'Coil', 'Boundary', "Re[A'-theta] [Wb]");
axis equal;
colorbar;
hold off;

% Imaginary part values
subplot(1,3,3);
hold on;
% Geometry
plot(core(:,1), core(:,2), 'g', 'LineWidth', 3);
plot(coil(:,1), coil(:,2), 'r', 'LineWidth', 3);
plot(coo(bound_nodes,1), coo(bound_nodes,2), 'm', 'LineWidth', 3);
% Contour potential plot
contour(R, Z, reshape(imag(u),[nn,nn]));
title("Magnetic vector potential A'-theta imaginary values [Wb]");
xlabel('R-axis [m]');
ylabel('Z-axis [m]');
legend('Core', 'Coil', 'Boundary', "Im[A'-theta] [Wb]");
axis equal;
colorbar;
hold off;

% Magnetic flux density field B contour plot
% Magnetic flux density components in element's nodes
[Br, Bz] = f_Bfield_axisymm_interp_T3(coo, coo, conn, u);

% Magnetic flux density B contour plot
figure('Name', 'Magnetic flux density abs(B)');
% Absolute value
Br_abs = abs(Br);
Bz_abs = abs(Bz);
B_norm = sqrt(Br_abs.^2 + Bz_abs.^2);
subplot(1,3,1);
hold on;
% Magnetic flux density B field magnitude plot
imagesc(bound_R, bound_Z, reshape(abs(B_norm), nn, nn), 'CDataMapping', 'scaled', 'Interpolation', 'bilinear');
% Geometry
plot(core(:,1), core(:,2), 'g', 'LineWidth', 3);
plot(coil(:,1), coil(:,2), 'r', 'LineWidth', 3);
plot(coo(bound_nodes,1), coo(bound_nodes,2), 'm', 'LineWidth', 3);
% Contour potential plot with A'_theta
contour(R, Z, reshape(abs(u),[nn,nn]), 'w');
title('Magnetic flux density B norm [T]');
xlabel('R-axis [m]');
ylabel('Z-axis [m]');
legend('Core', 'Coil', 'Boundary', "A'-theta norm [Wb]");
axis equal;
colorbar;
hold off;

% Real values
Br_real = real(Br);
Bz_real = real(Bz);
B_norm = sqrt(Br_real.^2 + Bz_real.^2);
subplot(1,3,2);
hold on;
% Magnetic flux density B field real part plot
imagesc(bound_R, bound_Z, reshape(B_norm, nn, nn), 'CDataMapping', 'scaled', 'Interpolation', 'bilinear');
% Geometry
plot(core(:,1), core(:,2), 'g', 'LineWidth', 3);
plot(coil(:,1), coil(:,2), 'r', 'LineWidth', 3);
plot(coo(bound_nodes,1), coo(bound_nodes,2), 'm', 'LineWidth', 3);
% Contour potential plot with A'_theta
contour(R, Z, reshape(real(u),[nn,nn]), 'w');
title('Magnetic flux density B real values [T]');
xlabel('R-axis [m]');
ylabel('Z-axis [m]');
legend('Core', 'Coil', 'Boundary', "Re[A'-theta] [Wb]");
axis equal;
colorbar;
hold off;

% Imaginary values
Br_imag = imag(Br);
Bz_imag = imag(Bz);
B_norm = sqrt(Br_imag.^2 + Bz_imag.^2);
subplot(1,3,3);
hold on;
% Magnetic flux density B field magnitude plot
imagesc(bound_R, bound_Z, reshape(B_norm, nn, nn), 'CDataMapping', 'scaled', 'Interpolation', 'bilinear');
% Geometry
plot(core(:,1), core(:,2), 'g', 'LineWidth', 3);
plot(coil(:,1), coil(:,2), 'r', 'LineWidth', 3);
plot(coo(bound_nodes,1), coo(bound_nodes,2), 'm', 'LineWidth', 3);
% Contour potential plot with A'_theta
contour(R, Z, reshape(imag(u),[nn,nn]), 'w');
title('Magnetic flux density B imaginary values [T]');
xlabel('R-axis [m]');
ylabel('Z-axis [m]');
legend('Core', 'Coil', 'Boundary', "Im[A'-theta] [Wb]");
axis equal;
colorbar;
hold off;

% Magnetic flux density B and magnetic vector pot A'-theta on core line
% Line in the middle of core
horz_line = [linspace(0, 5e-2, 1000)' zeros(1000,1)];
magflux_horz_line = f_pot_interp_T3(horz_line, coo, conn, u);
% Magnetic flux density components
[Br, Bz] = f_Bfield_axisymm_interp_T3(horz_line, coo, conn, u);

% Magnetic vector potential A'-theta on horizontal line
% Import COMSOL output [r_comp z_comp norm(Bz)]
comsol_Atheta = readtable("./comsol_Atheta_coreline.txt", 'HeaderLines', 9);
comsol_Atheta = table2array(comsol_Atheta);

figure('Name', "Magnetic vector potential A'-theta on horizontal line");
% Absolute values
subplot(3,1,1);
hold on;
plot(horz_line(:,1), abs(magflux_horz_line), 'b', 'DisplayName', "Magnitude A'-theta [Wb]");
plot(comsol_Atheta(:,1), comsol_Atheta(:,3), '--m', 'DisplayName', 'Comsol data');
hold off;
title("Magnitude A'-theta [Wb] on horizontal line");
xlabel('R-axis [m]');
ylabel("Abs\{A'-theta\} [Wb]");
legend;
% Real values
subplot(3,1,2);
hold on;
plot(horz_line(:,1), real(magflux_horz_line), 'g', 'DisplayName', "Re\{A'-theta\} [Wb]");
plot(comsol_Atheta(:,1), comsol_Atheta(:,4), '--m', 'DisplayName', 'Comsol data');
hold off;
title("Real A'-theta [Wb] on horizontal line");
xlabel('R-axis [m]');
ylabel("Re\{A'-theta\} [Wb]");
legend();
% Imaginary values
subplot(3,1,3);
hold on;
plot(horz_line(:,1), imag(magflux_horz_line), 'k', 'DisplayName', "Imag\{A'-theta\} [Wb]");
plot(comsol_Atheta(:,1), comsol_Atheta(:,5), '--m', 'DisplayName', 'Comsol data');
hold off;
title("Imaginary A'-theta [Wb] on horizontal line");
xlabel('R-axis [m]');
ylabel("Imag\{A'-theta\} [Wb]");
legend;

% Magnetic flux density B on horizontal line
% Import COMSOL output [r_comp z_comp norm(Bz)]
comsol_Bnorm = readtable("./comsol_Bnorm_coreline.txt", 'HeaderLines', 8);
comsol_Bnorm = table2array(comsol_Bnorm);

figure('Name', 'Magnetic flux density B on horizontal line');
% Absolute values
B_norm = sqrt(abs(Br).^2 + abs(Bz).^2);
subplot(3,1,1);
hold on;
plot(horz_line(:,1), abs(B_norm), 'b', 'DisplayName', 'Magnitude Bnorm [T]');
plot(comsol_Bnorm(:,1), comsol_Bnorm(:,2), '--m', 'DisplayName', 'Comsol data');
hold off;
title('Magnitude B [T] on horizontal line');
xlabel('R-axis [m]');
ylabel('Abs\{B\} [T]');
legend;
% Real values
B_norm = sqrt(real(Br).^2 + real(Bz).^2);
subplot(3,1,2);
plot(horz_line(:,1), real(B_norm), 'g');
title('Real B [T] on horizontal line');
xlabel('R-axis [m]');
ylabel('Re\{B\} [T]');
legend('Re\{B\} [T]');
% Imaginary values
B_norm = sqrt(imag(Br).^2 + imag(Bz).^2);
subplot(3,1,3);
plot(horz_line(:,1), imag(B_norm), 'k');
title('Imaginary B [T] on horizontal line');
xlabel('R-axis [m]');
ylabel('Imag\{B\} [T]');
legend('Imag\{B\} [T]');

% Magnetic flux density Br component on horizontal line
% Import COMSOL output [r_comp z_comp abs(Bz) real(Bz) imag(Bz)]
comsol_Br = readtable("./comsol_Br_coreline.txt", 'HeaderLines', 9);
comsol_Br = table2array(comsol_Br);

figure('Name', 'Magnetic flux density Br component on horizontal line');
% Absolute values
subplot(3,1,1);
hold on;
plot(horz_line(:,1), abs(Br), 'b', 'DisplayName', 'Magnitude Br [T]');
plot(comsol_Br(:,1), comsol_Br(:,3), '--m', 'DisplayName', 'Comsol data');
hold off;
title('Magnitude Br component [T] on horizontal line');
xlabel('R-axis [m]');
ylabel('Abs\{Br\} [T]');
legend;
% Real values
subplot(3,1,2);
hold on;
plot(horz_line(:,1), real(Br), 'g', 'DisplayName', 'Re\{Br\} [T]');
plot(comsol_Br(:,1), comsol_Br(:,4), '--m', 'DisplayName', 'Comsol data');
hold off;
title('Real Br component [T] on horizontal line');
xlabel('R-axis [m]');
ylabel('Re\{Br\} [T]');
legend;
% Imaginary values
subplot(3,1,3);
hold on;
plot(horz_line(:,1), imag(Br), 'k', 'DisplayName', 'Imag\{Br\} [T]');
plot(comsol_Br(:,1), comsol_Br(:,5), '--m', 'DisplayName', 'Comsol data');
hold off;
title('Imaginary Br component [T] on horizontal line');
xlabel('R-axis [m]');
ylabel('Imag\{Br\} [T]');
legend;

% Magnetic flux density Bz component on horizontal line
% Import COMSOL output [r_comp z_comp abs(Bz) real(Bz) imag(Bz)]
comsol_Bz = readtable("./comsol_Bz_coreline.txt", 'HeaderLines', 9);
comsol_Bz = table2array(comsol_Bz);

figure('Name', 'Magnetic flux density Bz component on horizontal line');
% Absolute values
subplot(3,1,1);
hold on;
plot(horz_line(:,1), abs(Bz), 'b', 'DisplayName', 'Magnitude Bz [T]');
plot(comsol_Bz(:,1), comsol_Bz(:,3), '--m', 'DisplayName', 'Comsol data');
hold off;
title('Magnitude Bz component [T] on horizontal line');
xlabel('R-axis [m]');
ylabel('Abs\{Bz\} [T]');
legend;
% Real values
subplot(3,1,2);
hold on;
plot(horz_line(:,1), real(Bz), 'g', 'DisplayName', 'Re\{Bz\} [T]');
plot(comsol_Bz(:,1), comsol_Bz(:,4), '--m', 'DisplayName', 'Comsol data');
hold off;
title('Real Bz component [T] on horizontal line');
xlabel('R-axis [m]');
ylabel('Re\{Bz\} [T]');
legend;
% Imaginary values
subplot(3,1,3);
hold on;
plot(horz_line(:,1), imag(Bz), 'k', 'DisplayName', 'Imag\{Bz\} [T]');
plot(comsol_Bz(:,1), comsol_Bz(:,5), '--m', 'DisplayName', 'Comsol data');
hold off;
title('Imaginary Bz component [T] on horizontal line');
xlabel('R-axis [m]');
ylabel('Imag\{Bz\} [T]');
legend;
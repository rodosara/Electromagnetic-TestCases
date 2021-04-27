%% TestCase 6A - HEAT DYNAMICS STEADY-STATE
clear all;
close all;
clc;

%% PRE-PROCESSING
% Elements along axes
nn = 300; % elements

% Billet
billet_R_start = 0; % [m]
billet_R_end = 5e-2; % [m]
billet_Z_start = 0; % [m]
billet_Z_end = 10e-2; % [m]

% Side discretization
r = linspace(billet_R_start, billet_R_end, nn);
z = linspace(billet_Z_start, billet_Z_end, nn);

% Meshgrid combination x and y
[R,Z] = meshgrid(r, z);

% Coordinate matrix
coo = [R(:), Z(:)];
% Connectivity matrix
conn = delaunay(coo(:,1), coo(:,2));

% Number of nodes
Nn = size(coo,1);
% Number of elements
Ne = size(conn,1);

% Triangle object
Triangle = triangulation(conn, coo);
% Compute center coordinate
center = incenter(Triangle);

% Materials parameters
% Mass density
ro = 2700; % [Kg/m^3]
% Specific heat capacity
Cp = 900; % [J/(Kg*K)]
% Thermal conductivity
lambda = 237; % [W/(k*m)]
% Heat generation rate
g = 1e6; % [W/m^3]
% Heat trasfer coeeficient
h = 15; % [W/(k*m^2)];
% Temperature unperturbed bulk
T_inf = 303.15; % [K]

%% PROCESSING
% FEM DYNAMIC SECTION

% Scaling axysymmetric factor
r = center(:,1);

% Thermal conductivity in all domain
lambda = ones(Ne,1)*lambda;
lambda = lambda .*r;

% Forcing function in all domain
forcing_func = ones(Ne,1)*g.*r;

% Compute the LHS and RHS global matices with fast method
tic;
[STIFF_matrix, B_matrix] = f_FEM_dynamic_fast(conn, coo, lambda, forcing_func);

% Impose convective BC
% Compute boundary edge index and elements edge
[bndedg,edg] = f_quickedges(conn);
ibndedg = bndedg;
% Remove symmetry axes from boundary edges
ibndedg(coo(edg(bndedg,1),1) == 0) = [];
[G_lhs, G_rhs] = f_BC_convective(conn, coo, h, T_inf, edg, ibndedg);

% Final RHS and LHS matrices
LHS_matrix = STIFF_matrix + G_lhs;
RHS_matrix = B_matrix + G_rhs;

% Potential computation
u = LHS_matrix \ RHS_matrix;

fprintf('Computed time of global matrices with fast method: %.3f sec\n', toc);
fprintf('\n');

%% POST-PROCESSING
% Mesh plot
figure('Name', 'Mesh');
triplot(conn, R, Z, 'c', 'DisplayName', 'Mesh');
title('Mesh figure');
hold on;
% Billet
billet = [[billet_R_start billet_R_end billet_R_end billet_R_start billet_R_start]', ...
   [billet_Z_start billet_Z_start billet_Z_end billet_Z_end billet_Z_start]'];
plot(billet(:,1), billet(:,2), 'g', 'LineWidth', 3, 'DisplayName', 'Billet');
scatter(coo(edg(ibndedg,1),1), coo(edg(ibndedg,2),2), '*', 'r', 'DisplayName', 'Convective BCs');
hold off;
legend;
xlabel('R-axis [m]');
ylabel('Z-axis [m]');
axis equal;

% Contour temperature plot
figure('Name', 'Contour plot temperature [K]');
contour(R, Z, reshape(u, nn, nn));
title('Contour plot temperature [K]');
% legend('Temperature [T]');
colorbar;
xlabel('R-axis [m]');
ylabel('Z-axis [m]');
axis equal;

% Temperature T in middle of the billet plot
% Line in the middle of the billet
middle_billet = (billet_Z_end - billet_Z_start)/2;
horz_line = [linspace(0, 5e-2, 1000)' ones(1000,1)*middle_billet];
T_horz_line = f_pot_interp_T3(horz_line, coo, conn, u);

% Import COMSOL output
comsol_temp_middle = readtable('./comsol_T_middleline.txt', 'HeaderLines', 8);
comsol_temp_middle = table2array(comsol_temp_middle);

figure('Name', 'Temperature T [K] on horizontal line');
hold on;
plot(horz_line(:,1), T_horz_line, 'b', 'DisplayName', 'Temperature T [K]');
plot(comsol_temp_middle(:,1), comsol_temp_middle(:,2), '--r', 'DisplayName', 'Comsol data');
hold off;
title('Temperature T [K] on horizontal line');
xlabel('R-axis [m]');
ylabel('Temperature T [K]');
legend();
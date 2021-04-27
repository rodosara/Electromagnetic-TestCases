%% TestCase 7 - COUPLED PROBLEM
clear all;
close all;
clc;

%% PRE-PROCESSING
% Elements along axes
nn = 100; % elements (higher number of elements cause problem for f_quickedges)

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

%% Electrical problem mesh
% Coordinate matrix for electrical problem
coo_el = [R(:), Z(:)];
% Connectivity matrix
conn_el = delaunay(coo_el(:,1), coo_el(:,2));

% Bound nodes
bound_nodes = convhull(coo_el(:,1), coo_el(:,2))';

% Triangle object for electrical problem
Triangle_el = triangulation(conn_el, coo_el);
% Compute center coordinate
center = incenter(Triangle_el);

% Number of nodes electrical mesh
Nn_el = size(coo_el,1);
% Number of elements electrical mesh
Ne_el = size(conn_el,1);

% Find elements of sub-domains using the baricenter points
% Core
core_elems = find(center(:,1) >= core_R_start & center(:,1) <= core_R_end & ...
    center(:,2) >= core_Z_start & center(:,2) <= core_Z_end);

% Coil
coil_elems = find(center(:,1) >= coil_R_start & center(:,1) <= coil_R_end & ...
    center(:,2) >= coil_Z_start & center(:,2) <= coil_Z_end);

% Bounding-box
bound_elems = setdiff(1:Ne_el, [core_elems; coil_elems])';

% %% Thermal problem mesh
% % Coil sides
% coil_R = linspace(coil_R_start, coil_R_end, nn);
% coil_Z = linspace(coil_Z_start, coil_Z_end, nn);
% 
% % Meshgrid combination x and y
% [R,Z] = meshgrid(coil_R, coil_Z);
% % Coordinate matrix for thermal problem
% coo_th = [R(:), Z(:)];
% % Connectivity matrix
% conn_th = delaunay(coo_th(:,1), coo_th(:,2));
% 
% % Triangle object for electrical problem
% Triangle_th = triangulation(conn_th, coo_th);
% % Compute center coordinate
% % center_th = incenter(Triangle_th);
% 
% % Number of nodes thermal mesh
% Nn_th_2 = size(coo_th,1);
% % Number of elements thermal mesh
% Ne_th_2 = size(conn_th,1);
% Number of nodes thermal mesh
Nn_th = length(unique(conn_el(core_elems,:)));
% Number of elements thermal mesh
Ne_th = length(conn_el(core_elems,:));

%% Electrical materials parameters
% Permeability vacuum
mu_0 = pi*4e-7; % [H/m]
% Assign all sub-domain mu_0
mu = ones(Ne_el,1) .* mu_0;
% Conducibility
sigma_core = 37.7e6; % [S/m]
% sigma = zeros(Ne_el,1);
% % Core
% sigma(core_elems) = sigma_core;

% Frequency
freq = 100; % [Hz]

% Source
I_coil_RMS = 4000; % [A]

%% Thermal material parameters
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
% Reference temperature
T0 = 293.15; % [K]
% Temperature coefficient
alpha = 4.29e-3;

% Theta method parameters
theta = 0.5;
Tmax = 10e3; % [s]
timestep = 100; % [s]
itnum = round(Tmax/timestep);

%% PROCESSING

%% Electrical problem configuration
% Scaling parameters to r factor (distance of element centroid from the z-axis)
r = center(:,1);
% Vector of material's different reluctivity rescaled
vi = (1./mu)./r;

% Source
section_coil = (coil_R_end - coil_R_start) * (coil_Z_end - coil_Z_start);
% Current density
Jz = I_coil_RMS / section_coil; % [A/m^2]
forcing_func_el = zeros(Ne_el,1);
forcing_func_el(coil_elems) = Jz;

%% Thermal problem configuration
% Diffusion coefficient in all domain
q = ones(Ne_th,1)*ro*Cp;
% Scaling axysymmetric factor
q = q.*r(core_elems);

% Thermal conductivity in all domain
lambda = ones(Ne_th,1)*lambda;
lambda = lambda.*r(core_elems);

% Impose convective BC
% Compute boundary edge index and elements edge
[bndedg,edg] = f_quickedges(conn_el(core_elems,:));
ibndedg = bndedg;
% Remove symmetry axes from boundary edges
ibndedg(coo_el(edg(bndedg,1),1) == 0) = [];

%% Theta method
times = zeros(itnum-1,1);
temp_centre = zeros(itnum-1,1);

%% Initial condition
% Initial temperature
oldsol_th = T_inf*ones(Nn_th,1);

% Initial current density
sigma = sigma_core*ones(Ne_el,1);
% Electrical problem matrix assembling
[LHS_global, RHS_global] = f_FEM_dynamic_fast(conn_el, coo_el, vi, forcing_func_el);
Mass_matrix = f_MASS_matrix(conn_el, coo_el, sigma./r);
LHS_global = LHS_global + (1i*2*pi*freq).*Mass_matrix;
% Dirichlet BC electrical exploiting
free_vars = setdiff(1:Nn_el, bound_nodes);
A = zeros(Nn_el,1);
% Solution vector
A(free_vars) = LHS_global(free_vars,free_vars) \ RHS_global(free_vars);
oldsol_el = A;




% Forcing function in core sub-domain
forcing_func_th = ones(Ne_th,1).*r(core_elems);

tic;
for k = 1:itnum
    times(k) = (k-1)*timestep;

    %% Electrical problem
    % Update the sigma values of each element
    sigma = f_sigma_update(conn_el, coo_el, oldsol_th, sigma_core, alpha, T0, core_elems);
    
    % Electrical problem matrix assembling
    [LHS_global, RHS_global] = f_FEM_dynamic_fast(conn_el, coo_el, vi, forcing_func_el);
    Mass_matrix = f_MASS_matrix(conn_el, coo_el, sigma./r);
    LHS_global = LHS_global + (1i*2*pi*freq).*Mass_matrix;

    %% Thermal problem
    % Compute Joule power density
    Pa = f_power_density(Triangle_el, A, sigma, freq);
    % Scaling r factor
    Pa = Pa.*coo_el(:,1);
    
    % Thermal problem matrix assembling
    % Compute matrices for convective BC
%     [G_lhs, G_rhs] = f_BC_convective(conn_el(core_elems,:), coo_el(unique(conn_el(core_elems,:)),:), h, T_inf, edg, ibndedg, Triangle_el);
    [G_lhs, G_rhs] = f_BC_convective(conn_el(core_elems,:), coo_el, h, T_inf, edg, ibndedg);
    [STIFF_matrix, B_matrix] = f_FEM_dynamic_fast(conn_el(core_elems,:), coo_el, lambda, forcing_func_th);
    MASS_matrix = f_MASS_matrix(conn_el(core_elems,:), coo_el, q);
    % Matrices assembly
    R_matrix = MASS_matrix./timestep + theta.*(STIFF_matrix + G_lhs);
    L_matrix = MASS_matrix./timestep + (theta-1).*(STIFF_matrix + G_lhs);
    % Forcing function is constant in time
    f_vector = Pa + G_rhs;
    
    %% Solving electrical problem
    % Dirichlet BC electrical exploiting
    free_vars = setdiff(1:Nn_el, bound_nodes);
    A = zeros(Nn_el,1);
    % Solution vector
    A(free_vars) = LHS_global(free_vars,free_vars) \ RHS_global(free_vars);
    
    %% Solving thermal problem
    bk = L_matrix*oldsol_th+f_vector;
    sol_th = R_matrix\bk;
    % Billet centre temperature
    temp_centre(k) = f_pot_interp_T3([0 0], coo_el, conn_el, sol_th);
    
    %% Update solution
    oldsol_th = sol_th;
end

fprintf('Computed time theta method: %.3f sec\n', toc);
fprintf('\n');

%% POST-PROCESSING

%% Mesh plot
figure('Name', 'Mesh');
triplot(conn_el, R, Z, 'c', 'DisplayName', 'Triangle mesh');
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
plot(coo_el(bound_nodes,1), coo_el(bound_nodes,2), 'm', 'LineWidth', 3, 'DisplayName', 'Boundary');
xlabel('R-axis [m]');
ylabel('Z-axis [m]');
legend;
hold off;

%% Time-varying temperature T in the billet's centre
% % Import COMSOL output
% comsol_temp = readtable('./comsol_time_temp.txt', 'HeaderLines', 8);
% comsol_temp = table2array(comsol_temp);
% 
% figure('Name', 'Time-varying temperature T [K]');
% hold on;
% plot(times, temp_centre, 'b', 'DisplayName', 'Temperature T [K]');
% plot(comsol_temp(:,1), comsol_temp(:,2), '--r', 'DisplayName', 'Comsol data');
% hold off;
% title("Tym-varying temperature T [K] on billet's centre");
% xlabel('Time [s]');
% ylabel('Temperature T [K]');
% legend();
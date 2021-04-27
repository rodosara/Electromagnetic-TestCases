%% TestCase 6B - HEAT DYNAMICS TRANSIENT
clear all;
close all;
clc;

%% PRE-PROCESSING
% Elements along axes
nn = 100; % elements

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

% Bound nodes
bound_nodes = convhull(coo(:,1), coo(:,2))';

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
% Theta method with FEM
% Volume of all billet
radius_billet = billet_R_end-billet_R_start;
height_billet = billet_Z_end-billet_Z_start;
V = pi*radius_billet^2*height_billet;
% Area of all billet boundary
S = 2*pi*radius_billet*height_billet + 2*pi*radius_billet^2;
% Time costant [s]
tau = ro*Cp*V / (h*S);

% Theta method parameters
theta = 0.5;
Tmax = 8*tau; % [s]
timestep = 100; % [s]
itnum = round(Tmax/timestep);

% Diffusion coefficient in all domain
q = ones(Ne,1)*ro*Cp;
% Scaling axysymmetric factor
r = center(:,1);
q = q.*r;

% Thermal conductivity in all domain
lambda = ones(Ne,1)*lambda;
lambda = lambda.*r;

% Forcing function in all domain
forcing_func = ones(Ne,1)*g.*r;

% Compute the LHS and RHS global matices with fast method
tic;
[STIFF_matrix, B_matrix] = f_FEM_dynamic_fast(conn, coo, lambda, forcing_func);

% Mass matrix computation
MASS_matrix = f_MASS_matrix(conn, coo, q);

% Impose convective BC
% Compute boundary edge index and elements edge
[bndedg,edg] = f_quickedges(conn);
ibndedg = bndedg;
% Remove symmetry axes from boundary edges
ibndedg(coo(edg(bndedg,1),1) == 0) = [];
[G_lhs, G_rhs] = f_BC_convective(conn, coo, h, T_inf, edg, ibndedg);

% Matrices assembly
R_matrix = MASS_matrix./timestep + theta.*(STIFF_matrix + G_lhs);
L_matrix = MASS_matrix./timestep + (theta-1).*(STIFF_matrix + G_lhs);
% Forcing function is constant in time
f_vector = B_matrix + G_rhs;

% Theta method loop
times = zeros(itnum-1,1);
temp_centre = zeros(itnum-1,1);
% Initial temperature
oldsol = T_inf*ones(Nn,1);
for k = 1:itnum
    times(k) = (k-1)*timestep;
    bk = L_matrix*oldsol+f_vector;
    sol = R_matrix\bk;
    % Billet centre temperature
    temp_centre(k) = f_pot_interp_T3([0.05 0], coo, conn, sol);
    % Update solution
    oldsol = sol;
end

fprintf('Computed time of theta method with fast method matrices assembling: %.3f sec\n', toc);
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

% Time-varying temperature T in the billet's centre
% Import COMSOL output
comsol_temp = readtable('./comsol_time_temp.txt', 'HeaderLines', 8);
comsol_temp = table2array(comsol_temp);

figure('Name', 'Time-varying temperature T [K]');
hold on;
plot(times, temp_centre, 'b', 'DisplayName', 'Temperature T [K]');
plot(comsol_temp(:,1), comsol_temp(:,2), '--r', 'DisplayName', 'Comsol data');
hold off;
title("Time-varying temperature T [K] on billet's centre");
xlabel('Time [s]');
ylabel('Temperature T [K]');
legend();
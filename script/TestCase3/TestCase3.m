%% TestCase 3 - ELECTROSTATIC WITH DIFFERENT PERMITTIVITY
clear all;
close all;
clc;

%% PRE-PROCESSING
% Elements along axes
nn = 100; % elements

% External medium
% x dimension start
ext_x_start = -4e-2; % [m]
% x dimension end
ext_x_end = 4e-2; % [m]
% x vector
x = linspace(ext_x_start, ext_x_end, nn);
% y dimension start
ext_y_start = -2e-2; % [m]
% y dimension end
ext_y_end = 2e-2; % [m]
% y vector
y = linspace(ext_y_start, ext_y_end, nn);

% Internal medium
% x dimension start
int_x_start = -1e-2; % [m]
% x dimension end
int_x_end = 1e-2; % [m]
% y dimension start
int_y_start = -1e-2; % [m]
% y dimension end
int_y_end = 1e-2; % [m]

% Permittivity vacuum
epsilon_0 = 8.86e-12; % [F/m]
% Medium permittivity
% External medium
p_1 = epsilon_0; % [F/m]
% Factor multiply between the two permittivity medium
% pf= 4;
pf = 16; 
% Internal medium 
p_2 = pf*epsilon_0; % [F/m]

% Forcing function Laplacian PDE problem
forcing_func = @(x,y) 0.*x.*y;

% Voltage of top plate
Vtop = +1; % [V]
% Voltage of bottom plate
Vbottom = -1; % [V]

% Meshgrid combination x and y
[X,Y] = meshgrid(x,y);

% Coordinate matrix
coo = [X(:), Y(:)];
% Connectivity matrix
conn = delaunay(coo(:,1), coo(:,2));

% Number of nodes
Nn = size(coo,1);
% Number of elements
Ne = size(conn,1);

% Bound nodes
bound_nodes = convhull(coo(:,1), coo(:,2))';

% Dirichlet boundary conditions parameters
% Find nodes of top and bottom capacitor's plates
bottom_plate = bound_nodes(1:nn);
top_plate = bound_nodes(2*nn-1:end-nn+1);

% Found internal dielectric dishomogenus medium
int_nodes = find(coo(:,1) >= int_x_start & coo(:,1) <= int_x_end & ...
    coo(:,2) >= int_y_start & coo(:,2) <= int_y_end);

Triangle = triangulation(conn, coo);
center = incenter(Triangle);
% Find elements of internal medium using the baricenter points
int_elems = find(center(:,1) >= int_x_start & center(:,1) <= int_x_end & ...
    center(:,2) >= int_y_start & center(:,2) <= int_y_end);

% Define new p vector considering two medium
p = ones(Ne,1) * p_1;
% Insert second medium
p(int_elems) = p_2;

%% PROCESSING

% Compute the LHS and RHS global matices with fast method
tic;
[LHS_global, RHS_global] = f_FEM_global_fast(conn, coo, p, forcing_func);
fprintf('Computed time of global matrices with fast method: %.3f sec\n', toc);
fprintf('\n');

% Dirichlet condition exploiting
free_vars = setdiff(1:Nn, [top_plate, bottom_plate]);
u = zeros(Nn,1);
u(top_plate) = Vtop;
u(bottom_plate) = Vbottom;
RHS_global = RHS_global - LHS_global*u;
% Solution vector
u(free_vars) = LHS_global(free_vars,free_vars) \ RHS_global(free_vars);

%% POST-PROCESSING
% Electric field magnitude
[Ex, Ey] = f_Efield_interp_T3(coo, coo, conn, u);
% Compute the magnitude
E = sqrt(Ex.^2 +Ey.^2);

% Capacitance
% Voltage applied at plates [V]
V0 = Vtop-Vbottom;
[C, W] = f_capacitance(conn, coo, u, p, V0);
% Print output
format shortEng;
fprintf('Capacitance: %d [F]\n', C);
fprintf('Electrostatic energy: %d [J]\n', W);
fprintf('\n');


% Sparsity plot
figure('Name', 'Sparsity LHS Matrix');
spy(LHS_global);
title('Sparsity LHS Matrix');

% Mesh plot
figure('Name', 'Mesh');
triplot(conn,X,Y);
title('Mesh figure');
hold on;
% Capacitor's plates - Dirichlet bounds
plot(coo(top_plate,1), coo(top_plate,2), 'k', 'LineWidth', 4);
plot(coo(bottom_plate,1), coo(bottom_plate,2), 'k', 'LineWidth', 4);
% Second medium bounds
square = [[int_x_start int_x_end int_x_end int_x_start int_x_start]', ...
   [int_y_start int_y_start int_y_end int_y_end int_y_start]'];
plot(square(:,1), square(:,2), 'r', 'LineWidth', 2);
xlabel('x [m]');
ylabel('y [m]');
axis equal;
legend('Mesh', 'Top plate', 'Bottom plate', 'Internal medium');
hold off;

% Potential contour plot
figure('Name', 'Potential');
hold on;
% Contour potential plot
contour(X, Y, reshape(u,[nn,nn]));
title('Potential');
% Capacitor's plates
plot(coo(top_plate,1), coo(top_plate,2), 'k', 'LineWidth', 3);
plot(coo(bottom_plate,1), coo(bottom_plate,2), 'k', 'LineWidth', 3);
% Second medium bounds
square = [[int_x_start int_x_end int_x_end int_x_start int_x_start]', ...
   [int_y_start int_y_start int_y_end int_y_end int_y_start]'];
plot(square(:,1), square(:,2), 'r', 'LineWidth', 1.5);
% Label permittivity
text(-0.002, 0.005,['\epsilon_{r}=' num2str(pf) '\epsilon_{0}'], 'FontSize', 20);
xlabel('x [m]');
ylabel('y [m]');
legend('Potential [V]', 'Top plate', 'Bottom plate', 'Internal medium');
axis equal;
colorbar;

% Electric field magnitude plot
figure('Name', 'Electric field magnitude');
hold on;
% Electric field magnitude plot
imagesc(x, y, reshape(E, nn, nn), 'CDataMapping', 'scaled', 'Interpolation', 'bilinear');
title('Electric field magnitude');
% Contour E field plot in white
contour(X, Y, reshape(E,[nn,nn]), 'w');
% Capacitor's plates
plot(coo(top_plate,1), coo(top_plate,2), 'k', 'LineWidth', 3);
plot(coo(bottom_plate,1), coo(bottom_plate,2), 'k', 'LineWidth', 3);
% Second medium bounds
plot(square(:,1), square(:,2), 'r', 'LineWidth', 1.5);
% Label permittivity
text(-0.002, 0.005, ['\epsilon_{r}=' num2str(pf) '\epsilon_{0}'], 'FontSize', 20, 'Color', 'w');
xlabel('x [m]');
ylabel('y [m]');
axis equal;
colorbar;
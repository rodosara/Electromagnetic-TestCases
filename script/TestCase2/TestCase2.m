%% TestCase 2 - 2-D ELECTROSTATIC
clear all;
close all;
clc;

%% PRE-PROCESSING
% Elements along axes
nx = 100;
% x dimension start
x_start = 0;
% x dimension end
x_end = 1;
% x vector
x = linspace(x_start, x_end, nx);
% Square domain
y = x;
% Function 1
u_func = @(x,y) x.*(1-x).*y.*(1-y);
forcing_func = @(x,y) 2.*(x+y-x.^2 -y.^2);
% Function 2
u_func = @(x,y) sin(2*pi.*x).*sin(2*pi.*y);
forcing_func = @(x,y) 8*pi^2.*sin(2*pi.*x).*sin(2*pi.*y);
% Meshgrid combination x and y
[X,Y] = meshgrid(x, y);
% Permittivity (homogeneus medium ex=ey)
p = 1;

% Coordinate matrix
coo = [X(:), Y(:)];
% Connectivity matrix
conn = delaunay(coo(:,1), coo(:,2));

% Number of nodes
Nn = size(coo,1);
% Number of elements
Ne = size(conn,1);

% Dirichlet boundary conditions parameters
% Dirichlte nodes
bound_nodes = convhull(coo(:,1), coo(:,2))';
bound_nodes = bound_nodes(1:end-1);

%% PROCESSING
% Compute the RHS and LHS global matices with for loop method
tic;
[LHS_global, RHS_global] = f_FEM_global_for(conn, coo, p, forcing_func);
fprintf('TIME METHODS COMPARISON\n');
fprintf('Computed time of global matrices with for loop method: %.3f sec\n', toc);

% Compute the LHS and RHS global matices with fast method
tic;
[LHS_global, RHS_global] = f_FEM_global_fast(conn, coo, p, forcing_func);
fprintf('Computed time of global matrices with fast method: %.3f sec\n', toc);
fprintf('\n\n');

% Dirichlet condition
free_vars = setdiff(1:Nn, bound_nodes);
u = zeros(Nn,1);
u(bound_nodes) = u_func(coo(bound_nodes,1), coo(bound_nodes,2));
RHS_global = RHS_global - LHS_global*u;
% Solution vector
u(free_vars) = LHS_global(free_vars,free_vars) \ RHS_global(free_vars);

%% POST PROCESSING

% Electric field
% Create object triangulation
Triangle = triangulation(conn, coo);
% Incenter (bisectors center) of triangles
center = incenter(Triangle);

[Ex, Ey] = f_Efield_interp_T3(center, coo, conn, u);

% Error
u_analytical = u_func(coo(:,1), coo(:,2));
err_perct = max(abs(u-u_analytical))/max(abs(u_analytical))*100;
fprintf('COMPUTED FEM SOLUTION U\n- Node number: %d\n- Error potential U: %e\n', Nn, err_perct);

% Sparsity A_global matrix plot
figure('Name', 'Sparsity A_global Matrix');
spy(LHS_global);
title('Sparsity A global Matrix');
figure('Name', 'Mesh figure');
triplot(conn,X,Y);
title('Mesh figure');
hold on;
% Boundary nodes
scatter(coo(bound_nodes',1), coo(bound_nodes',2), '*', 'LineWidth', 1.5);
% Center nodes
scatter(center(:,1), center(:,2), '+', 'LineWidth', 1.5);
xlabel('x dimension');
ylabel('y dimension');
axis equal;
axis([-0.5 +1.5 -0.5 +1.5]);
legend('mesh','boundary nodes');

% Computed and analytical solution comparison plot
figure('Name', 'Computed and analytical solution comparison');
% Analytical solution
subplot(1,3,1);
trisurf(conn, X, Y, reshape(u_analytical,[nx,nx]));
title('Analytical solution');
xlabel('x');
ylabel('y');
zlabel('V/m');
% Computed solution
subplot(1,3,2);
trisurf(conn, X, Y, reshape(u,[nx,nx]));
title('Computed solution');
xlabel('x');
ylabel('y');
zlabel('V/m');
% Error
subplot(1,3,3);
error = u_analytical - u;
trisurf(conn, X, Y, reshape(error,[nx,nx]));
title('Error');
xlabel('x');
ylabel('y');
zlabel('V/m');

% Electric field plot
figure('Name', 'Electric field');
% Surface plot
subplot(1,2,1);
trisurf(conn, X, Y, reshape(u,[nx,nx]));
title('Potential');
xlabel('x');
ylabel('y');
zlabel('V/m');
% Isolines plot
subplot(1,2,2);
contour(X, Y, reshape(u,[nx,nx]));
title('Electric field');
hold on;
quiver(center(:,1), center(:,2), Ex, Ey);
xlabel('x');
ylabel('y');
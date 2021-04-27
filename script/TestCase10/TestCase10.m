%% TestCase 10 - BEM
clear all;
close all;
clc;

%% PRE-PROCESSING
% Elements along axes
nx = 50;
% x dimension start
x_start = 0;
% x dimension end
x_end = 1;
% x vector
x = linspace(x_start, x_end, nx);
% Square domain
y = x;
% Function 1
% u_func = @(x,y) x.*(1-x).*y.*(1-y);
% forcing_func = @(x,y) 2.*(x+y-x.^2 -y.^2);
% Function 2
u_func = @(x,y) sin(2*pi.*x).*sin(2*pi.*y);
forcing_func = @(x,y) 8*pi^2.*sin(2*pi.*x).*sin(2*pi.*y);
% Meshgrid combination x and y
[X,Y] = meshgrid(x,y);
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

% Boundary nodes
bound_nodes = convhull(coo(:,1), coo(:,2));
% Boundary edge matrix (ExE)
bnd_edg = [bound_nodes(1:end-1), bound_nodes(2:end)];
% Potential values for each node
dir_v = zeros(size(bound_nodes));

%% PROCESSING

% FEM section
tic;
[LHSfem, RHSfem] = f_FEM_global_fast(conn, coo, p, forcing_func);
% Boundary Dirichlet condition
free_vars = setdiff(1:Nn, bound_nodes);
Ufem = zeros(Nn,1);
Ufem(bound_nodes) = u_func(coo(bound_nodes,1), coo(bound_nodes,2));
RHSfem = RHSfem - LHSfem*Ufem;
% Solution vector
Ufem(free_vars) = LHSfem(free_vars,free_vars)\RHSfem(free_vars);
time_fem = toc;

% BEM section

% Compute matrices for processing using Gaussian quadrature
% in order to have the same points of Collocation method with FEM mesh
tic;
% G matrix computation (ExE)
Gbem = f_assemble_Gmat_BEM_2d(coo, bnd_edg);
% Forcing function matrix computation (Ex1)
Fbem = f_assemble_Svec_BEM_2d(coo, conn, bnd_edg, forcing_func);
% Rewrite in form Hu-GDnu=F
Amat = -Gbem;
sol_BEM = Amat\Fbem;
bnd_Neumann = sol_BEM;
% Find the solution of internal points from the boundary points
[Fvec_post, Gmat_post] = f_postpro_BEM_2d_zeropot(coo, conn, bnd_edg, forcing_func);
Ubem = Fvec_post + Gmat_post * bnd_Neumann;
time_bem = toc;

%% POST-PROCESSING
% Result analysis

% Analytical solution
u_analytical = u_func(coo(:,1), coo(:,2));

% FEM error
err_perct = max(abs(Ufem-u_analytical))/max(abs(u_analytical))*100;
fprintf('COMPUTED FEM SOLUTION U\n');
fprintf('- Computed time: %.3f sec\n', time_fem);
fprintf('- Node number: %d\n', Nn);
fprintf('- Error potential U: %e\n', err_perct);
fprintf('\n');

% BEM error
err_perct = max(abs(Ubem-u_analytical))/max(abs(u_analytical))*100;
fprintf('COMPUTED BEM SOLUTION U\n');
fprintf('- Computed time: %.3f sec\n', time_bem);
fprintf('- Node number: %d\n', Nn);
fprintf('- Error potential U: %e\n', err_perct);

% Sparsity matrices comparison plot
figure('Name', 'Sparsity matrices comparison');
subplot(1,2,1);
spy(LHSfem);
title('Sparsity K FEM matrix');
subplot(1,2,2);
spy(Gbem);
title('Non-sparsity G BEM matrix');

% Computed and analytical solution comparison plot
figure('Name', 'Computed and analytical solution comparison');
% Analytical solution
subplot(1,3,1);
trisurf(conn, X, Y, reshape(u_analytical,[nx,nx]));
title('Analytical solution');
xlabel('x');
ylabel('y');
zlabel('V/m');
% Computed FEM solution
subplot(1,3,2);
trisurf(conn, X, Y, reshape(Ufem,[nx,nx]));
title('Computed FEM solution');
xlabel('x');
ylabel('y');
zlabel('V/m');
% Computed BEM solution
subplot(1,3,3);
trisurf(conn, X, Y, reshape(Ubem,[nx,nx]));
title('Computed BEM solution');
xlabel('x');
ylabel('y');
zlabel('V/m');

% Error between BEM and FEM plot
figure('Name', 'Error between BEM and FEM');
u_comparison = Ubem - Ufem;
trisurf(conn, X, Y, reshape(u_comparison,[nx,nx]));
title('Error between BEM and FEM');
xlabel('x');
ylabel('y');
zlabel('V/m');
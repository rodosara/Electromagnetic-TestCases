%% TestCase 1 - ELECTROSTATIC
clear all;
close all;
clc;

%% PRE-PROCESSING
% Total size [cm]
L = 1;
% Number of elements
Ne = 4;
% Vacuum permeability [F/m]
e_0 = 8.86e-12;
% Density charge
ro = 0;
% Scattering parameter
q = 0;

% Size of elements
le = L/Ne;
% Number of nodes
Nn = Ne+1;

% Dirichlet boundary conditions parameters
% Dirichlte nodes
dir_node = [1 Nn];
% Potential values for each node
dir_v = [1 0];

% Mixed boundary conditions !! ATTENTION --> with sign of alpha and beta
% give the flux direction (in computation all node are summed)

% Mixed nodes
M_node = [];
% Alpha values for each node
alpha = [];
% Beta values for each node
beta = [];

%% PROCESSING
% Connectivity matrix
conn = [(1:Nn-1)' (2:Nn)'];

% Global assembly matrix
A = sparse(Nn,Nn); B = sparse(Nn,1);
for el = 1:Ne
    % Omogeneus medium
    epsilon = e_0;
    % Inomogeneus medium
    epsilon = e_0*(1+el*L^2/Ne);
    A_e = (epsilon/L) .* sparse([1,1,2,2],[1,2,1,2],[1,-1,-1,1]) + (q*L/6) * sparse([1,1,2,2],[1,2,1,2],[2,+1,+1,2]);
    B_e = (0.5*ro*L) .* ones(Nn,1);
    for k = 1:size(conn,2)
        kg = conn(el,k);
        for h = 1:size(conn,2)
            hg = conn(el,h);
            A(kg,hg) = A(kg,hg) + A_e(k,h);
        end
        B(kg) = B(kg) + B_e(k);
    end
end

% Dirichlet condition
for ii = 1:length(dir_node)
    node = dir_node(ii);
    B = B - A(:,node) * dir_v(ii);
end
A(:, dir_node) = [];
A(dir_node, :) = [];
B(dir_node) = [];

% Mixed condition
for ii = 1:length(M_node)
    node = M_node(ii);
    A(node,node) = A(node,node) + alpha(ii);
    B(node,1) = B(node,1) + beta(ii);
end

% Solution vector
u = A\B;

% Add Dirichlet condition to the solution
u = cat(1, dir_v(1), u);
u = cat(1, u, dir_v(2));

% Compute electric field
E = -diff(u) / (le/100);

%% POST-PROCESSING
% Analytical solution
ll = linspace(0, L, Nn);
u_analytical = dir_v(1) + ((dir_v(2)-dir_v(1))/log(2)).*log(1+ll);
E_analytical = ((dir_v(1)-dir_v(2))/log(2)) .* (1./(ll+L))*100;

% Computed and analytical solution comparison plot
figure('Name', 'Computed and analytical solution comparison');
subplot(2,1,1);
title('Electric potential [V]');
hold on;
plot(0:le:L, u, 'o', 'DisplayName', 'Computed solution');
plot(0:le:L, u_analytical, 'r', 'DisplayName', 'Analytical solution');
hold off;
xlabel('x-axis [cm]');
ylabel('Volt [V]');
legend;
grid on;

subplot(2,1,2);
title('Electric field [V/m]');
hold on;
stairs(0:le:L, [E; E(end)], 'DisplayName', 'Computed solution');
plot(0:le:L, E_analytical, 'r', 'DisplayName', 'Analytical solution');
hold off;
xlabel('x-axis [cm]');
ylabel('E [Volt/m]');
legend;
grid on;
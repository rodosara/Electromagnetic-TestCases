function [nod, edg, ph_edg] = discretize_busbar(L, d, Ne)
% Phase centre with coordinate in 3 axis X, Y and Z
Yph = [[-d 0 d]; zeros(1,3); zeros(1,3)];

% Start point of busbar
A = [Yph(1,:); Yph(2,:); (L/2)*ones(1,3)];
% End point of busbar
B = [Yph(1,:); Yph(2,:); -(L/2)*ones(1,3)];

nod = zeros(3,0);
for ii = 1:3
    % Discretized one bus_bar
    wire = linspace(L/2, -L/2, Ne+1);
    
    % Node matrix [3xNn] with coordinates of each nodes
    nod_wire = [Yph(1,ii)*ones(1, length(wire)); zeros(1, length(wire)); wire];
    nod = [nod, nod_wire];
end

% Edges matrix [Nx3] numbering each nodes in column
edg = [(1:Ne)', (2:Ne+1)'; Ne+1+[(1:Ne)', (2:Ne+1)']; 2*(Ne+1)+[(1:Ne)', (2:Ne+1)']];

% Matrix to label the phase 1 for U, 2 for V and 3 for W
ph_edg = [ones(1,Ne), 2*ones(1,Ne), 3*ones(1,Ne)];
end
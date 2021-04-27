function sigma = f_sigma_update(conn, coo, T, sigma_core, alpha, T0, medium_elems)

% INPUT
% conn: connection matrix (Nex3) with global numbering node connection
% coo: coordinates matrix (Nnx2) with coordinates of all nodes
% T: vector (Nn_thx1) of temperature for each node from thermal problem
% sigma_core: value of sigma in the medmium considered
% medium_elems = index of elements in the medium where compute sigma
% T0: reference temperature
% alpha: temperature coefficient

% OUTPUT
% sigma = vector (Ne_elx1) sigma values for each element for  electrical meshes

% Compute temperature in the centre of element from node values
T_medium = @(T_el) (1/3) .* sum(T_el(conn(medium_elems,:)),2);

% Number of nodes electrical mesh
Nn_el = size(coo,1);
% Number of elements electrical mesh
Ne_el = size(conn,1);

% Assign electrical mesh nodes the value of temperature obtained from
% thermal problem
T_el = zeros(Nn_el,1);
T_el(unique(conn(medium_elems,:))) = T;

% Compute sigma medium in each elements from electrical mesh nodes
T_elements = T_medium(T_el);

sigma = zeros(Ne_el,1);
% Update sigma values
sigma(medium_elems) = (1/sigma_core) + alpha.*(T_elements-T0);

end
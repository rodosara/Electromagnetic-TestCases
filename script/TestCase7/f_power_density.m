function Pa = f_power_density(Triangle, A, sigma, freq)

% INPUT
% Triangle: Matlab object of mesh with connectivity and coordinate matrix
% A: vector (Nn_thx1) of density current from electrical problem
% sigma: array (Nn_thx1) of sigma values for each electrical mesh elements
% freq: frequency of electrical problem

% OUTPUT
% Pa = vector (Ne_elx1) of Joule losses power density for each nodes

% Number of electrical mesh nodes
Nn_el = length(Triangle.Points);

sigma_nodes = zeros(Nn_el,1);
for kk = 1:Nn_el
    % Elements connected to a node (vertex)
    elems_attch = vertexAttachments(Triangle, kk);
    elems_attch = cell2mat(elems_attch);

    % Value of sigma in the node from medium of all elements
    num_elems = length(elems_attch);
    sigma_nodes(kk) = (1/num_elems) * sum(sigma(elems_attch));
end

% Compute Joule losses power density array
Pa = 0.5*(2*pi*freq)^2.*sigma_nodes.*A.^2;
end
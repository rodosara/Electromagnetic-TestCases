function joule_losses = f_joule_losses(coo, conn, sigma, freq, u, medium_elems)

%% Function to compute Joule losses in particoular medium

% INPUT
% conn: connection matrix (Nex3) with global numbering node connection
% coo: coordinates matrix (Nnx2) with coordinates of all nodes
% sigma: conducibility considering an homogeneus medium
% freq: frequency [Hz] (scalar)
% u: magnetic vector potential A'_theta (Nnx1)
% medium_elemes: array of elements where compute Joule losses (Mx1)

% OUTPUT
% joule_losses = value of Joule losses compute in the medium (scalar)

% In-line function areatri
areatri = @(nodes_elem) .5*norm(cross([nodes_elem(3,:)-nodes_elem(2,:) 0], ...
    [nodes_elem(1,:)-nodes_elem(3,:) 0]));

% Triangle object
Triangle = triangulation(conn, coo);
% Compute center coordinate
center = incenter(Triangle);
% Axisymmetric scaling factor
r = center(:,1);

% W
w = 2*pi*freq;

% Consider only the medium where compute Joule losses
center = center(medium_elems, :);
r = r(medium_elems, :);
sigma = sigma(medium_elems, :);
conn = conn(medium_elems, :);

% Mesh is composed by same area triangles
area = (areatri(coo(conn(1,:),:)));

% Compute the magnetic vector potential A' in the center
A_medium = (1/3) .* sum(u(conn),2);
% Calc the Joule losses with analytical formula
joule_losses = pi*w^2 * sum((area.*sigma.*abs(A_medium).^2) ./ r);

% % Magnetic vector potential medium in the elements (sum of nodes value)
% summation = 0;
% for elem = 1:length(conn)
%     nodes_elem = conn(elem,:);
%     A_medium = (1/3) * sum(u(nodes_elem));
%     summation = summation + (areatri(coo(nodes_elem,:))*sigma(elem)*A_medium^2) / r(elem);
% end
%
% joule_losses = pi*w^2 * summation;

end
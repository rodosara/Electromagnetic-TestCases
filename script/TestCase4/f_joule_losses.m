function [joule_losses] = f_joule_losses(coo, conn, sigma, freq, u)

% In-line function areatri
areatri = @(nodes_elem) .5*norm(cross([nodes_elem(3,:)-nodes_elem(2,:) 0], [nodes_elem(1,:)-nodes_elem(3,:) 0]));

% Triangle object
Triangle = triangulation(conn, coo);
% Compute center coordinate
center = incenter(Triangle);
% Axisymmetric scaling factor
r = center(:,1);

% W
w = 2*pi*freq;

% Magnetic vector potential medium in the elements (sum of nodes value)
for elem = 1:length(conn)
   nodes_elem = conn(elem,:);
   A_medium = (1/3) * sum(u(tricell));
   summation = (area*sigma(elem)*A_medium) / r(elem);
end

   joule_losses = pi*w^2 * summation;
   
%       A_medium = (1/3) .* sum(u(tricell),2);
%       joule_losses = sum(pi*w^2*(area.*sigma.*A_medium) ./ r);

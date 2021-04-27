function [C, W] = capacitance(conn, coo, u, p, V0)

%% Function to compute capacitance C and electrostatic energy W

% INPUT
% conn: connection matrix (Nex3) with global numbering node connection
% coo: coordinates matrix (Nnx2) with coordinates of all nodes
% u: vector potential (Nnx1)
% p: permettivity considering an homogeneus medium
% V0: voltage applied at plates [V]

% OUTPUT
% C = capacitance [F]
% W = electrostatic energy [J]

% Number of nodes
Nn = size(coo,1);
% Number of elements
Ne = size(conn,1);

Triangle = triangulation(conn, coo);
% Normal face vector
ez = repmat(faceNormal(Triangle), [1,1,3]);
% Incenter (bisectors center) of triangles
center = incenter(Triangle);

% Edge from vertex 3 to 2
l1 = cat(2, coo(conn(:,3),:) - coo(conn(:,2),:), zeros(Ne,1));
% Edge from vertex 1 to 3
l2 = cat(2, coo(conn(:,1),:) - coo(conn(:,3),:), zeros(Ne,1));
% Edge from vertex 2 to 1
l3 = cat(2, coo(conn(:,2),:) - coo(conn(:,1),:), zeros(Ne,1));
% All edges in a 3 dim matrix
edges = cat(3, l1, l2, l3);
% Area
area = .5 * norm(cross(l1(1,:),l2(1,:)));
% Gradient (W) matrix (Nex3x3)
gradient = cross(ez, edges) ./ (2*area);
% Gradient of U vector potential
grad_u = gradient(:,:,1).*u(conn(:,1)) + gradient(:,:,2).*u(conn(:,2)) ...
    + gradient(:,:,3).*u(conn(:,3));
% Capacitace formula [F]
C = p.*(grad_u(:,1).^2 + grad_u(:,2).^2 + grad_u(:,3).^2)*area;
C = 1/V0^2*sum(C);

% Electrostatic energy [J]
W = 0.5*C*V0^2;

end
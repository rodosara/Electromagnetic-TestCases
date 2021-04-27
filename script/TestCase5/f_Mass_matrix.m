function M_global = Mass_matrix(conn, coo, q)

%% Function to create both RHS and LHS global matrices
% Use the for loop method (more computational cost) and sparse matrices

% INPUT
% conn: connection matrix (Nex3) with global numbering node connection
% coo: coordinates matrix (Nnx2) with coordinates of all nodes
% p: permettivity considering an homogeneus medium
% forcing_func: forcing function provided

% OUTPUT
% LHS_global = sparse matrix (NnxNn) of right hand side (stifness matrix)
% RHS_global = sparse matrix (Nnx1) of left hand side

% In-line function areatri
areatri = @(nodes_elem) .5*norm(cross([nodes_elem(3,:)-nodes_elem(2,:) 0], ...
    [nodes_elem(1,:)-nodes_elem(3,:) 0]));
% Mesh is composed by same area triangles
area = (areatri(coo(conn(1,:),:)));

% Number of nodes
Nn = size(coo,1);
% Number of elements
Ne = size(conn,1);

% Matrix of incident for compute M_local
M_incident = spones(3) + spdiags([1; 1; 1], 0, 3, 3);

M_global = sparse(Nn,Nn);
for elem = 1:Ne
    % Local M matrix computation
    M_local = (area*q(elem)/12) .* M_incident;
    
    % Global assembly matrix
    for kk = 1:3
        kg = conn(elem,kk);
        for hh = 1:3
            hg = conn(elem,hh);
            M_global(kg,hg) = M_global(kg,hg) + M_local(kk,hh);
        end
    end
end

end
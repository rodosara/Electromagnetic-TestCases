function [LHS_global, RHS_global] = FEM_global_fast(conn, coo, p, forcing_func)

%% Function to create both RHS and LHS global matrices
% Use the fast method (less computational cost) and sparse matrices

% INPUT
% conn: connection matrix (Nex3) with global numbering node connection
% coo: coordinates matrix (Nnx2) with coordinates of all nodes
% p: permettivity considering an homogeneus medium
% forcing_func: forcing function provided

% OUTPUT
% LHS_global = sparse matrix (NnxNn) of right hand side (stifness matrix)
% RHS_global = sparse matrix (Nnx1) of left hand side

% Number of nodes
Nn = size(coo,1);
% Number of elements
Ne = size(conn,1);

% Compute triangle features (in this particoular case are equilateral,
% so it could be put output for loop)
% Create object triangulation
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

% Matrix of all values
value = zeros(9*Ne,1);
% Matrix of row indexes
ir = zeros(9*Ne,1);
% Matrix of column indexes
ic = zeros(9*Ne,1);
count = 0;
for mm = 1:3
    for nn = 1:3
        first = count*Ne + 1;
        last = first+Ne - 1;
        ind = first:last;
        ir(ind) = conn(:,mm);
        ic(ind) = conn(:,nn);
        value(ind) = p.*area .* dot(gradient(:,:,mm), gradient(:,:,nn), 2);
        count = count+1;
    end
end

% Fast assembly LHS matrix A (NnxNn)
LHS_global = sparse(ir, ic, value, Nn, Nn);

% Fast assembly RHS matrix B (Nnx1)
bloc = 1/3*area*forcing_func;%(center(:,1), center(:,2));
for ii = 1:Nn
    [row, col] = find(conn==ii);
    RHS_global(ii,1) = sum(bloc(row));
end

end
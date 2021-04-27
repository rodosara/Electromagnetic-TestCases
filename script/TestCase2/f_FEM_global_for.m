function [LHS_global, RHS_global] = FEM_global_for(conn, coo, p, forcing_func)

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

% Number of nodes
Nn = size(coo,1);
% Number of elements
Ne = size(conn,1);

% Compute triangle features (in this particoular case are equilateral,
% so it could be put output for loop)
% Create object triangulation
Triangle = triangulation(conn, coo);
% Normal face vector
ez = faceNormal(Triangle);
% Incenter (bisectors center) of triangles
center = incenter(Triangle);

LHS_global = sparse(Nn,Nn);
RHS_global = sparse(Nn,1);
for el = 1:Ne
    % Local A matrix computation
    
    % Edge from vertex 2 to 3
    l1 = cat(2, coo(conn(el,3),:) - coo(conn(el,2),:), 0);
    % Edge from vertex 3 to 1
    l2 = cat(2, coo(conn(el,1),:) - coo(conn(el,3),:), 0);
    % Edge from vertex 1 to 2
    l3 = cat(2, coo(conn(el,2),:) - coo(conn(el,1),:), 0);
    % All edges in a vector
    edges = cat(1, l1, l2, l3);
    % Area
    area = .5 * norm(cross(l1, l2));
    
    % Gradient computation
    gradient = zeros(3);
    for mm = 1:3
        gradient(mm,:) = cross(ez(el,:), edges(mm,:)) / (2*area);
    end
    
    % LHS local
    a_local = zeros(3);
    for mm = 1:3
        for nn = 1:3
            a_local(mm,nn) = area*p*dot(gradient(mm,:), gradient(nn,:));
        end
    end
    
    % RHS local
    density = forcing_func(center(el,1), center(el,2));
    b_local = 1/3*area*density * ones(Nn,1);
    
    % Global assembly matrix
    for kk = 1:3
        kg = conn(el,kk);
        for hh = 1:3
            hg = conn(el,hh);
            LHS_global(kg,hg) = LHS_global(kg,hg) + a_local(kk,hh);
        end
        RHS_global(kg) = RHS_global(kg) + b_local(kk);
    end
end

end
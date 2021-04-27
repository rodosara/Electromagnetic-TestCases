function  [G_lhs_global, G_rhs_global] = f_BC_convective(conn, coo, h, T_inf, edg, bound_edges)

%% Function to G_rhs and G_lhs matrices for conductive BC in axisymmetric problems

% INPUT
% conn: connection matrix (Nex3) with global numbering node connection
% coo: coordinates matrix (Nnx2) with coordinates of all nodes
% h: heat trasfer rate
% T_inf: unperturbed bulk temperature
% edges: matrix (Px2) with the start-nodes and end-nodes of element's edge
% bound_edges: vector (Ex1) with the index of boundary edges

% OUTPUT
% G_lhs = sparse matrix (NnxNn) of right hand side
% G_lhs = sparse matrix (Nnx1) of left hand side

% Lenght of element in-line function (perimenter equilateral triangle)
length_elem = @(edgs) norm(coo(edgs(1),:) - coo(edgs(2),:));

% Number of nodes
Nn = size(coo,1);
% Number of elements
Ne = size(conn,1);

% Create object triangulation
Triangle = triangulation(conn, coo);

% Identify the element of boundary edges
bound_elems = edgeAttachments(Triangle, edg(bound_edges,1), edg(bound_edges,2));
bound_elems = cell2mat(bound_elems);

G_lhs_global = sparse(Nn,Nn);
G_rhs_global = sparse(Nn,1);
for elem = 1:length(bound_edges)
    Le = length_elem(edg(bound_edges(elem),:));
    r1 = coo(edg(bound_edges(elem),1),1);
    r2 = coo(edg(bound_edges(elem),2),1);
    % Local G_lhs matrix computation
    G_lhs_local = (Le*h/12) .* [3*r1+r2 r1+r2; r1+r2 r1+3*r2];
    % Local G_rhs matrix computation
    G_rhs_local = (Le*h*T_inf/6) .* [2*r1+r2; r1+2*r2];
    
    % Global assembly matrixes
    for kk = 1:2
        kg = conn(bound_elems(elem),kk);
        for hh = 1:2
            hg = conn(bound_elems(elem),hh);
            G_lhs_global(kg,hg) = G_lhs_global(kg,hg) + G_lhs_local(kk,hh);
        end
        G_rhs_global(kg) = G_rhs_global(kg) + G_rhs_local(kk);

    end
end

end
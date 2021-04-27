function [bndedg,edg]=f_quickedges(tri)

% =========================================================
% author: Federico Moro (federico.moro@unipd.it)
% affiliation: Dip. Ing. Industriale, Universit� di Padova
% date: 31-03-20
% notes: compute the edge matrix and the indices of the edges
%        lying on the boundary of a triangle mesh.
% input: conn. matrix tri (Tx3)
% output: bnd edge indices bndedg (Ix1), edge matrix edg (Ex2) 
% =========================================================

% unsorted boundary (OLD vs)
% edg=[tri(:,[2 3]); tri(:,[3 1]); tri(:,[1 2])];
% [edg,m,n] = unique(sort(edg,2),'rows');
% counts = accumarray(n(:),1);
% bnd = edg(counts == 1,:);
% bndedg = find(ismember(edg,bnd,'rows'));

% sorted boundary!! (exterior normals)
edg=[tri(:,[2 3]); tri(:,[3 1]); tri(:,[1 2])];
[~,m,n] = unique(sort(edg,2),'rows');
edg=edg(m,:); 
counts = accumarray(n(:),1);
bnd = edg(counts == 1,:);
bndedg = find(ismember(edg,bnd,'rows'));

end
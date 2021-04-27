function Fvec = assemble_Svec_BEM_2d(nod,tri,bndedg,fun)

%=====================================================================
% file: "assemble_Svec_BEM_2d.m"
% author: F. Moro
% affiliation: Universit� degli Studi di Padova, Italy
% date: 14-05-20
% note: source array assembly
% input: nod (coordinate matrix: N x 3), tri (element matrix: T x 3), 
% bndedg (boundary edge matrix: E x 2), fun (source density function)
% output: Fvec (source vector: E x 1)
% ref. C.W.Trowbridge, "An Introduction to CAE Analysis", Vector Fields Ltd.
%=====================================================================

Fvec = zeros(size(bndedg,1),1);

% source term
for k=1:size(bndedg,1)
    edge = bndedg(k,:); % any row is a pair of nodal indices
    centre = 0.5*(nod(edge(1),:)+nod(edge(2),:));
    % collocation point is at the element centre
    for m=1:size(tri,1)
        cell_m = tri(m,:);
        nodtri=nod(cell_m,:);
        %         area=areatri(nodtri);
        area = .5 * norm(cross([nodtri(3,:)-nodtri(2,:) 0], [nodtri(1,:)-nodtri(3,:) 0]));
        % nine-point quadrature formula
        gaussnods=1/18*[1 4 13; 4 1 13; 1 13 4; 13 1 4; 4 7 7;...
            7 4 7; 7 7 4; 4 13 1; 13 4 1]*nodtri;
        distvecs = repmat(centre,9,1)-gaussnods;
        dists = sqrt(sum(distvecs.^2,2));
        weight=area/9;
        % integration on the whole domain
        Fvec(k) = Fvec(k)+weight/(2*pi)*sum(log(1./dists)...
            .*fun(gaussnods(:,1),gaussnods(:,2)));
    end
end

end


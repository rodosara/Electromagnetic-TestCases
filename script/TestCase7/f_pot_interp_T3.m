function [intu]=pot_interp_T3(points,nod,tri,u)

%=====================================================================
% file: "pot_interp_T3.m"
% author: Federico Moro
% affiliation: Universit� degli Studi di Padova
% date: 20-03-20
% note: interpolate potentials on a triangle mesh (T3 elemns)
% input: query coord. matrix points (Px2), coord. matrix nod (Nx2), 
%        conn. matrix tri (Tx3), potentials u (Nx1)
% output: potentials intu (Px1) 
%=====================================================================

% In-line function areatri
areatri = @(nodtri) .5*norm(cross([nodtri(3,:)-nodtri(2,:) 0], [nodtri(1,:)-nodtri(3,:) 0]));


TR = triangulation(tri,nod);
indtri = pointLocation(TR,points); 
indtri(isnan(indtri)) = 0; 
intu=zeros(size(points,1),1);

for k=1:size(points,1)
    if indtri(k)~=0
        point=points(k,:);
        tricell=tri(indtri(k),:);
        nodtri=nod(tricell,:);
        area=areatri(nodtri);
        % nodal shape functions (T3 elemns)
        N1=areatri([nodtri([2 3],:); point])/area;
        N2=areatri([nodtri([3 1],:); point])/area;
        N3=areatri([nodtri([1 2],:); point])/area;
        % interpolated potential value
        intu(k) = N1*u(tricell(1))+N2*u(tricell(2))+N3*u(tricell(3));
    end
end

end

function [intBr,intBz] = Bfield_axisymm_interp_T3(points,nod,tri,magflux)

%=====================================================================
% file: "Bfield_interp_T3.m"
% author: F. Moro
% affiliation: Universit� degli Studi di Padova, Italy
% date: 21-03-20
% note: compute magnetic flux density for 2-D axisymm. FEM with T3 elements
% input: query coord. matrix "points" (Px2), coord. matrix "nod" (Nx2), 
%        conn. matrix "tri" (Tx3), magnetic flux r.Atheta "magflux" (Nx1)
% output: B-field interpolated components "intBr" (Px1), "intBz" (Px1)
%=====================================================================

TR = triangulation(tri,nod);
indtri = pointLocation(TR,points); 
indtri(isnan(indtri)) = 0; 
grads = zeros(size(points,1),2); 

for k=1:size(points,1)
    if indtri(k)~=0
        tricell=tri(indtri(k),:); 
        nodtri=nod(tricell,:);
%         area=areatri(nodtri);
        area = .5 * norm(cross([nodtri(3,:)-nodtri(2,:) 0], [nodtri(1,:)-nodtri(3,:) 0]));
        centre=mean(nodtri,1); 
        % edge vectors
        L1 = (nodtri(3,:)-nodtri(2,:)); 
        L2 = (nodtri(1,:)-nodtri(3,:)); 
        L3 = (nodtri(2,:)-nodtri(1,:)); 
        % shape function gradients
        w1 = [-L1(2),L1(1)]/(2*area);
        w2 = [-L2(2),L2(1)]/(2*area);
        w3 = [-L3(2),L3(1)]/(2*area);
        grads(k,:) = w1*magflux(tricell(1))+w2*magflux(tricell(2))+...
            w3*magflux(tricell(3));
        % compute 1/r*\nabla(r.Atheta)
        grads(k,:) = grads(k,:)*1/centre(1); 
    end
end
% magnetic flux density: Br = -1/r*d(r.Atheta)/dz, Bz = 1/r*d(r.Atheta)/dtheta
intBr = -grads(:,2); intBz = grads(:,1); 

end
function [intEx,intEy] = Efield_interp_T3(points,nod,tri,u)

%=====================================================================
% file: "Efield_interp_T3.m"
% author: F. Moro
% affiliation: Universit� degli Studi di Padova, Italy
% date: 21-03-20
% note: compute electric field for 2-D FEM with T3 elements
% input: query coord. matrix points (Px2), coord. matrix nod (Nx2), 
%        conn. matrix tri (Tx3), potentials u (Nx1)
% output: E-field interpolated components intEx (Px1), intEy (Px1)
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
        % edge vectors
        L1 = (nodtri(3,:)-nodtri(2,:)); 
        L2 = (nodtri(1,:)-nodtri(3,:)); 
        L3 = (nodtri(2,:)-nodtri(1,:)); 
        % shape function gradients
        w1 = [-L1(2),L1(1)]/(2*area);
        w2 = [-L2(2),L2(1)]/(2*area);
        w3 = [-L3(2),L3(1)]/(2*area);
        grads(k,:) = w1*u(tricell(1))+w2*u(tricell(2))+w3*u(tricell(3));
    end
end
intEx = -grads(:,1); intEy = -grads(:,2); 

end

function Gmat = assemble_Gmat_BEM_2d(nod,bndedg)

%=====================================================================
% file: "assemble_Gmat_BEM_2d.m"
% author: F. Moro
% affiliation: Universit� degli Studi di Padova, Italy
% date: 14-05-20
% note: single-layer matrix assembly for solution
% input: nod (coordinate matrix: N x 3), bndedg (boundary edge matrix: E x 2)
% output: Gmat (single-layer matrix: E x E)
% ref. C.W.Trowbridge, "An Introduction to CAE Analysis", Vector Fields Ltd. 
%=====================================================================

Gmat = zeros(size(bndedg,1),size(bndedg,1)); 

% gauss points in the interval [-1,+1]
Ngauss = 4; % number of Gauss points
a = -1; b = 1; % interval end points
[xi,wi]=f_lgwt(Ngauss,a,b); % nodes and weights

% self-terms
for k=1:size(bndedg,1)
    edge = bndedg(k,:); % any row is a pair of nodal indices
    edgvec = nod(edge(2),:)-nod(edge(1),:); % edge vector 
    Ledg = norm(edgvec); 
    Lhalf = 0.5*Ledg; % half-length
    % collocation point is at the element centre
    Gmat(k,k) = Lhalf/pi*(1-log(Lhalf)); 
end

% mutual-terms 
for ii = 1:size(bndedg,1)
   edge_ii = bndedg(ii,:); 
   centre_ii = 0.5*(nod(edge_ii(1),:)+nod(edge_ii(2),:)); 
   for jj = 1:size(bndedg,1) 
       if ii ~= jj
           % oriented edge from A to B 
           edge_jj = bndedg(jj,:); 
           edgvec = nod(edge_jj(2),:)-nod(edge_jj(1),:); 
           % note: counter-clockwise oriented boundary; unit normal n points 
           % inwardly of computational domain; rotation in counterclockwise direct.   
           Ledg = norm(edgvec); 
           Lhalf = 0.5*Ledg; % half-length
           % note: xg = (1-xi)/2*A+(1+xi)/2*B  
           gpoints = (1-xi)/2.*repmat(nod(edge_jj(1),:),Ngauss,1)+...
               (1+xi)/2.*repmat(nod(edge_jj(2),:),Ngauss,1);
           distvecs = gpoints-repmat(centre_ii,Ngauss,1); 
           dists = sqrt(sum(distvecs.^2,2));
           Gmat(ii,jj) = Lhalf/(2*pi)*sum(wi.*log(1./dists)); 
       end
   end
end

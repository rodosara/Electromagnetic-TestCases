function [Fvec,Gmat] = postpro_BEM_2d_zeropot(nod,tri,bndedg,fun)

%=====================================================================
% file: "postpro_BEM_2d_zeropot.m"
% author: F. Moro
% affiliation: Universit� degli Studi di Padova, Italy
% date: 14-05-20
% note: source and single-layer matrix assembly for post-processing
% input: nod (coordinate matrix: N x 3), tri (element matrix: T x 3), 
%        bndedg (boundary edge matrix: E x 2), fun (source density function)    
% output: Fvec (source vector: E x 1), Gmat (single-layer matrix: E x E)
% ref. C.W.Trowbridge, "An Introduction to CAE Analysis", Vector Fields Ltd. 
%=====================================================================

Fvec = zeros(size(nod,1),1); 
Gmat = zeros(size(nod,1),size(bndedg,1)); 

% gauss points in the interval [-1,+1]
Ngauss = 4; % number of Gauss points
a = -1; b = 1; % interval end points
[xi,wi]=f_lgwt(Ngauss,a,b); % nodes and weights

% interior nodes
bndnods = unique(bndedg(:)); 
intnods = setdiff(1:size(nod,1),bndnods); 

% mutual-terms 
for ii = 1:length(intnods)
   nod_ii = nod(intnods(ii),:); 
   for jj = 1:size(bndedg,1) 
       % oriented edge from A to B 
       edge_jj = bndedg(jj,:); 
       edgvec = nod(edge_jj(2),:)-nod(edge_jj(1),:); 
       Ledg = norm(edgvec); 
       Lhalf = 0.5*Ledg; % half-length
       % note: xg = (1-xi)/2*A+(1+xi)/2*B  
       gpoints = (1-xi)/2.*repmat(nod(edge_jj(1),:),Ngauss,1)+...
           (1+xi)/2.*repmat(nod(edge_jj(2),:),Ngauss,1);
       vecs = gpoints-repmat(nod_ii,Ngauss,1); 
       distances = sqrt(sum(vecs.^2,2));
       Gmat(intnods(ii),jj) = Lhalf/(2*pi)*sum(wi.*log(1./distances)); 
   end
   % source term
   for m=1:size(tri,1)
       cell_m = tri(m,:);
       nodtri=nod(cell_m,:);
%        area=areatri(nodtri);
       area = .5 * norm(cross([nodtri(3,:)-nodtri(2,:) 0], [nodtri(1,:)-nodtri(3,:) 0]));
       % nine-point quadrature formula
       gaussnods=1/18*[1 4 13; 4 1 13; 1 13 4; 13 1 4; 4 7 7;...
           7 4 7; 7 7 4; 4 13 1; 13 4 1]*nodtri;
       distvecs = repmat(nod_ii,9,1)-gaussnods;
       dists = sqrt(sum(distvecs.^2,2));
       weight=area/9;
       % integration on the whole domain
       Fvec(intnods(ii)) = Fvec(intnods(ii))+weight/(2*pi)*sum(log(1./dists).*fun(gaussnods(:,1),gaussnods(:,2)));
   end
end

function Afield = f_Afield_segment(iph,nodA,nodB,points)

% =========================================================
% author: Federico Moro (federico.moro@unipd.it)
% affiliation: Dip. Ing. Industriale, Universit� di Padova
% date: 20-04-18
% =========================================================

% magnetic vector potential for a thin wire (oriented from A to B)
% formula: A = mu*I/(4*pi)*(asinh(BP/dist)+sin(PA/dist))*u
% u : wire direction
% P : field point projection on segment 
% dist : distance between field point and P

% iph: current phasor (reference from nodA to nodB)
% nodA, nodB (3x1): segment end point coordinates
% points (3xN): field point coordinates
% Afield (3xN): magnetic vector potential components computed in points

normcol = @(X) sqrt(sum(X.*X,1)); % compute norms along matrix columns
mu0 = 4*pi*1e-7; % vacuum magnetic permeability (H/m)

% field computed in all field points simultaneously
% P is the projection of field point on segment AB 

N = size(points,2); 
vecA = repmat(nodA,1,N)-points; 
vecB = repmat(nodB,1,N)-points; 
vecAB = repmat(nodB-nodA,1,N); 
uAB = vecAB./normcol(vecAB); 
PA = -dot(uAB,vecA,1); % projections along AB
BP = dot(uAB,vecB,1); 
dist = sqrt(normcol(vecA).^2-PA.^2); % segment distance from field point
Anorm = mu0*iph./(4*pi).*(asinh(BP./dist)+asinh(PA./dist));
Afield = repmat(Anorm,3,1).*uAB;

end

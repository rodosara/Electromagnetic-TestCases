function Bfield = f_Bfield_segment(iph,nodA,nodB,points)

% =========================================================
% author: Federico Moro (federico.moro@unipd.it)
% affiliation: Dip. Ing. Industriale, Universit� di Padova
% date: 20-04-18
% =========================================================

% magnetic flux density for a thin wire (oriented from A to B)
% formula: B = mu*I/(4*pi*r)*(sin(theta_A)+sin(theta_B))*u_phi
% r : wire distance from the field point
% theta_A, theta_B : angle at segment ends A, B

% iph: current phasor (reference from nodA to nodB)
% nodA, nodB (3x1): segment end point coordinates
% points (3xN): field point coordinates
% Bfield (3xN): magnetic flux density components computed in points

normcol = @(X) sqrt(sum(X.*X,1)); % compute norms along matrix columns
mu0 = 4*pi*1e-7; % vacuum magnetic permeability (H/m)

% field computed in all field points simultaneously
% P is the projection of field point on segment AB 

N = size(points,2); 
vecA = repmat(nodA,1,N)-points; 
vecB = repmat(nodB,1,N)-points; 
vecAB = repmat(nodB-nodA,1,N); 
uAB = vecAB./normcol(vecAB); 
tA = vecA./normcol(vecA); 
vec_phi = cross(tA,uAB,1); % azimuthal vectors
u_phi = vec_phi./normcol(vec_phi); 
PA = dot(uAB,-vecA,1); % projections along AB
BP = dot(uAB,vecB,1); 
dist = sqrt(normcol(vecA).^2-PA.^2); % segment distance from field point
Bnorm = mu0*iph./(4*pi*dist).*(BP./normcol(vecB)+PA./normcol(vecA));
Bfield = repmat(Bnorm,3,1).*u_phi;

end
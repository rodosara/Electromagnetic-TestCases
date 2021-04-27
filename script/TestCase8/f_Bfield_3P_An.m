%% FUNCTION SECTION
function Bfield_3P_AN = f_Bfield_3P_An(mu0, Yph, I, field_points, normcol)

% Distance between cable phase 1 and 2
d12 = norm(Yph(:,1)-Yph(:,2));
% Distance between cable phase 2 and 3
d23 = norm(Yph(:,2)-Yph(:,3));
% Distance between cable phase 3 and 1
d31 = norm(Yph(:,3)-Yph(:,1));

% Radius from filed point x and source pint phase 1
r1 =  normcol(field_points - repmat(Yph(:,1), 1, length(field_points)));
% Radius from filed point x and source pint phase 2
r2 =  normcol(field_points - repmat(Yph(:,2), 1, length(field_points)));
% Radius from filed point x and source pint phase 3
r3 =  normcol(field_points - repmat(Yph(:,3), 1, length(field_points)));

% Analytical formula fro 3 phases system
Bfield_3P_AN = mu0*I/(2*sqrt(2)*pi) * sqrt((d12./(r1.*r2)).^2 + (d23./(r2.*r3)).^2 + (d31./(r3.*r1)).^2);
end

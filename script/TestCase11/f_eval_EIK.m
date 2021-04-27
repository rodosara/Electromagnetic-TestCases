% --------------------------------------------------------------
% Evaluate the complete elliptic integral of the first kind
% by means of a polynomial approximation [M Abramowitz and
% I A Stegun, Handbook of Mathematical Functions, National
% Bureau of Standards, 1965]
% --------------------------------------------------------------
function EIK = eval_EIK(x)
% Arguments:
% x = argument for K(x) in the interval 0 <= x < 1
% Returns:
% EIK = the value of the complete elliptic integral of
% the first kind (with an error less than 2e-8)
a = [0.01451196212; 0.03742563713; 0.03590092383; ...
0.09666344259; 1.38629436112];
b = [0.00441787012; 0.03328355346; 0.06880248576; ...
0.12498593597; 0.50000000000];
m1 = 1 - x;
EIK = polyval(a,m1) - polyval(b,m1).*log(m1);
end

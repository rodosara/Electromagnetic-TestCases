function f_Mph_check(M_comp, L, d, mu0)
% mutual-inductance between parallel wires (dist >> radius)
% ref. C.R. Paul, "Inductance: Loop & Partial", Wiley 2009 (p. 211,
% eq. 5.21b)

M_an = zeros(3,3); % analytical

mm = @(dist) mu0/(2*pi)*L*(asinh(L/dist)-sqrt(1+(dist/L)^2)+(dist/L)); 

M_an(1,2) = mm(d); M_an(2,1) = M_an(1,2); 
M_an(1,3) = mm(2*d); M_an(3,1) = M_an(1,3);
M_an(2,3) = mm(d); M_an(3,2) = M_an(2,3); 

% Validation
disp('Discrepancy between analytical and numerical mutual-inductance matrices')
err = norm(M_an-M_comp)/norm(M_an)*100

end


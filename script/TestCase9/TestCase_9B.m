%% TestCase 9B - PEEC
clear all;
close all;
clc;

%% PRE-PROCESSING
% RMS current for each phase
I = 100; % [A]
% Frequency
freq = 50; % [Hz]
% Distance between conductor
d = 20e-2; % [m]
% Distance between conductor and measurement line
h = .5; % [m]
% Length of conductors
L = 2; % [m]
% Phase centre with coordinate in 3 axis X, Y and Z
Yph = [[-d 0 d]; zeros(1,3); zeros(1,3)];
% Current phasor
I_phasor = I*exp(1i*[-2*pi/3 0 2*pi/3]);
% Vacuum magnetic permeability
mu0 = 4*pi*1e-7; % [H/m]
% Radius normalize
r0 = 1; % [m]
% Resistivity
ro = 20e-9; % [ohm*m]
% Cross radius of conductors
r = 2e-2; % [m]
% Cross section of conductors
S = pi*r^2; % [m^2]

% Number of edges for each busbar
Ne = 100;

%% PROCESSING
% Discretized measurement line
x = linspace(-2, +2, 1000); % [m]
% Vector of position along measurement line
field_points = [x; h*ones(1, length(x)); zeros(1, length(x))];

% Compute nodal ((3x3(N+1)) matrix contains nodal coordinates) and edges (
% (3(N+1)x2) matrix contain exdges nodes)
[nod, edg, ph_edg] = f_discretize_busbar(L, d, Ne);

% R matrix resistivity of all elements
% Partial resistance (consider equal length of all elements)
Ri = ro*(L/Ne)/S;
% Ri = ro*norm(nod(:,1)+nod(:,2))/S;
R_matrix = spdiags(Ri*I*ones(3*Ne,1), 0, 3*Ne, 3*Ne);

% M mutual inductance matrix (elements in different phase conductors)
M_matrix = sparse(Ne,Ne);
for ii = 1:length(edg)
    % Node for ii element
    nodA_ii = nod(:,edg(ii,1));
    nodB_ii = nod(:,edg(ii,2));
    for jj = ii:length(edg)
        % Node for jj element
        nodA_jj = nod(:,edg(jj,1));
        nodB_jj = nod(:,edg(jj,2));
        % Center point of element g
        point_jj = (nodA_jj+nodB_jj)/2;
        % Vector
        edg_vec_jj = nodB_jj-nodA_jj;
        % Verify elements below different phase
        if ph_edg(ii) ~= ph_edg(jj)
            % Compute a vector potential with 1 A of current
            Apot_ii = f_Afield_segment(1, nodA_ii, nodB_ii, point_jj);
            % Calc numerical aprozimation of analytical formula
            M_matrix(ii,jj) = dot(Apot_ii,edg_vec_jj);
            % Symmetric matrix
            M_matrix(jj,ii) = M_matrix(ii,jj);
        end
    end
end

% L_wire self-inductance matrix (elements in same phase conductors)
% Aspect radio radius/length used for numerical aproximation
a_r = r/L;
L_wire = mu0/(2*pi)*L*(asinh(1/a_r)-sqrt(1+a_r^2)+a_r+1/4);

% C incident matrix (3Nx3)
C_matrix = kron(speye(3),ones(Ne,1));

% PH matrix for each phase
Rph = C_matrix'*R_matrix*C_matrix;
Lph = L_wire*eye(3) + C_matrix'*M_matrix*C_matrix;
I_ph = I_phasor.';

% Check the matrix respect the analytical
Mph = C_matrix'*M_matrix*C_matrix;
f_Mph_check(Mph, L, d, mu0);

% Compute Uph vector of potential's phasor for each phase
Uph = (Rph + 1i*2*pi*freq*Lph)*I_ph

%% POST-PROCESSING
% Uph phasor plot
figure('Name', 'U and I phasor for each phase');
subplot(1,2,1);
for ii = 1:3
    hold on;
    quiver(0, 0, real(Uph(ii)), imag(Uph(ii)), 'Color', 'b');
end
hold off;
title('U phasor for each phase');
legend ({'Voltage U' 'Voltage V' 'Voltage W'});
grid on;
axis equal;
xlabel('Re[U]');
ylabel('Imag[U]');

subplot(1,2,2);
for ii = 1:3
    hold on;
    quiver(0, 0, real(I_phasor(ii)), imag(I_phasor(ii)), 'Color', 'r');
end
hold off;
title('I phasor for each phase');
legend ({'Current U' 'Current V' 'Current W'});
grid on;
axis equal;
xlabel('Re[I]');
ylabel('Imag[I]');
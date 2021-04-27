%% TestCase 11 - MoM
clear all;
close all;
clc;

%% PRE-PROCESSING
% Speed of light 
c_const = 3.0e+08; % [m/s]
% Intrinsic impedance of free-space
etha = 377; % [ohm]
% Operating frequency
freq = 3.0e+08; % [Hz]
% Applied voltage 
volt = 1; % [V]
% Field wavelength 
lamda = c_const/freq; % [m]
% Wire length
L_rel = 0.50; % [m]
% Antenna length in wavelengths
L = L_rel*lamda; 
% Wire radius
A_rel = 0.005; % [m]
% Antenna radius in wavelengths
A = A_rel*lamda;
% Free-space wavenumber
k = 2*pi/lamda; % [m^-1]
% Number of segments (center element is odded: 0-th segment centered in the middle)
N = 21;
% N = 201; % approximate kernel solution is unstable!!
% Segment length
delta = L/N; % [m]

%% PROCESSING
% Hallen equation

% Number of Gauss’s points
Nint = 31;
% Number of current samples on the upper-half of the antenna
M = (N-1)/2; % N = 2*M+1 (N is odd!)
% Hallen equation solving for benchmark
[Iz_orf,Z,cnd] = f_hallen_orfanidis(L_rel, A_rel, M, Nint, 1);

% Green's kernel methods

% Antenna discretization
% Collocation or mesh points
Zm = delta*(-M:M);

C_vector = cos(k*Zm)';
S_vector = sin(k*abs(Zm))';

Z_app = zeros(N,N);
Z_acc = zeros(N,N);
for n = 1:N
    for m = 1:N
    % Aproximate Green's kernel inline function
    R = @(z) sqrt(A.^2+(z-Zm(m)).^2);
    % Radius inline function
    G_app = @(R) exp(-1i*k*R)./R;
    
    % Accurate Green's kernel
    xi = @(z) Zm(m)-z;
    % Parametric coordinate xi inline function
    G_acc = @(xi) 2/pi./sqrt(xi.^2+4*A^2)... % eval_EIK external file function for elliptic integral K
        .*f_eval_EIK(4*A^2./(xi.^2+4*A^2))...% elliptic integral
        +(exp(-1i*k*sqrt(xi.^2+A^2))-1)./sqrt(xi.^2+A^2);

    % Integral approximate inline function
    Gint_app = integral(@(z) G_app(R(z)), Zm(n)-delta/2, Zm(n)+delta/2, ...
        'AbsTol', 1e-6, 'RelTol', 1e-3);
    
    % Integral accurate inline function
    Gint_acc = integral(@(z) G_acc(xi(z)), Zm(n)-delta/2, Zm(n)+delta/2, ...
        'AbsTol', 1e-6, 'RelTol', 1e-3);
    % Impedence matrix Z assemblying
    Z_app(n,m) = 1i*etha/(2*pi) * Gint_app;
    Z_acc(n,m) = 1i*etha/(2*pi) * Gint_acc;
    
    end
end

% Boundary condition - impose zero current at both ends of antenna
u = zeros(1,N);
u([1 end]) = 1;

% Costant C1 computation without inverse
C1_app = -volt*((u*(Z_app\S_vector))/(u*(Z_app\C_vector)));
C1_acc = -volt*((u*(Z_acc\S_vector))/(u*(Z_acc\C_vector)));

% Current density vector computation
Iz_app = Z_app \ (C1_app*C_vector + volt*S_vector);
Iz_acc = Z_acc \ (C1_acc*C_vector + volt*S_vector);

%% POST-PROCESSING
% Solution comparison plot
figure('Name', 'Comparison plot');
hold on;
plot(-M:M, real(Iz_acc), 'k-', -M:M, imag(Iz_acc), 'r-');
plot(-M:M, real(Iz_app), 'k--', -M:M, imag(Iz_app), 'r--');
title('Comparison plot');
hold off;
xlabel('Segment number');
ylabel('Current [A]');
legend('Accurate Re[Iz]', 'Accurate Im[Iz]', 'Approximate Re[Iz]', 'Approximate Im[Iz]');

% Orfanidi's book solution plot
figure('Name', 'Orfanidi solution');
plot(-M:M, real(Iz_orf), 'k-', -M:M, imag(Iz_orf), 'r-');
title('Orfanidi solution');
xlabel('Segment number');
ylabel('Current [A]');
legend('Orfanidi Re[Iz]', 'Orfanidi Im[Iz]');
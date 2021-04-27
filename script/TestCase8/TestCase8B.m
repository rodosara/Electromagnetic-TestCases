%% TestCase 8B - BIOT-SAVART LAW with finite conductors
clear all;
close all;
clc;

%% PRE-PROCESSING
% RMS current for each phase
I = 100; % [A]
% Distance between conductor
d = 20e-2; % [m]
% Distance between conductor and measurement line
h = .5; % [m]
% Length of conductors
L = 2; % [m]
% Phase centre with coordinate in 3 axis X, Y and Z
Yph = [[-d 0 d]; zeros(1,3); zeros(1,3)];
% Current phasor
Iph = I*exp(1i*[-2*pi/3 0 2*pi/3]);
% Vacuum magnetic permeability
mu0 = 4*pi*1e-7; % [H/m]

% Discretized measurement line
x = linspace(-2, +2, 100); % [m]
% Vector of position along measurement line
points = [x; h*ones(1, length(x)); zeros(1, length(x))];

% Start point of busbar
A = [Yph(1,:); Yph(2,:); (L/2)*ones(1,3)];
% End point of busbar
B = [Yph(1,:); Yph(2,:); -(L/2)*ones(1,3)];

Bfield = zeros(3,length(points));
for ii = 1:3
    b_field = f_Bfield_segment(Iph(ii), A(:,ii), B(:,ii), points);
    Bfield = Bfield + b_field;
end

%% PROCESSING
% Function for compute the norm for each column of a vector
normcol = @(X) sqrt(sum(X.*X,1));
% Compute magnitude of vector by each point (by column)
Bfield = abs(Bfield);
Bfield = normcol(Bfield);

% Re-compute Bfield with L=10 [m] to simulate infinite lenght
L = 10; % [m]
% Start point of busbar
A = [Yph(1,:); Yph(2,:); (L/2)*ones(1,3)];
% End point of busbar
B = [Yph(1,:); Yph(2,:); -(L/2)*ones(1,3)];

% Bfield computation
Bfield_10 = zeros(3,length(points));
for ii = 1:3
    b_field = f_Bfield_segment(Iph(ii), A(:,ii), B(:,ii), points);
    Bfield_10 = Bfield_10 + b_field;
end

% Compute magnitude of vector by each point (by column)
Bfield_10 = abs(Bfield_10);
Bfield_10 = normcol(Bfield_10);


% Function to compute analytical solution with formula
Bfield_3P_AN = f_Bfield_3P_An(mu0, Yph, I, points, normcol);

%% POST-PROCESSING
% Plot B field solution comparison
figure('Name', 'B field solution comparison');
hold on;
% Plot computed solution L = 2 [m]
plot(x, Bfield, '--');
title('B field solution comparison');
% Plot analytical solution L = 10 [m]
plot(x, Bfield_3P_AN, 'b', 'LineWidth', 1.5);
% Plot computed solution L = 10 [m]
scatter(x, Bfield_10, 'r');
% Plot phase conductors
scatter(Yph(1,:), Yph(2,:), 'o', 'MarkerFaceColor', 'r');
hold off;
legend ({'Analytical solution L=2 [m]' 'Analytical solution infinite' 'Computed solution L=10 [m]' ' Phase conductors'});
grid on;
xlabel('x [m]');
ylabel('B_{RMS} [T]');

% FEM 2dim plot B field solution comparison
% Import COMSOL output
comsolBfield = readtable('./comsol_Bfield.txt', 'HeaderLines', 8);
comsolBfield = table2array(comsolBfield);

figure('Name', 'FEM 2dim B field solution comparison');
hold on;
% Plot FEM 2dim with COMSOL solution L = 10 [m]
plot(comsolBfield(:,1), comsolBfield(:,2), 'b');
title('FEM 2dim B field solution comparison');
% Plot analytical solution L = 10 [m]
scatter(x, Bfield_3P_AN, 'k', 's');
% Plot computed solution L = 10 [m]
scatter(x, Bfield_10, 'r');
% Plot phase conductors
scatter(Yph(1,:), Yph(2,:), 'o', 'MarkerFaceColor', 'r');
hold off;
legend ({'FEM 2dim COMSOL solution L=10 [m]' 'Analytical solution infinite' 'Computed solution L=10 [m]' ' Phase conductors'});
grid on;
xlabel('x [m]');
ylabel('B_{RMS} [T]');

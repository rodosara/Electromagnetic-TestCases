%% TestCase 8A - BIOT-SAVART LAW with infinite conductors
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

%% PROCESSING
% Function for compute the norm for each column of a vector
normcol = @(X) sqrt(sum(X.*X,1));

% Iterate for each phase
Bfield = zeros(3,length(points));
for ii = 1:3
    % Vector position for each field points
    vec = [x-repmat(Yph(1,ii), 1, length(points)); h*ones(1,length(points))-repmat(Yph(2,ii), 1, length(points)); zeros(1,length(points))];
    % Versor vector position
    e_r = vec./normcol(vec);
    % Versor ez according current sense
    e_z = repmat([0; 0; -1], 1, length(points));
    % Versor e phi
    e_phi = cross(e_z,e_r,1);
    % Biot-Savart law in a 3-dim matrix
    b_field = mu0/(2*pi) * Iph(ii)./normcol(vec) .* e_phi;
    
    % Sum each phase
    Bfield = Bfield + b_field;
end

% Compute magnitude of vector by each point (by column)
Bfield = abs(Bfield);
Bfield = normcol(Bfield);

% Function to compute analytical solution with formula
Bfield_3P_AN = f_Bfield_3P_An(mu0, Yph, I, points, normcol);

%% POST-PROCESSING
% Plot B field solution comparison
figure('Name', 'B field solution comparison');
hold on;
% Plot analytical solution
plot(x, Bfield_3P_AN, 'b', 'DisplayName', 'Analytical solution', 'LineWidth', 1.5);
title('B field solution comparison');
% Plot computed solution
scatter(x, Bfield, 'Marker', 'o', 'LineWidth', 1.5);
% Plot phase conductors
scatter(Yph(1,:), Yph(2,:), 'o', 'MarkerFaceColor', 'r');
hold off;
legend ({'Analytical solution' 'Computed solution' ' Phase conductors'});
grid on;
xlabel('x [m]');
ylabel('B_{RMS} [T]');

% Plot vector behavoir
figure('Name', 'Vector');
color = ['m' 'g' 'b'];
for ii = 1:3
    plot([Yph(1,ii) points(1,1)], [Yph(2,ii) points(2)], 'Color', color(ii));
    hold on;
    scatter(Yph(1,ii), Yph(2,ii), color(ii));
    vec = [x-repmat(Yph(1,ii), 1, length(points)); h*ones(1,length(points))-repmat(Yph(2,ii), 1, length(points)); zeros(1,length(points))];
    % Versor vector position
    e_r = vec./normcol(vec);
    % Versor ez according current sense
    e_z = repmat([0; 0; -1], 1, length(points));
    % Versor e phi
    e_phi = cross(e_z,e_r,1);
    quiver(points(1,1:10:end), points(2,1:10:end), e_phi(1,1:10:end), e_phi(2,1:10:end), 'Color', color(ii), 'AutoScale', 'on');
    axis = 'equal';
end
plot([points(1,1) points(1,end)], [points(2,1) points(2,end)], 'r', 'LineWidth', 1.5);
scatter(points(1,1:10:end), points(2,1:10:end), 'k', '+', 'LineWidth', 1.5);
hold off;
legend ({'Vector r' 'Conductor' 'Versor phi' 'Measurement line' 'Field points'});
xlabel('x [m]');
ylabel('y [m]');

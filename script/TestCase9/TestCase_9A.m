%% TestCase 9A - MT-BT BUS BAR
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
% Radius normalize
r0 = 1; % [m]

% Discretized measurement line
x = linspace(-2, +2, 1000); % [m]
% Vector of position along measurement line
field_points = [x; h*ones(1, length(x)); zeros(1, length(x))];

% Start point of busbar
A = [Yph(1,:); Yph(2,:); (L/2)*ones(1,3)];
% End point of busbar
B = [Yph(1,:); Yph(2,:); -(L/2)*ones(1,3)];

Afield = zeros(3,length(field_points));
for ii = 1:3
    a_vector = f_Afield_segment(Iph(ii), A(:,ii), B(:,ii), field_points);
    Afield = Afield + a_vector;
end

%% PROCESSING
% Function for compute the norm for each column of a vector
normcol = @(X) sqrt(sum(X.*X,1));
% Compute magnitude of vector by each point (by column)
Afield = abs(Afield);
Afield = normcol(Afield);

% Recompute the A vector with analytical formula
% using L=10 [m] to simulate the infinite length
L = 10; % [m]

% Start point of busbar
A = [Yph(1,:); Yph(2,:); (L/2)*ones(1,3)];
% End point of busbar
B = [Yph(1,:); Yph(2,:); -(L/2)*ones(1,3)];

Afield_10 = zeros(3,length(field_points));
for ii = 1:3
    a_vector = f_Afield_segment(Iph(ii), A(:,ii), B(:,ii), field_points);
    Afield_10 = Afield_10 + a_vector;
end

% Function for compute the norm for each column of a vector
normcol = @(X) sqrt(sum(X.*X,1));
% Compute magnitude of vector by each point (by column)
Afield_10 = abs(Afield_10);
Afield_10 = normcol(Afield_10);

% Avector analytical for three phases
Afield_3P_AN = zeros(3,length(field_points));
for ii = 1:3
    a_vector = -mu0/(2*pi) * Iph(ii) * log(normcol(field_points-Yph(:,ii))./r0);
    a_vector = [zeros(1,length(field_points)); zeros(1,length(field_points)); a_vector];
    Afield_3P_AN = Afield_3P_AN + a_vector;
end
% Compute magnitude of vector by each point (by column)
Afield_3P_AN = abs(Afield_3P_AN);
Afield_3P_AN = normcol(Afield_3P_AN);

%% POST-PROCESSING
% B field solution comparison plot
figure('Name', 'A field solution comparison');
hold on;
% Plot computed solution
plot(x, Afield, '--');
plot(x, Afield_10, 'LineWidth', 1.5);
title('A field solution comparison');
% Plot analytical solution
plot(x, Afield_3P_AN, 'b');
% Plot phase conductors
scatter(Yph(1,:), Yph(2,:), 'o', 'MarkerFaceColor', 'r');
hold off;
legend ({'Computed solution L=2 [m]' 'Computed solution L=10 [m]' 'Analytical solution infinite' 'Phase conductors'});
grid on;
xlabel('x [m]');
ylabel('A_{RMS} [T]');
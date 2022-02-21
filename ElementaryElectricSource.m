close all;
clear;
clc;

%% Constants
% Index
N = 500;                        % Number of points
% Field
f = 30e9;                       % Frequency of source [Hz]
kx_max = 3;                     % Maximum x-component of the wave number
% Medium
er = 1;                         % Relative permittivity
c = physconst('LightSpeed');    % Speed of light [m/s]

%% Parameters
wlen = c / f;                   % Wavelength [m]
k0 = 2*pi / wlen;               % Magnitude of wave number [rad/m]

%% x and y-Components of Wave Number
kx = k0 * linspace(eps, kx_max, N);
ky = zeros( 1, length(kx) );
[ KX, KY ] = meshgrid(kx, ky);

%% Calculate Spectral Green's Function (SGF)
ej_SGF = EJ_SGF(er, k0, KX, KY);

%% Plot x-Component of SGF
figure();
plot(kx / k0, real( ej_SGF(1, :, 1, 1) ), 'LineWidth', 3);
hold on;
plot(kx / k0, imag( ej_SGF(1, :, 1, 1) ), 'LineWidth', 3);
grid on;
xlabel('k_{x} / k_{0}');
ylabel('G^{ej}_{xx}');
legend('Real', 'Imag');
xticks((0: 0.25 : 3));
yticks((-200 : 100 : 600));

%% Plot y-Component of SGF
figure();
plot(kx / k0, real( ej_SGF(1, :, 2, 2) ), 'LineWidth', 3);
hold on;
plot(kx / k0, imag( ej_SGF(1, :, 2, 2) ), 'LineWidth', 3);
grid on;
xlabel('k_{x} / k_{0}');
ylabel('G^{ej}_{yy}');
legend('Real', 'Imag');
ylim([-700 100]);
xticks((0: 0.25: 3));
yticks((-700 : 100: 100));

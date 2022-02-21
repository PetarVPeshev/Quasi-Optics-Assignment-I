close all;
clear;
clc;

%% Constants
% Index
N = 500;                        % Number of points
P0_indx = 1;                    % Index corresponding to Phi 0 degrees
P45_indx = 64;                  % Index corresponding to Phi 45 degrees
P90_indx = 126;                 % Index corresponding to Phi 90 degrees
% Field
f = 30e9;                       % Frequency of source [Hz]
R = 5;                          % Radial distance [m]
kx_max = 3;                     % Maximum x-component of the wave number
% Dipole
h = 7.5e-3;                     % Height of dipole [m]
% Height sweep
hmin = 0e-3;                    % Min height of dipole [m]
hmax = 15e-3;                   % Max height of dipole [m]
% Medium
er = 1;                         % Relative permittivity
c = physconst('LightSpeed');    % Speed of light [m/s]

%% Parameters
wlen = c / f;                   % Wavelength [m]
k0 = 2*pi / wlen;               % Magnitude of wave number [rad/m]
l = wlen / 2;                   % Length of antenna dipole [m]
w = wlen / 40;                  % Width of antenna dipole [m]

%% Height Sweep
hh = linspace(hmin, hmax, N);

%% Theta and Phi-Components of Spherical Coordinates
th = linspace(0, pi / 2, N);
dth = th(2) - th(1);
ph = linspace(0, 2 * pi, N);
dph = ph(2) - ph(1);
[ TH, PH ] = meshgrid(th, ph);

%% x, y, and z-Components of Wave Number
kx = k0 * linspace(eps, kx_max, N);
ky = zeros( 1, length(kx) );
kz = -1j * sqrt( -(k0^2 - kx.^2 - ky.^2) );
[ KX, KY ] = meshgrid(kx, ky);

%% Calculate Spectral Green's Function (SGF)
ej_SGF = EJ_SGF(er, k0, KX, KY);

%% Calculate Fourier Transform (FT) of Current Distribution
Jx = FTCurrent(k0, er, KX, KY, l, w) .* 2j .* sin(kz * h);

%% Calculate Electric Far-Field
E = farfield(k0, R, TH, PH, kz, ej_SGF, Jx);
Et = sqrt( abs( E(:, :, 2) ).^2 + abs( E(:, :, 3) ).^2 );

%% Plot the Total Electric Field in Phi 0, 45, 90 degree planes
figure();
plot(th * 180/pi, 20*log10( Et(P0_indx, :) ) - ...
     max( 20*log10( Et(P0_indx, :) ) ), 'LineWidth', 3);
hold on;
plot(th * 180/pi, 20*log10( Et(P45_indx, :) ) - ...
     max( 20*log10( Et(P45_indx, :) ) ), 'LineWidth', 3);
hold on;
plot(th * 180/pi, 20*log10( Et(P90_indx, :) ) - ...
     max( 20*log10( Et(P90_indx, :) ) ), 'LineWidth', 3);
grid on;
legend(['\phi = 0' char(176)], ['\phi = 45' char(176)], ...
       ['\phi = 90' char(176)]);
xlabel('\theta [deg]');
ylabel('E_{t} [dB]');
ylim([-40 0]);

%% Interpolate Data for Smoother Surface
th_inter = linspace(0, pi / 2, N * 5);
ph_inter = linspace(0, pi, N * 5);
[ TH_inter, PH_inter ] = meshgrid(th_inter, ph_inter);
Eth_inter = interp2(TH, PH, E(:, :, 2), TH_inter, PH_inter, 'spline');

%% UV Representation Coordinates
U = sin(TH_inter) .* cos(PH_inter);
V = sin(TH_inter) .* sin(PH_inter);

%% Plot UV Representation of Far-Field
figure();
surface(U, V, 20*log10( abs(Eth_inter) ) - ...
        max( max( 20*log10( abs(Eth_inter) ) ) ), 'LineStyle', 'none' );
grid on;
colormap('jet');
colorbar;
caxis([-10, 0]);
view(-37.5, 30);
xlabel('U');
ylabel('V');
zlabel('|E_{\theta}| [dB]');
zlim([-10 0]);
xticks((-1 : 0.25 : 1));
yticks((0 : 0.25 : 1));

%% Evaluate Directivity as a Function of Frequency at Broadside
Dir = zeros(1, length(hh));
for i = 1 : length(hh)   
    % Calculate Fourier Transform (FT) of Current Distribution
    Jx_h = FTCurrent(k0, er, KX, KY, l, w) .* 2j .* sin( kz * hh(i) );

    % Calculate Electric Far-Field
    E_h = farfield(k0, R, TH, PH, kz, ej_SGF, Jx_h);
    Et_h = sqrt( abs( E_h(:, :, 2) ).^2 + abs( E_h(:, :, 3) ).^2 );

    % Calculate Directivity
    [ TDir, ~ ] = Directivity(Et_h, TH, dth, dph, er, R);
    Dir(i) = TDir(1, 1);
end

%% Plot Directivity at Broadside as Function of Frequency
figure();
plot(hh * 1e3, 10*log10( Dir ) - max( 10*log10( Dir ) ), 'LineWidth', 3);
grid on;
xlabel('h [mm]');
ylabel(['D( \theta = 0' char(176) ' , \phi = 0' char(176) ' )']);
ylim([-40 0]);
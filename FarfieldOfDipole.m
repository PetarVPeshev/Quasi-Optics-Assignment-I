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
% Frequency sweep
fmin = 20e9;                     % Min frequency for directivity [Hz]
fmax = 40e9;                    % Max frequency for directivity [Hz]
% Medium
er = 1;                         % Relative permittivity
c = physconst('LightSpeed');    % Speed of light [m/s]

%% Parameters
wlen = c / f;                   % Wavelength [m]
k0 = 2*pi / wlen;               % Magnitude of wave number [rad/m]
l = wlen / 2;                   % Length of antenna dipole [m]
w = wlen / 40;                  % Width of antenna dipole [m]

%% Frequency Sweep
ff = linspace(fmin, fmax, N);

%% Theta and Phi-Components of Spherical Coordinates
th = linspace(eps, pi, N);
dth = th(2) - th(1);
ph = linspace(eps, 2 * pi, N);
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
Jx = FTCurrent(k0, er, KX, KY, l, w);

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
th_inter = linspace(eps, pi, N * 5);
ph_inter = linspace(eps, 2 * pi, N * 5);
[ TH_inter, PH_inter ] = meshgrid(th_inter, ph_inter);
Eth_inter = interp2(TH, PH, E(:, :, 2), TH_inter, PH_inter, 'spline');
Ephi_inter = interp2(TH, PH, E(:, :, 3), TH_inter, PH_inter, 'spline');
E_inter = sqrt( abs( Eth_inter ).^2 + abs( Ephi_inter ).^2 );

%% UV Representation Coordinates
U = sin(TH_inter) .* cos(PH_inter);
V = sin(TH_inter) .* sin(PH_inter);

%% Plot UV Representation of Far-Field
figure();
surface(U, V, 20*log10( abs(E_inter) ) - ...
        max( max( 20*log10( abs(E_inter) ) ) ), 'LineStyle', 'none' );
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
yticks((-1 : 0.25 : 1));

%% Evaluate Directivity as a Function of Frequency at Broadside
Dir = zeros(1, length(ff));
for i = 1 : length(ff)
    % Parameters
    wlen_f = c / ff(i);          % Wavelength [m]
    k0_f = 2*pi / wlen_f;     % Magnitude of wave number [rad/m]
    
    % x, y, and z-Components of Wave Number
    kx_f = k0_f * linspace(eps, kx_max, N);
    ky_f = zeros( 1, length(kx_f) );
    kz_f = -1j * sqrt( -(k0_f^2 - kx_f.^2 - ky_f.^2) );
    [ KX_f, KY_f ] = meshgrid(kx_f, ky_f);

    % Calculate Spectral Green's Function (SGF)
    ej_SGF_f = EJ_SGF(er, k0_f, KX_f, KY_f);

    % Calculate Fourier Transform (FT) of Current Distribution
    Jx_f = FTCurrent(k0_f, er, KX_f, KY_f, l, w);

    % Calculate Electric Far-Field
    E_f = farfield(k0_f, R, TH, PH, kz_f, ej_SGF_f, Jx_f);
    Et_f = sqrt( abs( E_f(:, :, 2) ).^2 + abs( E_f(:, :, 3) ).^2 );

    % Calculate Directivity
    [ TDir, ~ ] = Directivity(Et_f, TH, dth, dph, er, R);
    Dir(i) = TDir(1, 1);
end

%% Plot Directivity at Broadside as Function of Frequency
figure();
plot(ff * 1e-9, 10*log10( Dir ) - max( 10*log10( Dir ) ), 'LineWidth', 3);
grid on;
xlabel('f [GHz]');
ylabel(['D( \theta = 0' char(176) ' , \phi = 0' char(176) ' )']);
ylim([-11 0]);
xlim([20 40]);

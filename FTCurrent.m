function [ Jx ] = FTCurrent( k0, er, KX, KY, L, W )
%FTCurrent This function computes the Fourier Transform of the electric
%current
%   The function takes the relative permittivity of the medium, the
%   magnitude of the wave number, meshgrid of the x and y-components of the
%   wave number, length, and width of a dipole antenna as inputs.
%   It outputs the Fourier Transform (FT) of the electric current density
%   of the dipole antenna in 2D matrix corresponding to the x and
%   y-component meshgrid of the wave number.
    %% Extract x and y-Coordinates Wave Number Length
    kx_l = size(KX, 2);
    ky_l = size(KY, 1);
    %% Calculate keq
    keq = k0 * sqrt(er);
    %% Calculate FT
    T = sinc(KY * W / (2*pi));
    L = 2 * keq .* (cos(KX * L/2) - cos(keq * L/2)) ./...
        ( (keq.^2 - KX.^2) .* sin(keq * L/2) );
    Jx = zeros(kx_l, ky_l, 3);
    Jx(:, :, 1) = T .* L;
end
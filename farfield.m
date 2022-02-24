function [ E ] = farfield( k0, R_FF, TH, PH, KZ, ej_SGF, Jx )
%farfield This function computes the electric far-field of an antenna
%   The function takes the magnitude of the wave number, an array of the
%   radial distance, a meshgrid of the theta and phi spherical coordinates,
%   an array of the z-component of the wave number, a 4D matrix of the
%   Spectral Green's Function (SGF) of the electric field due to electric
%   current (dimensions 1 and 2 corresponding to a meshgrid of the x and
%   y-components of the wave number, and dimensions 3 and 4 representing
%   the dyadic tensor product), and a meshgrid of the Fourier Transform of
%   the electric current as inputs.
%   It outputs the electric field in spherical coordinates in a 4D matrix
%   with dimensions 1 and 2 corresponding to the theta and phi meshgrid,
%   and dimensions 3 and 4 representing the dyadic tensor product.
    %% Calculate Product of SGF and Spatial Current Distribution
    SGF_JX = zeros( size(ej_SGF) );
    for i = 1:3
        for n = 1:3
            SGF_JX(:, :, i, n) = 1j * KZ .* ej_SGF(:, :, i, n) .* Jx(:, :, n);
        end
    end
    %% Convert to Spherical Coordinates
    SGF_JX_sph = cart2sphereV(squeeze( SGF_JX(:, :, 1, :) ), TH, PH);
    %% Calculate Total Electric Field
    E = SGF_JX_sph .* exp(-1j * k0 * R_FF) ./ (2*pi * R_FF);
end
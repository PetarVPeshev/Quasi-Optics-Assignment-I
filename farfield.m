function [ E ] = farfield( k0, R_FF, TH, PH, KZ, ej_SGF, Jx )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    %% Calculate product of SGF and spatial current distribution
    SGF_JX = zeros( size(ej_SGF) );
    for i = 1:3
        for n = 1:3
            SGF_JX(:, :, i, n) = 1j * KZ .* ej_SGF(:, :, i, n) .* Jx;
        end
    end
    %% Convert to spherical coordinates
    SGF_JX_sph = cart2sphereV(squeeze( SGF_JX(:, :, 1, :) ), TH, PH);
    %% Calculate total electric field
    E = SGF_JX_sph .* exp(-1j * k0 * R_FF) ./ (2*pi * R_FF);
end
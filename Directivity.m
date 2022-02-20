function [ Dir, Prad ] = Directivity( E_tot_abs, TH, dth, dph, er, R_FF )
%Directivity This function computes the directivity and the total radiated
%power of an antenna
%   The function takes the absolute value of the total electric field in a
%   2D matrix with dimensions corresponding to the theta and phi meshgrid
%   in spherical coordinates, the theta meshgrid, the step sizes of the
%   theta and phi arrays, the relative permittivity of the medium, and an
%   array of the radial distance as inputs.
%   It outputs the directivity in a 2D matrix corresponding to the theta
%   and phi meshgrid, and the total radiated power of the antenna producing
%   the absolute value of the total electric field.
    %% Constants
    e0 = 8.8541878128e-12;      % Permittivity of free space [F/m]
    u0 = 4*pi*1e-7;             % Permeability of free space [H/m]
    %% Calculate wave impedance
    Z = sqrt(u0 / (e0*er));
    %% Calculate radiated power
    Prad = 2 * sum( sum( E_tot_abs .* sin(TH) )) * dth * dph;
    %% Calculate radiation intensity
    U = (E_tot_abs.^2) .* (R_FF.^2) / (2 * Z);
    %% Calculate directivity
    Dir = U / (Prad / (4*pi));
end
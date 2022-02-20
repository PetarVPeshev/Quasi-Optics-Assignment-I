function [ Jx ] = FTCurrent( k0, er, KX, KY, L, W )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   Note: fix documentation of FTCurrent routine
    %% Calculate keq
    keq = k0 * sqrt(er);
    %% Calculate FT
    T = sinc(KY * W / (2*pi));
    L = 2*keq * (cos(KX * L/2) - cos(keq * L/2)) ./...
        ( (keq.^2 - KX.^2) .* sin(keq * L/2) );
    Jx = T .* L;
end
function [c_Op, c0_Op, L_Op] = sphericalSplineOperators(coords, m)

% calculate Operators for Spherical Spline Interpolation and Laplacian
%
% [c_Op, c0_Op, L_Op] = sphericalSplineOperators(coords, m)
%
% coords:   cartesian coordinates of the electrode positions on the unit sphere
% m:        order of the spherical spline
% c_Op:     operator to calculate the spherical spline coefficients c_i
% c0_Op:    ", for the coefficient c_0
% L_Op:     operator to calculate the laplacian
%
% see Perrin et al., Spherical Splines for Scalp Potential and Current Density Mapping
%
% Copyright (C) 2003 Carsten Allefeld
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

    cosines = coords * coords';
    N = size(coords, 1);
    
    G = sphericalSpline(cosines, m);
    Gi = inv(G);

    T = ones(N, 1);
    I = eye(N);
    
    c0_Op = T' * Gi / (T' * Gi * T);
    c_Op  = Gi * (I - T * c0_Op);
    
    H = - sphericalSpline(cosines, m - 1);
    
    L_Op = H * c_Op;
    
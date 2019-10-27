function g = sphericalSpline(x, m)

% spherical spline interpolation function g_m(x)
%
% g = sphericalSpline(x, m)
%
% x:    argument, cosine of the angle
% m:    order of the spline
%
% for m = 2, maximal numerical error is of the order 1e-7, smaller for higher m
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

    mlock
    persistent pp
    
    if (size(pp, 2) < m) | isempty(pp{m})
        fprintf(2, 'sphericalSpline: generating piecewise polynomial (for m = %d)\n', m)
        xs = cos((0 :0.005: 1) * pi);
        gs = ss(xs, m);
        pp{m} = spline(xs, gs);
    end
    
    g = ppval(pp{m}, x);
    


function g = ss(x, m)

    g = 0;
    for n = 1 : 100
        g = g + (2 * n + 1) / (n * (n + 1)) ^ m * legendreP(n, 0, x);
    end
    g = g / (4 * pi);
    
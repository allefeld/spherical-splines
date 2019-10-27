function P = legendreP(n, m, X)

% interface to the MATLAB legendre function, specifying the m index parameter
%
% P = legendreP(n, m, X)
%
% Evaluates the associated Legendre polynomial P_n^m at X;
% the resulting P has the same dimensions as X.
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

    if n == 0
        P = legendre(0, X);
    else
        P = legendre(n, X);
        P = reshape(P(m + 1, :), size(X));
    end
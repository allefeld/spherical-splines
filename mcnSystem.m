function [mcnCoords, mcnNames] = mcnSystem(opt)

% calculates idealized EEG electrode positions according to the 10-20 MCN system 
%
% [mcnCoords, mcnNames] = mcnSystem(opt)
%
% opt:        '2d' or 'projected' for two-dimensional coordinates (default: three-dimensional)
% mcnCoords:   cartesian coordinates of the electrodes (n x 3-matrix)
% mcnNames:    names of the electrodes (cell array with strings)
%
% The result of the calculation may be retrieved by the output arguments
% or by the global variables mcn_coords and mcn_names.
%
% The 3D-positions are constructed on the surface of the unit sphere according
% to rules inspired by those of the "International 10-20 System". The
% 2D-positions are obtained by the same rules applied within the plane. The projected
% positions are derived from the three-dimensional coordinates by transforming the
% elevation angle (theta) into a planar radius (sphere2plane). The given electrode
% names follow the "Modified Combinatorial Nomenclature".
%
% cf. Lagerlund et al., Determination of 10-20 Electrode Locations... and
% Perrin et al., Spherical Splines for Scalp Potential and Current Density Mapping
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

    global mcn_names mcn_coords mcn_2d
    
    mcn_names = {};
    mcn_coords = [];
    
    if (nargin == 1) & strcmp(opt, '2d'), mcn_2d = 1; else, mcn_2d = 0; end
    
    setPos('NZ',    0,  62.5)    % nasion
    setPos('IZ',  100,  62.5)    % inion
    
    setPos('FPZ',   0,  50)
    setPos('FP1', -10,  50)
    setPos('FP2',  10,  50)
    setRib(3, {'AF7', 'AF3', 'AFZ', 'AF4', 'AF8'}, -1 :1/2: 1)
    setRib(2, {'F9', 'F7', 'F5', 'F3', 'F1', 'FZ', 'F2', 'F4', 'F6', 'F8', 'F10'}, -5/4 :1/4: 5/4)
    setRib(1, {'FT9', 'FT7', 'FC5', 'FC3', 'FC1', 'FCZ', 'FC2', 'FC4', 'FC6', 'FT8', 'FT10'}, -5/4 :1/4: 5/4)
    setPos({'T9', 'T7', 'C5', 'C3', 'C1', 'CZ', 'C2', 'C4', 'C6', 'T8', 'T10'}, 50, -62.5 :12.5: 62.5)
    setRib(-1, {'TP9', 'TP7', 'CP5', 'CP3', 'CP1', 'CPZ', 'CP2', 'CP4', 'CP6', 'TP8', 'TP10'}, -5/4 :1/4: 5/4)
    setRib(-2, {'P9', 'P7', 'P5', 'P3', 'P1', 'PZ', 'P2', 'P4', 'P6', 'P8', 'P10'}, -5/4 :1/4: 5/4)
    setRib(-3, {'PO7', 'PO3', 'POZ', 'PO4', 'PO8'}, -1 :1/2: 1)
    setPos('O2',  90,  50)
    setPos('O1', 110,  50)
    setPos('OZ', 100,  50)
    
    if (nargin == 1) & strcmp(opt, 'projected')
        mcn_coords = sphere2plane(mcn_coords);
    end
    
    if mcn_2d
        mcn_coords = mcn_coords(:, 1 : 2);
    end
    
    if nargout > 0
        mcnNames = mcn_names;
        mcnCoords = mcn_coords;
    end

    
    
function setPos(name, angle, radius)

% set electrode positions by direct definition

    if ~iscell(name)
        name = {name};
    end
    for i = 1 : size(radius, 2)
        p = position(angle, radius(i));
        register(name{i}, p)
    end

    
    
function setRib(ribval, name, ratio)

% set electrode positions by interpolation along a cirle

    % sagittal electrode position (Z), x0 = 0
    p0 = position(0, ribval * 12.5);
    y0 = p0(2);
    z0 = p0(3);
    % circumference electrode position (7, corresp. 8)
    p1 = position(50 - ribval * 10, 50);%
    x1 = p1(1);
    y1 = p1(2);
    z1 = p1(3);
    % center of the circle through 7-Z-9 electrodes, x = 0
    y = -((-y0^3 + y1*y0^2 + x1^2*y0 + y1^2*y0 + 2*z0*z1*y0 - y1^3 ...
        - x1^2*y1 - y0*z0^2 - y1*z0^2 + 2*y1*z0*z1 - y0*z1^2 - y1*z1^2) ...
        /(2*(y0^2 + y1^2 + z0^2 + z1^2 - 2*y0*y1 - 2*z0*z1)));
    z = -(((z1 - z0)*(-x1^2 + y0^2 - y1^2 + z0^2 - z1^2) ...
        - (2*y1 - 2*y0)*(y1*z0 - y0*z1)) ...
        /((z1 - z0)*(2*z1 - 2*z0) - (y0 - y1)*(2*y1 - 2*y0)));
    % orthogonal vectors in the circle plane and angle between electrodes
    d0 = [0, y0 - y, z0 - z];
    d1 = [x1, y1 - y, z1 - z];
    sp = (d0 * d1') / (norm(d0) * norm(d1));
    alpha = acos(sp);
    ds = (d1 - d0 * (sp)) / sin(alpha);
    % linear interpolation along the circle
    for i = 1 : size(ratio, 2)
        p = [0, y, z] + d0 * cos(alpha * ratio(i)) + ds * sin(alpha * ratio(i));
        register(name{i}, p);
    end


    
function p = position(angle, radius)

% calculate cartesian coordinates corresponding to the position specification

    global mcn_2d
    
    if mcn_2d
        p = radius / 50 * [sin(angle * pi / 100), cos(angle * pi / 100), 0];
    else
        r = sin(radius * pi / 100);
        p = [r * sin(angle * pi / 100), r * cos(angle * pi / 100), cos(radius * pi / 100)];
    end

    
    
function register(name, p)

% register electrode name and coordinates in the global variables

    global mcn_names mcn_coords
    
    mcn_names = [mcn_names; {name}];
    mcn_coords = [mcn_coords; p];

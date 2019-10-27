# Spherical Splines for Scalp Potential and Current Density Mapping

This repository contains code implementing the method published by
Perrin, Pernier, Bertrand, Echallier, Spherical splines for scalp potential and current density mapping, *Electroencephalography and Clinical Neurophysiology* 72(2), 184-187, 1989.

	[c_Op, c0_Op, L_Op] = sphericalSplineOperators(coords, m)

The important output is `L_Op`, which is a square matrix that transforms the
raw data matrix of EEG channels into the matrix of SCD channels. The input
`coords` gives the 3d-coordinates of the electrodes. `m` is the order of the
spherical spline, I normally use `m = 4`.

	g = sphericalSpline(x, m)

This is a helper function for `sphericalSplineOperators` that implements the
spherical spline function.

	[mcnCoords, mcnNames] = mcnSystem(opt)

If `opt` is omitted, this function calculates theoretical 3d-coordinates for
the electrodes of the Modified Combinatorial Nomenclature; the corresponding
names are given in the second output. This can be used to provide coordinates
for `sphericalSplineOperators`.


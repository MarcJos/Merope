Known and expected bugs

# Non-intended use
- [ ] Mérope interface is implemented in Python3, not in Python 2. Hence, always execute python scripts with command `python3`.
- [ ] Intersecting `MicroInclusions` should have the same phase. Otherwise, the result is not predictible. But does not throw an error.

# Approximations
- [ ] Composite Voxels are :
    - approximated on curved surfaces in 3D (lines in 2D),
    - false on edges and vertices in 3D (vertices in 2D).

# Errors related to computations on double
- [ ] Additional matrix phase may occur for polycrystals for rule `Voxel::Center`
- [ ] The criterion for intersection between shapes is sensitive to the floating error. Moving from gcc 9 to gcc 12 breaks most non-regression tests, presumably for that reason. Whether this is due to the compiler alone or the code is not clear yet.

# Informatic bugs
- :warning: Intrications between Python, C++ and OpenMP are dangerous, *cf* [#issue9](https://www-git-cad.intra.cea.fr/DEC/pleiades/merope/merope_bibliotheque/-/issues/9). The *GIL* has sometimes to be disabled, which is dangerous as well.
    - [ ] Enabling openMP for FFT functions causes a sever loss of performance, due to this nontrivial interaction between these 3 components.

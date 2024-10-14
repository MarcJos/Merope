

**Remark** : All the Python command correspond 1 to 1 to C++ commands. Consider reading the doxygen documentation for more details.

**Remark** : All the classes and functions can be subscripted `_3D` or `_2D` depending on the dimension of the considered space.

This page is devoted to the abstract framework for voxelation in Mérope. For concrete examples, please have a look at the [gallery](/doc/Gallery.md).


# Principles of the voxelation

The voxelation procedure aims at turning an abstract description of a structure into a cubic (or squared) regular mesh in order to perform a FFT computation with the underlying model phase-coefficient.
There are different ways to perform this transformation, depending on the information that is stored inside each voxel, which one one the following :
- a `Pure` phase,
- an `Iso`-tropic mixture of different phases,
- an `AnIso`-tropic mixture of different phases.
Finally, the aim is to put the voxelation in the TMFFT format or the [AMITEX format](http://www.maisondelasimulation.fr/projects/amitex/user_guide/_build/html/input_files.html).

We explain first the principles, and then, how to appeal to the adequate Python functions.

:warning: Unfortunately, a typing error lies in the code, where the incorrect ''voxellation'' is used instead of the correct ''voxelation''. For retrocompatibility reasons, we maintain the word ``voxellation'' in the code.

# Voxelation principles and functionalities

## Different formats and computation pipeline

We show the various voxelation formats managed by Mérope, and the natural transformations between them.
Hereafter are explained each manipulation.
<img src="/doc/Pictures/Voxellation_scheme.png" alt="drawing" width="1000"/>

## Python and C++ types

Mérope's code is written in C++.
Here, we have chosen to resort to template methods in order to implement the various transforms from one format to another.
This logic is not easily tractable to Python; hence, the interface differs slightly between the internal C++ code and the Python user interface.

## From a Structure/FieldStructure to a voxelation

### From a Structure to a voxelation

To compute the representation in each voxel, one `VoxelRule` is used:
- `Center` : each voxel is of `Pure` type. That is, each voxel contains a single phase, which is identified from the geometry as the phase of the center of the voxel.
- `Average` : each voxel is of `Iso` type. It represents an `Iso`-tropic mixture by storing a list of phase index associated to a volume fraction inside the voxel.
  - The volume fractions are approximated by considering that all geometrical surfaces might be approximated as planes inside the voxel. (This approximation is satisfactory for spheres of radius far larger than the voxel size, but it is false when attempting to voxelate surface irregularities, such as edges of polyhedron.)
  - Renormalization is performed so that it is guaranteed that the volume fractions sum up to 1 (with an error of order $10^{-6}$).
  - Boolean operations are performed in a coherent way w.r.t. the averaging process. Nevertheless, since they are performed on the level of the voxelations, information is lost in the process, yielding further approximations.
    - :warning: More precisely, the user has to be cautious when appealing to `VoxelRule.Average` in combination with complex microstructures. For example, the following is **not recommended** : define a new `Structure` by `Structure_3D(structure_3D_1, structure_3D_2, mask)` with `mask=structure_3D_2`. Indeed, assume that `structure_3D_1` has only phase 0, and that the `structure_3D_2` has only phases 1 and 0. Then the volume fraction of phase 1 inside a voxel will be approximated as the square of its actual value, which is not accurate for all intermediate values.
- `Laminate` : each voxel is `AnIso` type. It represents an `Iso`-tropic mixture by storing a normal associated with each phase besides an `Iso`-tropic representation. Algorithmically, it appeals to the same functions and concepts as `Average` (hence, it suffers from the same drawbacks).

Here, the phase is naturally an integer, and we think of `Pure(int)`, `Iso(int)` and `AnIso(int)` voxelations.

### From a FieldStructure to a voxelation

These voxelations, which are well-suited to represent discrete geometries, where the phase indicates to which solid a point belongs to, is generalized to real-valued functions on the torus (identified as `FieldStructure` in Mérope).
In particular :
- `Pure(real)` voxelations of such functions boils down to a center-valued $P_0$ discretization,
- `Iso(real)` and `AnIso(real)` voxelations correspond to transforming phases stored `Iso(int)` / `AnIso(int)` into a real value (by allocating a coefficient depending on the phase, for example).

## Composite field, from phase to real values

### Apply coefficients

This operation consists in replacing the phase value by a real value.
It is useful when a phase is associated with a special value (representing for example material constants, such as thermal conductivity).

### Apply texture

This operation consists in transforming the phase value $i$ at a given point $x$ into a real coefficient $f(x, i)$.
Here, $f$ is termed the *texture*.
See an example below, where a polycrystal (phase representation) is turnt into a scalar field with sinusoidal dependance in $x$, but with orientation depending on the crystal phase.

<img src="/doc/Pictures/before_texture.png" alt="drawing" width="300"/>
<img src="/doc/Pictures/after_texture.png" alt="drawing" width="300"/>

See [Z_texture.cxx](modules/merope_core/src/Z_texture.cxx).

## Apply homogenization rule

When considering a voxelation of `Iso`-tropic real-valued voxels, an effective coefficient field might be obtained by appealing to classical homogenization rules.
Mérope proposes the following homogenization rules : `Voigt` (arithmetic average), `Reuss` (harmonic average), `Largest`, `Smallest`.
Each is applied on each `Iso(real)` voxel of coefficients associated with volume fractions $(a_i, \theta_i)$ and results in an effective isotropic coefficent $\bar{a}$ computed as follows:
- `Voigt` : $\bar{a} = \langle a \rangle $
- `Reuss` : $\bar{a} = \langle a^{-1} \rangle^{-1} $
- `Largest` : $\bar{a} = \max_i(a_i) $
- `Smallest` : $\bar{a} = \min_i(a_i) $

## Convert to stl format

The internal format of Mérope's voxelation cannot be easily exported.
Hence, the different grids can be exported as combination of vectors and tuples of the C++ standard library in order to be further manipulated by other codes/interfaces.

## Convert to numpy format

A numpy representation of the grid can be obtained for `Pure(int)` and `Pure(double)` voxel formats.

## Print

As a rule, Mérope's voxelation outputs are written in the `.vtk` format of the [vtk library](https://vtk.org/).
These can be displayed by using [paraview](https://www.paraview.org/) or [Salomé](https://www.salome-platform.org/?lang=fr).

### Standard print

Mérope's only offers printing facilities for grids made of `Pure(real)` or `Pure(int)` voxels.
These are printed as real-valued and integer-valued fields, respectively, by means of the function `printVTK`.

### Segmentation and format for AMITEX/TMFFT

AMITEX and TMFFT only accept phase voxelations where :
- each voxel contains an integer $n \in [0, N]$,
- each phase $n \in [0, N]$ is present in at least one voxel (for $N$ sufficiently small).

Hence before using such a solver :
- real-valued fields should be segmented (thus leading to integer-valued fields). Then, the voxelation comes with a list of coefficients associated to phases.
- integer-valued fields should have their phases potentially changed so that they pave a integral segment $\{0, 1, 2, \dots, N\}$.

Such operations are performed by Mérope by means of functions `printVTK_segmented` and `printVTK_removeUnusedPhase`, respectively.

# Python Interface

All the methods under consideration are stored inside the package `merope.vox`.

## Grid size

Consider a grid of lengths $(L_1, L_2, L_3)$ with $L_i = {\rm{nbNodes}}_i * dx_i$.
Mérope implements a subgrid of voxels of centers $(n_1 * dx_1, n_2*dx_2, n_3 * dx_3)$ for ${\rm nMin}_i \leq n_i < {\rm nMax}_i$.

- `merope.vox.GridParameters_3D` : Simple class for defining a grid.
- `merope.vox.create_grid_parameters_N_L_3D(nbNodes, nMin, nMax, L)` : Return a subgrid `merope.vox.GridParameters_3D` of given dimensions.
- `merope.vox.create_grid_parameters_N_L_3D(nbNodes, L) = merope.vox.create_grid_parameters_N_L_3D(nbNodes, nMin = [0, 0, 0], nMax = nbNodes, L)` : Return a subgrid `merope.vox.GridParameters_3D` of given dimensions. Notice the default parameters for the whole grid.


## Voxel rule

- `merope.vox.VoxelRule.Center` : for building `Pure` voxels
- `merope.vox.VoxelRule.Average` : for building `Iso`-tropic mixtures of different phases
- `merope.vox.VoxelRule.Laminate` : for building `AnIso`-tropic mixtures of different phases

## Handle for grid of composite voxels

A universal handle for grids of composite voxels is given.
It performs inplace all the transfroms requested by the user.

- `merope.vox.GridRepresentation_3D` : Universal hand for a grid of composite voxels.
  - `merope.vox.GridRepresentation_3D(structure_3D, gridParameters_3D, voxelRule)` : Constructor. Builds a composite voxel grid. Depends on:
    - parameter `structure_3D` : of type `merope.Structure_3D`, structure to be voxellized
    - parameter `gridParameters_3D` : of type `merope.vox.GridParameters_3D`, grid dimensions
    - parameter `voxelRule` : chosen among `merope.vox.VoxelRule.Center`, `merope.vox.VoxelRule.Average`, `merope.vox.VoxelRule.Laminate`
  - `merope.vox.GridRepresentation_3D(fieldStructure_3D, gridParameters_3D, voxelRule)` : Constructor. Builds a composite voxel grid. Depends on:
    - parameter `fieldStructure_3D` : of type `merope.FieldStructure_3D`, structure to be voxellized
    - parameter `gridParameters_3D` : of type `merope.vox.GridParameters_3D`, grid dimensions
    - parameter `voxelRule` : chosen among `merope.vox.VoxelRule.Center`, `merope.vox.VoxelRule.Average`, `merope.vox.VoxelRule.Laminate`
  - `gridRepresentation_3D.apply_coefficients(coefficients)` : replace the internal field by replacing phase `n` by given `coefficients[n]`.
    - require : internal representation should be of type `Pure(int)`, `Iso(int)`, `AnIso(int)`
  - `gridRepresentation_3D.apply_homogRule(homogRule)` : apply on the internal field the chosen homogenization rule
    - parameter `homogRule` : chosen among `merope.HomogenizationRule.Voigt` (arithmetic average), `merope.HomogenizationRule.Reuss` (harmonic average), `merope.HomogenizationRule.Largest`, `merope.HomogenizationRule.Smallest`
    - require : internal representation should be of type `Iso(real)`
  - `gridRepresentation_3D.apply_homogRule(homogRule, coefficients)` : equivalent to successively appealing to
    - `gridRepresentation_3D.apply_coefficients(coefficients)`
    - `gridRepresentation_3D.apply_homogRule(homogRule)`
  - `gridRepresentation_3D.convert_to_stl_format()` : converts the internal field into the associated stl format
  - `gridRepresentation_3D.removeUnusedPhase()` : change the phase id so that all phases are inside [0, N] with each phase present at least in one voxel while preserving the phase order.
  - `gridRepresentation_3D.get_PurePhaseField()`: if internal representation is of type `Pure(int)`, return it as a Cartesian grid.
  - `gridRepresentation_3D.get_IsoPhaseField()`: if internal representation is of type `Iso(int)`, return it as a Cartesian grid.
  - `gridRepresentation_3D.get_as_list()` : return the internal representation of the grid as a Python list.
  - `gridRepresentation_3D.__str__()` : string describing the internal state of the handle

## Analyzer

- `merope.vox.GridAnalyzer_3D` : object for analyzing voxelations. No internal state.
  - `merope.vox.GridAnalyzer_3D()` : default constructor
  - `gridAnalyzer.compute_percentages(gridRepresentation)` : computes the volume fraction of each phase in the gridRepresentation
    - result is a map `{phase: volume_fraction}`
    - parameter `gridRepresentation` : of type `merope.vox.GridRepresentation_3D`
    - require : internal state of `gridRepresentation` of type `Pure(int)`, `Iso(int)`, `AnIso(int)`
  - `gridAnalyzer.print_percentages(gridRepresentation)` : print the result of `gridAnalyzer.compute_percentages(gridRepresentation)`


## Printer

- `merope.vox.vtk_printer_3D` : object for printing voxelations. No internal state.
  - `merope.vox.vtk_printer_3D()` : default constructor
  - `printer_3D.printVTK(gridRepresentation, fileVTK, nameValue = "MaterialId")` :
    - parameter `gridRepresentation` : of type `merope.vox.GridRepresentation_3D`, grid to be printed
      - :warning: the internal state of the gridRepresentation should be either Pure(int) or Pure(real)
    - parameter `fileVTK` : name of the .vtk file to be printed
    - parameter `nameValue` : name of the field in the .vtk representation. Default is compatible with AMITEX/TMFFT.
  - `printer_3D.printVTK_segmented(gridRepresentation, fileVTK, fileCoeff, nameValue = "MaterialId")` :
    - parameter `gridRepresentation` : of type `merope.vox.GridRepresentation_3D`, grid to be printed
      - :warning: the internal state of the gridRepresentation should be Pure(real)
    - parameter `fileVTK` : name of the .vtk file to be printed
    - parameter `fileCoeff` : name of the .txt file containing segmented values
    - parameter `nameValue` : name of the field in the .vtk representation. Default is compatible with AMITEX/TMFFT.
  - `printer_3D.printVTK_removeUnusedPhase(gridRepresentation, fileVTK, fileCoeff, nameValue = "MaterialId")` :
    - parameter `gridRepresentation` : of type `merope.vox.GridRepresentation_3D`, grid to be printed
      - :warning: the internal state of the gridRepresentation should be Pure(int)
    - parameter `fileVTK` : name of the .vtk file to be printed
    - parameter `fileCoeff` : name of the .txt file containing phase indices $n_i$ such that the phases $i$ of the .vtk output correspond to the phase $n_i$ of the `gridRepresentation`.
    - parameter `nameValue` : name of the field in the .vtk representation. Default is compatible with AMITEX/TMFFT.

## Export as numpy array

- `merope.vox.NumpyConverter_3D` : object for exporting grids as numpy arrays. No internal state.
  - `merope.vox.NumpyConverter_3D()`: default constructor
  - `numpyConverter_3D.compute_RealField(gridRepresentation)` : return the numpy array(real) from a gridRepresentation in a suitable internal state
  -  `numpyConverter_3D.compute_PhaseField(gridRepresentation)` return the numpy array(int) from a gridRepresentation in a suitable internal state

# Deprecated Python interface

:warning: DEPRECATED. This will not be further developped and is only maintained for non-regression.

- Constructors and setters
  - `merope.Voxellation_3D(multiInclusions_3D)` : constructor from a `MultiInclusions_3D`.
  - `merope.Voxellation_3D(structure_3D)` : constructor from a `Structure_3D`.
  - `merope.Voxellation_3D(fieldFtructure_3D)` : constructor from a `FieldStructure_3D`.
  - `setPureCoeffs(coeffsList)` : associate to each phase `i` a coefficient `coeffList[i]`
    - :warning: For nonpositive coefficient values, does not work
    - :warning: The `coeffsList` is not verified to contain the `coeffList[i]`. It is highly recommanded to make use of method `getAllPhases` to verify it a priori
  - `setVoxelRule(voxelRule)` : define how to represent the content of each voxel. `voxelRule` is an object `VoxelRule.X` for `X` being `Average` or `Center`.
  - `setHomogRule(homogRule)` : define the rule to compute the coefficient of each voxel, from the percentage of each phase in it. Parameter `homogRule` is a object `HomogenizationRule.X` for `X` being `Reuss`, `Voigt`, `Smallest` or `Largest`.
- Functions
  - `proceed([N1, N2, N3], NMin=[0,0,0], NMax = [N1, N2, N3])` : build a voxelation of width `[N1, N2, N3]` voxels in the 3 directions. If the  optional parameters NMin and NMax are activated (both or none should be activated), extracts from the grid the voxels of indices `NMin[i]<=index[i]<NMax[i]`. It does it in an efficient way, avoiding to compute the whole grid.
- Ouptuts  
  - `getField_Numpy()` : return, if applicable, the output in the form of a numpy.array containing a scalar field.
  - `getField()` : return, if applicable, the output as a Mérope object `GridField`, which contains the output in the form of a scalar field.
  - `computePhaseGrid([N1, N2, N3])` :  build a voxelation of width `[N1, N2, N3]` voxels in the 3 directions and returns it in the form of a list describing in each voxel the percentage of each phase.
  - `printFile(Zone.vtk, Coeffs.txt)` : write 2 output files. The first one contains the *.vtk* file containing the zones of the material. The second one contains the associated coefficients (important for segmented image) in the ASCII format. These two outputs are expected by the `FFT` solvers.
  - `printFieldFile(File.vtk)` : write a single output file, containing the field (with real values).


# Effects of discretization and composite voxels

In some situations (in particular, for high contrasts), composite voxels can save a lot of computational resources, see \[Schneider, 2021\].
We **highly recommend** to the user to proceed with a small (empirical) study of the discretization error before making use of large RVEs.
[Here is an example.](/studies/Voxels_composites)

Here is an example in which we consider lead spheres coated with gold inside water, and we search for the effective conductivity. The contrast is quite high ; indeed, the three materials have the following conductivities :
$\lambda_{\rm gold} = 317 Wm^{-1}K^{-1}$, $\lambda_{\rm lead} = 35.3Wm^{-1}K^{-1}$, $\lambda_{\rm water} = 0.606Wm^{-1}K^{-1}$.

We study the convergence of the effective conductivity computed by the RVE method, by only varying the discretization parameters.
The RVE under concern is fixed and quite small (with only 19 spheres).
We make use of the 5 different methods, namely `merope.vox.VoxelRule.Center`, and `merope.vox.VoxelRule.Average` with the following `homogenizationRule` : `merope.HomogenizationRule.Largest`,  `merope.HomogenizationRule.Smallest`, `merope.HomogenizationRule.Voigt`, `merope.HomogenizationRule.Reuss`.
When letting the discretization step go to 0, we observe :
- there is a high discrepancy in the results between the 5 methods, but no significant difference in terms of computation time,
- the `Reuss` method converges faster than the other ones. For 256 x 256 x 256 voxels (which is not small, since there are 1.6 1e7 voxels), the `Reuss` methods has a relative error of order 1%, the `Smallest` method has a relative error of order 10%, and the 3 other methods suffer from a relative error larger than 100%.

<img src="/studies/Voxels_composites/Relative_error.png" alt="drawing" width="500"/>
<img src="/doc/Pictures/Plomb_or.png" alt="drawing" width="400"/>



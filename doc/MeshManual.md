
# Description

Mérope translates its internal representation of microstructures into a `.geo` script file that can be read by `gmsh`.
Thus, this is only **mesh parametrization**.

## Principle

Unlike voxellation, the translation from a periodic RVE to a mesh is not straightforward ; indeed, special care should be dedicated to *periodic* boundaries of the mesh.

The question of *boundaries* in a periodic geometry is unnatural from a fundamental point of view : indeed, periodic cuboid does *not* have any boundary (only its projection in the Euclidean space, namely the usual cuboid, has boundaries). However, for practical concern (since, apart from CGAL, mesh tools proceed in the usual Euclidean space), we shall define boundaries. Since we have the freedom to choose them, we decide to use *periodic boundaries* such that no inclusion intersects them. Once this periodic boundary is built, the whole mesh is easily obtained by filling in the inner volume with the inclusions. This procedure comes with an obvious advantage : since no inclusion intersects the boundary, there is no artificial mesh refinement due to inclusions intersecting a boundary. Yet, this comes at a cost : the chosen boundaries are not anymore *simple* cuboids but rather polygonal surfaces.

For Laguerre tessellations, *periodic boundaries* can obviously be chosen as the boundary of the union of the projection of each tessel in the Euclidean space (this is a choice made for example by `Neper`).
For the case of objects $O_i$ that are separated from each other by a minimal distance $d$, the choice is not as obvious.
Yet, we may chose a strategy to get back to the previous case. We cover each object $O_i$ by a finite number of spheres $S_{i,j}$, sufficiently small so that spheres related to two different objects are separated by a minimal distance $d/2$, build the Laguerre tessellations associated to the spheres, and merge all the tessels related to the spheres $S_{i,j}$ associated to each object $O_i$, thus obtaining a polygonal volume $V_i$. Then, the object $O_i$ is contained into the volume $V_i$, and is at a distance of at least $d/4$ of its boundary. As a consequence, we may take the boundary of the union of the projections of the volumes $V_i$ in the Euclidean space as the *periodic boundary*, which enjoys the desirable property that no inclusion intersects this boundary.

## Functionalities 

The class `Merope.mesh.MeshGenerator` is used to produce a `.geo` script from the Mérope class `Merope.MultiInclusions_3D`.
The following features are supported :
- microstructures :
    - Laguerre tessellations,    
    - sphere inclusions, **but** there should be a minimal (nonzero) distance between the spheres,
    - Cylinder inclusions, **but** the cylinders should be included in spheres separated by a (nonzero) distance,
- inner layer definitions,
- holes (ignore some inclusions/layers when meshing),

## Limits

- This strategy is limited to well-separated inclusions or to tessellations.
- For technical reasons (at least), boolean operations cannot be implemented, due to the `gmsh` way of performing boolean operations, which interacts in a destructive way with the implementation of periodic boundary conditions.

# Manual

## `Merope.mesh.MeshGenerator`

`Merope.mesh.MeshGenerator` supports the following functions :
- `setMultiInclusions(multiInclusions_3D)` : acquire the microstructure to be meshed (:warning: the microstructure should satisfy the above implicit assumptions).
- `write(nameGeoFile)` : write the mesh commands for `gmsh`  in `.geo` format into the file `nameGeoFile`.
- `setMeshOrder(meshOrder)` / `getMeshOrder()` : set/get the mesh order (=order of the elements) to be equal to the `meshOrder` (default 2).
- `setMeshSize(meshSize)`/`getMeshSize()` : set the minimal mesh size to be equal to `meshSize` (default = 0.05).
- `setBinaryOutput(boolean)` : trigger whether the mesh file should be written in binary format or not (default `False`).
- `do_not_mesh(list_of_phases)` : the solid part corresponding to the phases will not be meshed.
- `set_nameOutput([name_output_1, ...])` : set which outputs gmsh should produce. Notice gmsh automatically recognize the extension ('.vtk', '.msh') from the name of the output.

The additional functions can be employed, but should not unless the programmer understands how they work :
- `setAdimMergeDistance0(epsilon_0)` : set the distance criterion for merging elements of the original Laguerre tessellation, before identifiying periodic faces, when building the periodic enveloppe (it is adimensional, that is, it will be multiplied by the minimal edge length of the periodic cuboid). (Default = 1e-5.)
    - :warning: This is a technical parameter, that should be small.
    - :warning: In any case, it is assumed that $\epsilon_0 \leq \epsilon_1$.
- `setAdimMergeDistance1(epsilon_1)` : set the distance criterion for merging elements of the original Laguerre tessellation, after identifiying periodic faces, when building the periodic enveloppe (it is adimensional, that is, it will be multiplied by the minimal edge length of the periodic cuboid).  (Default = 1e-5.)
    - :warning: This is a technical parameter, that should be small.
    - :warning: In any case, it is assumed that $\epsilon_0 \leq \epsilon_1$.

# Examples

## Spherical inclusions

<img src="/doc/Pictures/Mesh_200spheres.png" alt="drawing" width="500"/>
See [mesh_spheres_0.py](/tests/microstructures/mesh_spheres_0/mesh_spheres_0.py)

## Laguerre tessellations

<img src="/doc/Pictures/Mesh_Polyhedron.png" alt="drawing" width="500"/>
See [mesh_poly_0.py](/tests/microstructures/mesh_poly_0/mesh_poly_0.py)

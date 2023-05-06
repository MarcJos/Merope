[[_TOC_]]

# Simplified Interface for building geometric microstructures (obsolete)

This Interface allows for a less flexible but more straightforward way of building microstructures.
It is highly inspired from the previous software `VER`.

## For MultiInclusions

A more modular interface for `MultiInclusions_3D` is `SimpleMultiInclusions_3D`.
:construction:

**Methods**
- `setSpheres(sphereList)` : set the spheres defining the structure
- `getSpheres()` : return a sphereList
- `fromFile(fileName)` : set the sphere and length from a given file (format PaQhull)
- `setAspRatio([a1, a2, a3])` : set the aspect ratio
- `setTypeCrystal(typeCrystal)` : set the type of inclusions, among the following : `TypeCrystal.Voronoi`, `TypeCrystal.Laguerre`, `TypeCrystal.JohnsonMehl` (inactive), `TypeCrystal.Spheres`
- `setLayerList(layerList)` : set the layer list for each inclusion (not recommended)
- `getLayerList()` : return the list of layers

## For Structure

A non-recursive simplified interface `SimpleStructure_3D` is provided to build `Structures_3D`, in order to be closer to the ancient versions of VER of Mérope.
It is a structure limited to, at most, 2 `MultiInclusions_3D`.

**Main methods** :
- `SimpleStructure_3D()` : constructor
- `setLength([L1, L2, L3])` (*optional*, default = [1, 1, 1]) : set the lengths of the box
- `setColorization(colorMaterialID)` (*optional*, default = ColorMaterialID.Poly) : define what type of predefined template structure you wish.
  - `ColorMaterialID.Poly` : only the mainInclusions are considered, and all inclusions have different phases.
  - `ColorMaterialID.Erode` : only the mainInclusions are considered, all inclusions have a single layer. All the core inclusions have different phases, all the layers have the same phase (=`erosionPhase`)
  - `ColorMaterialID.Erode2mat` : only the mainInclusions are considered, all inclusions have a single layer. All the core inclusions have a the same phase (=`innerPhase`), all the layers have the same phase (=`erosionPhase`).
  - `ColorMaterialID.Erode3mat` : like `ColorMaterialID.Erode2mat`. Additionnaly, the layer are intersected with spheres, giving rise to a third domain, with phase `erosionInclusionsPhase`.
- `build()` : return the `Structure_3D` that has been parametrized.
- `setErosionWidth(width)` (*optional*, default = 0) : set the erosion width of the `mainInclusions`.
- `setInnerPhase(phaseNumber)` (*optional*, default=0) : set the phase number of the core of the inclusions.
- `setErosionPhase(phaseNumber)` (*optional*, default=1) : set the phase number of the eroded phase.
- `setErosionInclusionsPhase(phaseNumber)` (*optional*, default=2) : set the phase number of the spheres inside the inclusions phase.
- `frac2erosionWidth(frac)` : infer the erosion width (heuristically) such that the volume fraction of the layers is close to `frac`.
- `erosionWidth2frac()` : infer the volume fraction of the layers.

The `SimpleStructure_3D` has two `Pre_InterfaceMultiInclusions_3D` properties : 
- `mainInclusions` (used in all the cases),
- `secdInclusions` (only use for the case `ColorMaterialID.Erode3mat`, in order to define the porosity spheres).
These are of type *close* to `SimpleMultiInclusions_3D` (but do not have the `layerList` methods).

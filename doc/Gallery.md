
# Gallery

## Meshed structures

### Spherical inclusions (meshed)

<img src="/doc/Pictures/Mesh_200spheres.png" alt="drawing" width="500"/>
See [mesh_spheres_0.py](https://github.com/MarcJos/Merope/tests/microstructures/mesh_spheres_0/mesh_spheres_0.py)

### Spherical inclusions with holes/empty bubbles (meshed)

<img src="/doc/Pictures/mesh_bubbles.png" alt="drawing" width="500"/>
See [mesh_spheres_1.py](https://github.com/MarcJos/Merope/tests/microstructures/mesh_spheres_1/mesh_spheres_1.py)

### Polyhedral inclusions with holes (meshed)

<img src="/doc/Pictures/mesh_holes.png" alt="drawing" width="500"/>
See [mesh_poly_1.py](https://github.com/MarcJos/Merope/tests/microstructures/mesh_poly_1/mesh_poly_1.py)

### Cylinder inclusions with holes (meshed)
<img src="/doc/Pictures/mesh_cylinders.png" alt="drawing" width="500"/>
See [mesh_cylinder.py](https://github.com/MarcJos/Merope/tests/microstructures/mesh_cylinder/mesh_cylinder.py)


## Inclusions

### Overlapping coated spherical inclusions

<img src="/doc/Pictures/overlapping_coated_inclusions.png" alt="drawing" width="500"/>
See [overlapping_coated_inclusions.py](https://github.com/MarcJos/Merope/tests/microstructures/coated_inclusions/coated_inclusions.py)

### Inclusions of 2 types, with composite voxels.

<img src="/doc/Pictures/Voxels_composites.png" alt="drawing" width="500"/>
See [Inclusions_Voxellisation.py](https://github.com/MarcJos/Merope/tests/microstructures/inclusions/Inclusions_Voxellisation.py).

### Inclusions with layers

<img src="/doc/Pictures/VER_1.png" alt="drawing" width="500"/>
See [Inclusions_Voxellisation.py](https://github.com/MarcJos/Merope/tests/microstructures/multiLayer/Inclusions_Voxellisation.py)

### Cylinder inclusions

<img src="/doc/Pictures/cylinders.png" alt="drawing" width="500"/>
See [many_cylinders_random_orient.py](https://github.com/MarcJos/Merope/tests/microstructures/many_cylinders_random_orient/many_cylinders_random_orient.py)

## PolyCrystals

### Regular hexagonal crystals, with a layer of another phase.

<img src="/doc/Pictures/hexagones.png" alt="drawing" width="500"/>
See [hexagones.py](https://github.com/MarcJos/Merope/tests/microstructures/hexagones/hexagones.py)

### Simple PolyCrystals (Laguerre tesselation)

<img src="/doc/Pictures/Laguerre.png" alt="drawing" width="500"/>
See [polyCrystal.py](https://github.com/MarcJos/Merope/tests/microstructures/polyCrystal/polyCrystal.py)


### PolyCrystals with prescribed crystallite volumes

<div align="center">
<img src="/doc/Pictures/Original.png" alt="drawing" width="250"/>
<img src="/doc/Pictures/All_equal.png" alt="drawing" width="250"/>
<img src="/doc/Pictures/Center_small.png" alt="drawing" width="250"/>
</div>
See [optimize_Laguerre_tess](https://github.com/MarcJos/Merope/tests/microstructures/optimize_Laguerre_2D/optimize_Laguerre_tess.py).


### PolyCrystals with layers

<img src="/doc/Pictures/PolyCrystal.png" alt="drawing" width="350"/>
See [polyCrystal.py](https://github.com/MarcJos/Merope/tests/microstructures/multiLayer/polyCrystal_Voigt.py)

<img src="/doc/Pictures/Poly_2D.png" alt="drawing" width="350"/>
See [polyCrystal_2D.py](https://github.com/MarcJos/Merope/tests/microstructures/polyCrystal_2D/polyCrystal_2D.py)

## Gaussian fields

<img src="/doc/Pictures/Gauss.png" alt="drawing" width="350"/>
See [parallel_gaussian.py](https://github.com/MarcJos/Merope/tests/microstructures/parallel_gaussian/parallel_gaussian.py)

## Deterministic scalar field

<img src="/doc/Pictures/prescribedField.png" alt="drawing" width="350"/>
See [prescribedField.py](https://github.com/MarcJos/Merope/tests/microstructures/prescribedField/prescribedField.py)

## Complex structures

### PolyCrystals with layers and porosity

<img src="/doc/Pictures/Laguerre_Filamentaire.png" alt="drawing" width="500"/>
See [polyCrystal_filamentaire.py](https://github.com/MarcJos/Merope/tests/microstructures/polyCrystal_filamentaire/polyCrystal_filamentaire.py) (with `SimpleStructure_3D`)

### Use of mask

<img src="/doc/Pictures/Mask/Mask_Hexa.png" alt="drawing" width="300"/> and 
<img src="/doc/Pictures/Mask/Mask_Poly.png" alt="drawing" width="300"/> with mask 
<img src="/doc/Pictures/Mask/Mask_Mask.png" alt="drawing" width="300"/> gives  
<img src="/doc/Pictures/Mask/Mask_Final.png" alt="drawing" width="500"/> 
See [recursive_structure.py](https://github.com/MarcJos/Merope/tests/microstructures/recursive_structure/recursive_structure.py)

## Additional structures

Additional structures may be consulted in [CEA_Microstructures](https://www-git-cad.intra.cea.fr/DEC/pleiades/merope/merope_nucleaire/-/blob/master/doc/Microstructures.md).

# Didactic examples

- Build a **structure** and **voxellation** [buildVoxellation.py](https://github.com/MarcJos/Merope/tests/microstructures/buildVoxellation/buildVoxellation.py)
- Compute **thermal conductivity** with **amitex** [Thermal_amitex.py](https://github.com/MarcJos/Merope/tests/microstructures/buildVoxellation/Thermal_amitex.py)
- Compute **thermal conductivity** with **tmfft** [Thermal_tmfft.py](https://github.com/MarcJos/Merope/tests/microstructures/buildVoxellation/Thermal_tmfft.py)

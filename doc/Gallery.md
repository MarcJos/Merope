
# Gallery

## Meshed structures

### Spherical inclusions (meshed)

See [mesh_spheres_0.py](https://github.com/MarcJos/Merope/tests/microstructures/mesh_spheres_0/mesh_spheres_0.py)  
<img src="/doc/Pictures/Mesh_200spheres.png" alt="drawing" width="500"/>

### Spherical inclusions with holes/empty bubbles (meshed)

See [mesh_spheres_1.py](https://github.com/MarcJos/Merope/tests/microstructures/mesh_spheres_1/mesh_spheres_1.py)  
<img src="/doc/Pictures/mesh_bubbles.png" alt="drawing" width="500"/>

### Polyhedral inclusions with holes (meshed)

See [mesh_poly_1.py](https://github.com/MarcJos/Merope/tests/microstructures/mesh_poly_1/mesh_poly_1.py)  
<img src="/doc/Pictures/mesh_holes.png" alt="drawing" width="500"/>

### Cylinder inclusions with holes (meshed)

See [mesh_cylinder.py](https://github.com/MarcJos/Merope/tests/microstructures/mesh_cylinder/mesh_cylinder.py)  
<img src="/doc/Pictures/mesh_cylinders.png" alt="drawing" width="500"/>


## Inclusions

### Overlapping coated spherical inclusions

See [overlapping_coated_inclusions.py](https://github.com/MarcJos/Merope/tests/microstructures/coated_inclusions/coated_inclusions.py)  
<img src="/doc/Pictures/overlapping_coated_inclusions.png" alt="drawing" width="500"/>

### Inclusions of 2 types, with composite voxels.

See [Inclusions_Voxellisation.py](https://github.com/MarcJos/Merope/tests/microstructures/inclusions/Inclusions_Voxellisation.py).  
<img src="/doc/Pictures/Voxels_composites.png" alt="drawing" width="500"/>

### Inclusions with layers

See [Inclusions_Voxellisation.py](https://github.com/MarcJos/Merope/tests/microstructures/multiLayer/Inclusions_Voxellisation.py)  
<img src="/doc/Pictures/VER_1.png" alt="drawing" width="500"/>

### Cylinder inclusions

See [many_cylinders_random_orient.py](https://github.com/MarcJos/Merope/tests/microstructures/many_cylinders_random_orient/many_cylinders_random_orient.py)  
<img src="/doc/Pictures/cylinders.png" alt="drawing" width="500"/>

## PolyCrystals

### Regular hexagonal crystals, with a layer of another phase.

See [hexagones.py](https://github.com/MarcJos/Merope/tests/microstructures/hexagones/hexagones.py)  
<img src="/doc/Pictures/hexagones.png" alt="drawing" width="500"/>

### Simple PolyCrystals (Laguerre tesselation)

See [polyCrystal.py](https://github.com/MarcJos/Merope/tests/microstructures/polyCrystal/polyCrystal.py)  
<img src="/doc/Pictures/Laguerre.png" alt="drawing" width="500"/>

### PolyCrystals with prescribed crystallite volumes

See [optimize_Laguerre_tess](https://github.com/MarcJos/Merope/tests/microstructures/optimize_Laguerre_2D/optimize_Laguerre_tess.py).  
<div align="center">
<img src="/doc/Pictures/Original.png" alt="drawing" width="250"/>
<img src="/doc/Pictures/All_equal.png" alt="drawing" width="250"/>
<img src="/doc/Pictures/Center_small.png" alt="drawing" width="250"/>
</div>


### PolyCrystals with layers

See [polyCrystal.py](https://github.com/MarcJos/Merope/tests/microstructures/multiLayer/polyCrystal_Voigt.py)  
<img src="/doc/Pictures/PolyCrystal.png" alt="drawing" width="350"/>

See [polyCrystal_2D.py](https://github.com/MarcJos/Merope/tests/microstructures/polyCrystal_2D/polyCrystal_2D.py)  
<img src="/doc/Pictures/Poly_2D.png" alt="drawing" width="350"/>

## Gaussian fields

See [parallel_gaussian.py](https://github.com/MarcJos/Merope/tests/microstructures/parallel_gaussian/parallel_gaussian.py)  
<img src="/doc/Pictures/Gauss.png" alt="drawing" width="350"/>

## Deterministic scalar field

See [prescribedField.py](https://github.com/MarcJos/Merope/tests/microstructures/prescribedField/prescribedField.py)  
<img src="/doc/Pictures/prescribedField.png" alt="drawing" width="350"/>

## Complex structures

### PolyCrystals with layers and porosity

See [polyCrystal_filamentaire.py](https://github.com/MarcJos/Merope/tests/microstructures/polyCrystal_filamentaire/polyCrystal_filamentaire.py) (with `SimpleStructure_3D`)  
<img src="/doc/Pictures/Laguerre_Filamentaire.png" alt="drawing" width="500"/>

### Use of mask

See [recursive_structure.py](https://github.com/MarcJos/Merope/tests/microstructures/recursive_structure/recursive_structure.py)  
<img src="/doc/Pictures/Mask/Mask_Hexa.png" alt="drawing" width="300"/> and 
<img src="/doc/Pictures/Mask/Mask_Poly.png" alt="drawing" width="300"/> with mask 
<img src="/doc/Pictures/Mask/Mask_Mask.png" alt="drawing" width="300"/> gives  
<img src="/doc/Pictures/Mask/Mask_Final.png" alt="drawing" width="500"/> 

## Additional structures

Additional structures may be consulted in [CEA_Microstructures](https://www-git-cad.intra.cea.fr/DEC/pleiades/merope/merope_nucleaire/-/blob/master/doc/Microstructures.md).

# Didactic examples

- Build a **structure** and **voxellation** [buildVoxellation.py](https://github.com/MarcJos/Merope/tests/microstructures/buildVoxellation/buildVoxellation.py)
- Compute **thermal conductivity** with **amitex** [Thermal_amitex.py](https://github.com/MarcJos/Merope/tests/microstructures/buildVoxellation/Thermal_amitex.py)
- Compute **thermal conductivity** with **tmfft** [Thermal_tmfft.py](https://github.com/MarcJos/Merope/tests/microstructures/buildVoxellation/Thermal_tmfft.py)

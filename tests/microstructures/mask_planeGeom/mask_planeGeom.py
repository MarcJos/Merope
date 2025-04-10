import numpy as np
import merope

import sac_de_billes as sdb


sphere_radius = 0.5
L = np.ones(3)*6*sphere_radius
nbNodes = np.ones(3, dtype=int)*64


#### Generate the structure with the two core spheres
core_spheres   = np.array([sdb.Sphere_3D(np.ones(3)*2*sphere_radius, sphere_radius,1),
                           sdb.Sphere_3D(np.ones(3)*2*sphere_radius + np.array([0,1.8*sphere_radius,0]), sphere_radius, 2)])
si = merope.SphereInclusions_3D()
si.setLength(L)
si.setSpheres(core_spheres)
mi = merope.MultiInclusions_3D()
mi.setInclusions(si)
ids = mi.getAllIdentifiers()
mi.changePhase(ids,[1 for _ in range(len(ids))])
mi.setMatrixPhase(0)
core_structure = merope.Structure_3D(mi)


#### Same structure but with slightly bigger spheres, that will be the layer
layer_spheres = np.array([sdb.Sphere_3D(np.ones(3)*2*sphere_radius, sphere_radius + sphere_radius/5,3),
                           sdb.Sphere_3D(np.ones(3)*2*sphere_radius + np.array([0,1.8*sphere_radius,0]), sphere_radius+ sphere_radius/5, 4)])

si = merope.SphereInclusions_3D()
si.setLength(L)
si.setSpheres(layer_spheres)
mi = merope.MultiInclusions_3D()
mi.setInclusions(si)
mi.setMatrixPhase(0)
ids = mi.getAllIdentifiers()
mi.changePhase(ids,[2 for _ in range(len(ids))])
layer_structure = merope.Structure_3D(mi)

#### Structure with uniform background
mi=merope.MultiInclusions_3D()
mi.setMatrixPhase(5)
mask_structure = merope.Structure_3D(mi)

#### Final structure of spheres with additional layer
structure = merope.Structure_3D(layer_structure, mask_structure, core_structure)

#### Mappings for PolyGeom
all_phases = structure.getAllPhases()
con_phase = 9
print(con_phase)

mappings = {}
for i in range(30): 
    for j in range(30):
        mappings[(i,j)] = con_phase

#### Make the grid
def make_grid(structure, typeOfRule, nameFile):
    gridParameters = merope.vox.create_grid_parameters_N_L_3D(nbNodes, L)
    if typeOfRule == merope.vox.VoxelRule.PolyGeom:
        gridRepr = merope.vox.GridRepresentation_3D(structure, 
                                                gridParameters, 
                                                merope.vox.VoxelRule.PolyGeom,
                                                mappings)
        gridRepr.convert_to_Iso_format()
    else:
        gridRepr = merope.vox.GridRepresentation_3D(structure, 
                                                gridParameters, 
                                                typeOfRule)
    #### Apply real conductivities values and homogenize
    gridRepr.apply_coefficients([i for i in range(30)])
    gridRepr.apply_homogRule(merope.HomogenizationRule.Voigt)
    #### Print the grid
    my_printer = merope.vox.vtk_printer_3D()
    my_printer.printVTK(gridRepr, nameFile)

make_grid(structure, merope.vox.VoxelRule.PolyGeom, "planeGeom.vtk")
make_grid(structure, merope.vox.VoxelRule.Average, "average.vtk")

make_grid(layer_structure, merope.vox.VoxelRule.PolyGeom, "layer.vtk")
make_grid(mask_structure, merope.vox.VoxelRule.PolyGeom, "mask.vtk")
make_grid(core_structure, merope.vox.VoxelRule.PolyGeom, "core.vtk")


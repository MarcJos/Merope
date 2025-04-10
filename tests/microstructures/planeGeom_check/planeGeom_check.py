import sac_de_billes
import merope
import numpy as np


L = [5, 5, 5]
nbNodes=[64 for i in range(3)]

def compute_volume(structure, voxelRule, mapping = None):
    gridParameters = merope.vox.create_grid_parameters_N_L_3D(nbNodes, L)
    if mapping == None:
        gridRepr = merope.vox.GridRepresentation_3D(structure, 
                                                gridParameters, 
                                                voxelRule)
    else:
        gridRepr = merope.vox.GridRepresentation_3D(structure, 
                                                gridParameters, 
                                                mapping)
    gridRepr.convert_to_Iso_format()
    #### Apply real conductivities values and homogenize
    gridRepr.apply_coefficients([i for i in range(30)])
    gridRepr.apply_homogRule(merope.HomogenizationRule.Voigt)
    #### Print the grid

    volume = 0
    liste_phaseFrac = gridRepr.get_as_list()
    for ph in liste_phaseFrac:
        volume += ph
    volume /= nbNodes[0] * nbNodes[1] * nbNodes[2]
    volume *= L[0] * L[1] * L[2]
    return volume

def display_volume(structure, voxelRule, exp_volume, tolerance, mapping=None):
    comp_volume = compute_volume(structure, voxelRule, mapping)
    relative_error = np.abs(comp_volume - exp_volume) / exp_volume
    print("comp_volume : ", comp_volume, " exp_volume : ", exp_volume, " relative_error : ", relative_error)
    if (relative_error > tolerance):
        raise Exception('too large error')

# case 1 : single sphere with layer
def make_XP_1():
    delta_layer = 0.001
    radius = 1

    sphere1 = sac_de_billes.Sphere_3D([0, 0, 0], radius, 0)
    sphInc = merope.SphereInclusions_3D()
    sphInc.setLength(L)
    sphInc.setSpheres([sphere1])

    mi = merope.MultiInclusions_3D()
    mi.setInclusions(sphInc)
    mi.addLayer(mi.getAllIdentifiers(), [1 for _ in mi.getAllIdentifiers()], [delta_layer for _ in mi.getAllIdentifiers()])
    mi.setMatrixPhase(0)

    structure = merope.Structure_3D(mi)
    
    exp_volume = 4./3. * np.pi * (radius**3 - (radius-delta_layer)**3)
    print("Test 1 : planeGeom")
    display_volume(structure, merope.vox.VoxelRule.PolyGeom, exp_volume, 1e-2)
    print("Test 1 : Iso")
    display_volume(structure, merope.vox.VoxelRule.Average, exp_volume, 1e-2)

# case 2 : 2 spheres with layer
def make_XP_2():
    delta_layer = 0.001
    radius = np.pi/3

    sphere1 = sac_de_billes.Sphere_3D([0, 0, 0], radius, 0)
    sphere2 = sac_de_billes.Sphere_3D([2 * radius + 1e-6, 0, 0], radius, 0)
    sphInc = merope.SphereInclusions_3D()
    sphInc.setLength(L)
    sphInc.setSpheres([sphere1, sphere2])

    mi = merope.MultiInclusions_3D()
    mi.setInclusions(sphInc)
    mi.addLayer(mi.getAllIdentifiers(), [1 for _ in mi.getAllIdentifiers()], [delta_layer for _ in mi.getAllIdentifiers()])
    mi.setMatrixPhase(0)

    structure = merope.Structure_3D(mi)
    
    exp_volume = 2 * 4./3. * np.pi * (radius**3 - (radius-delta_layer)**3)
    print("Test 2 : planeGeom")
    display_volume(structure, merope.vox.VoxelRule.PolyGeom, exp_volume, 1e-2)
    print("Test 2 : Iso")
    display_volume(structure, merope.vox.VoxelRule.Average, exp_volume, 1e-2)

# case 3 : sphere with layer by boolean intersections
def make_XP_3():
    delta_layer = 0.001
    radius = np.pi/3

    sphere1 = sac_de_billes.Sphere_3D([0, 0, 0], radius, 1)
    sphere2 = sac_de_billes.Sphere_3D([0, 0, 0], radius - delta_layer, 0)
    sphInc1 = merope.SphereInclusions_3D()
    sphInc1.setLength(L)
    sphInc1.setSpheres([sphere1])
    sphInc2 = merope.SphereInclusions_3D()
    sphInc2.setLength(L)
    sphInc2.setSpheres([sphere2])

    mi1 = merope.MultiInclusions_3D()
    mi1.setInclusions(sphInc1)
    mi1.setMatrixPhase(0)

    mi2 = merope.MultiInclusions_3D()
    mi2.setInclusions(sphInc2)
    mi2.setMatrixPhase(1)

    mi3 = merope.MultiInclusions_3D()
    mi3.setMatrixPhase(0)
    
    s1 = merope.Structure_3D(mi1)
    s2 = merope.Structure_3D(mi2)
    s3 = merope.Structure_3D(mi3)
    structure = merope.Structure_3D(s3, s2, s1)
    
    exp_volume = 4./3. * np.pi * (radius**3 - (radius-delta_layer)**3)
    print("Test 3 : planeGeom")
    display_volume(structure, merope.vox.VoxelRule.PolyGeom, exp_volume, 1e-2)
    print("Test 3 : Iso")
    display_volume(structure, merope.vox.VoxelRule.Average, exp_volume, 1000)


# case 4 : polyhedron
def make_XP_4():
    delta_layer = 0.5 * L[0] / nbNodes[0]
    radius = np.pi/3
    liste_vertices = [[delta_layer + 0, delta_layer + 0, delta_layer + 0], [delta_layer + 1, delta_layer + 0, delta_layer + 0], [delta_layer + 0, delta_layer + 1, delta_layer + 0], [delta_layer + 0, delta_layer + 0, delta_layer + 1]]
    face_indices = [[0, 1, 3], [0, 2, 1], [0, 3, 2], [1, 2, 3]]
    polyFactory = merope.microInclusion.PolyhedronFactory_3D()
    my_poly = polyFactory.fromVertices(1, liste_vertices, face_indices)

    pinc = merope.PolyInclusions_3D()
    pinc.setLength(L)
    pinc.setInclusions([my_poly])

    minc = merope.MultiInclusions_3D()
    minc.setInclusions(pinc)
    minc.setMatrixPhase(0)
    structure = merope.Structure_3D(minc)
    exp_volume = 0.5 / 3
    print("Test 4 : planeGeom")
    display_volume(structure, merope.vox.VoxelRule.PolyGeom, exp_volume, 1e-10)
    print("Test 4 : Iso")
    display_volume(structure, merope.vox.VoxelRule.Average, exp_volume, 1)



make_XP_1()
make_XP_2()
make_XP_3()
make_XP_4()
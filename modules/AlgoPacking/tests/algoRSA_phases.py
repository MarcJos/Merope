### import the library (should accessible through PYTHON_PATH)
import sac_de_billes as sac_de_billes

### Select the big shape in which spheres are thrown
### Among : Tore, Cube (=Cuboid), Sphere, Cylinder
shape = sac_de_billes.NameShape.Tore
### Select the dimensions of the cuboid (sharply) containing the shape
### At the moment, indirectly defines the size of the shape for 
##### Sphere via : L = [2 * radius, 2 * radius, 2 * radius]
##### Cylinder via : L = [2 * radius, 2 * radius, height]
L = [16., 16., 16.]
### Select the desired {radius, volume fraction}
### Remark : if a packed configuration is desired, put 1 as the volume fraction of the smallest radius
desiredRPhi = [[1., 0.1], [1.,0.05], [0.5, 1]]
### Assign phases to the desired (r,phi)
### Remark desiredRPhi[i] corresponds to phases[i]
phases = [0, 1, 0]
###
dico = { 0:"ph1", 1:"ph2"}
### Select the exclusion distance between non-intersecting spheres
exclusionDistance = 0.1
### Select the boundary exclusion distance between spheres and boundary of the shape
boundaryExclusionDistance = 4
### Select the seed for random engine
seed = 2
### Select the algorithm (1 or 2, 1 is recommended)
### Among: 1, 2 (1 is recommended)
method = 1


###############################################################################
### Standard way for 3D
###############################################################################
## Parametrize
###################
algo_3D = sac_de_billes.AlgoRSA_3D()                          			### build the object
algo_3D.setExclusionDistance(exclusionDistance)     			### optional
algo_3D.setBoundaryExclusionDistance(boundaryExclusionDistance)		### optional
algo_3D.setBigShape(L, shape)                       			### 
algo_3D.setRadiusGenerator(desiredRPhi, phases)     			### define the radii of the spheres to be thrown, from the desired{radius, volume fraction}
algo_3D.setNamePhase(dico)                          			### optional, for writing the phases names as string in the outputs
### may be replaced by below :
### algo_3D.setRadiusGenerator(table_of_radii, table_of_phases) ### take directly the radii of the spheres, and the related phases as input
### algo_3D.setRadiusGenerator(name_file)           ### reads directly the radii of spheres from the file name_file (ASCII, 1 or 2 double per line)
###################
## Launch the algo
###################
algo_3D.proceed(seed, method)
###################
## Get the output
###################
algo_3D.printDump("sortie3D.dump")                    ### output for Ovito
algo_3D.printCSV("sortie3D.csv")                      ### output ASCII csv
algo_3D.printPos("sortie3D.pos")                      ### output for Combs
the_spheres = algo_3D.getPlacedSpheres()            ### standard python array [center[0], center[1], center[2], radius]
###############################################################################
### END
###############################################################################

###############################################################################
### Standard way for 2D
###############################################################################
L_2D = [64., 64.]
## Parametrize
###################
algo_2D = sac_de_billes.AlgoRSA_2D()                          ### builds the object
algo_2D.setExclusionDistance(exclusionDistance)     ### optional
algo_2D.setBigShape(L_2D, shape)                    ### L should be 2D
### may be replaced by below :
### algo_2D.setRadiusGenerator(table_of_radii, table_of_phases)  ### takes directly the radii of the spheres as input
### algo_2D.setRadiusGenerator(desiredRPhi, phases)     ### defines the radii of the spheres to be thrown
algo_2D.setRadiusGenerator("File_radius_phase.txt")           ### reads directly the radii of spheres from the file name_file (ASCII, 1 or 2 double per line)
###################
## Launch the algo
###################
algo_2D.proceed(seed, method)
###################
## Get the output
###################
algo_2D.printDump("sortie2D.dump")                    ### output for Ovito
algo_2D.printCSV("sortie2D.csv")                      ### output ASCII csv
the_spheres = algo_2D.getPlacedSpheres()            ### standard python array [center[0], center[1], center[2], radius]
###############################################################################
### END
###############################################################################
print("Fin phase 4")


### Compact way, will be deprecated
algo3D = sac_de_billes.AlgoRSA3D()
algo3D = sac_de_billes.AlgoRSA3D(L, desiredRPhi, exclusionDistance, seed, method, shape)




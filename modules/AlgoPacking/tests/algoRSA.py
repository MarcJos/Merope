### import the library (should accessible through PYTHON_PATH)
import sac_de_billes as sac_de_billes

### Select the big shape in which spheres are thrown
### Among : Tore, Cube (=Cuboid), Sphere, Cylinder
nameShape = sac_de_billes.NameShape.Tore
### Select the dimensions of the cuboid (sharply) containing the shape
### At the moment, indirectly defines the size of the shape for 
##### Sphere via : L = [2 * radius, 2 * radius, 2 * radius]
##### Cylinder via : L = [2 * radius, 2 * radius, height]
L = [32., 32., 32.]
### Select the desired {radius, volume fraction}
### Remark : if a packed configuration is desired, put 1 as the volume fraction of the smallest radius
desiredRPhi = [[1., 0.2], [0.5, 1]]
## Table of phases (each phase corresponding to each pair (radius, volume_fraction) of desiredRPhi
tabPhases = [0,1]
### Select the exclusion distance between non-intersecting spheres
exclusionDistance = 0.0
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
algo_3D.setBigShape(L, nameShape)                       			### 
algo_3D.setRadiusGenerator(desiredRPhi, tabPhases)   			### define the radii of the spheres to be thrown, from the desired{radius, volume fraction}. tabPhases is an optional argument
### may be replaced by below :
### algo_3D.setRadiusGenerator(table_of_radii)      ### take directly the radii of the spheres as input
### algo_3D.setRadiusGenerator(name_file)           ### read directly the radii of spheres from the file name_file (ASCII, 1 double per line)
###################
## Launch the algo
###################
output_message = algo_3D.proceed(seed, method)      ### returns a map
print(output_message)
###################
## Get the output
###################
algo_3D.printDump("sortie3D.dump")                    ### output for Ovito
algo_3D.printCSV("sortie3D.csv")                      ### output ASCII csv
algo_3D.printPos("sortie3D.pos")                      ### output for Combs
theSpheres = algo_3D.getSpheres()                     ### standard python array of sphere objects
###############################################################################
### END
###############################################################################

### Compact way
typeAlgo = sac_de_billes.TypeAlgo.RSA
theSpheres = sac_de_billes.throwSpheres_3D(typeAlgo, nameShape, L, seed, desiredRPhi, tabPhases, exclusionDistance)


###############################################################################
### Standard way for 2D
###############################################################################
L_2D = [64., 64.]
## Parametrize
###################
algo_2D = sac_de_billes.AlgoRSA_2D()                          ### builds the object
algo_2D.setExclusionDistance(exclusionDistance)     ### optional
algo_2D.setBoundaryExclusionDistance(boundaryExclusionDistance)		### optional
algo_2D.setBigShape(L_2D, nameShape)                       ### L should be 2D
algo_2D.setRadiusGenerator(desiredRPhi)             ### defines the radii of the spheres to be thrown
### may be replaced by below :
### algo_2D.setRadiusGenerator(table_of_radii)      ### takes directly the radii of the spheres as input
### algo_2D.setRadiusGenerator(name_file)           ### reads directly the radii of spheres from the file name_file (ASCII, 1 double per line)
###################
## Launch the algo
###################
output_message = algo_2D.proceed(seed, method)      ### returns a map
print(output_message)
###################
## Get the output
###################
algo_2D.printDump("sortie2D.dump")                    ### output for Ovito
algo_2D.printCSV("sortie2D.csv")                      ### output ASCII csv
theSpheres = algo_2D.getSpheres()                     ### standard python array of sphere objects
###############################################################################
### END
###############################################################################





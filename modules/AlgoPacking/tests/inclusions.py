### This python script generates a set of inclusions within
### a periodic cube. The inclusions are translated in order
### to guarantee that the origin (0,0,0) is not within an
### inclusion and at a minimal distance of one inclusion.
###
### Once you have launched this script, you should get a file
### named `Translated_spheres_3D.csv`.
### You can rename this file `sortie_csv.csv` and move it to
### `../Ecriture_gmsh directory`. Then, in executing the other
### python script `From_RSA_2_gmsh.py` you will get a gmsh input
### file that will allow you to generate a mesh.


### import the library sac_de_billes*.so (should accessible through PYTHON_PATH)
import sac_de_billes

## Select the big shape in which spheres are thrown
### Among : Tore, Cube (=Cuboid), Sphere, Cylinder
shape = sac_de_billes.NameShape.Tore

### Select the dimensions of the cuboid (sharply) containing the shape
### At the moment, indirectly defines the size of the shape for
##### Sphere via : L = [2 * radius, 2 * radius, 2 * radius]
##### Cylinder via : L = [2 * radius, 2 * radius, height]
L = [1., 1., 1.]

### Select the desired {radius, volume fraction}
### Remark : if a packed configuration is desired, put 1 as the volume fraction of the smallest radius
desiredRPhi = [[0.1, 0.05], [0.06, 0.2], [0.03, 1]]

### Select the exclusion distance between non-intersecting spheres
exclusionDistance = 0.0

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

algo_3D = sac_de_billes.AlgoRSA_3D()                          ### build the object
algo_3D.setExclusionDistance(exclusionDistance)     ### optional
algo_3D.setBigShape(L, shape)                       ###
algo_3D.setRadiusGenerator(desiredRPhi)             ### define the radii of the spheres to be thrown, from the desired{radius, volume fraction}
algo_3D.proceed(seed, method)

###################
## Get the output
###################

algo_3D.printDump("sortie3D.dump")                    ### output for Ovito
algo_3D.printCSV("sortie3D.csv")                      ### output ASCII csv

the_spheres = algo_3D.getPlacedSpheres()            ### standard python array [center[0], center[1], center[2], radius]

sphereManip_3D = sac_de_billes.SphereManipulator_3D(algo_3D.getPlacedSpheres(), algo_3D.getLength())
print("You cannot achieve a distmin better than")
print(sphereManip_3D.upper_bound_on_best_distmin())

### Randomly translates (periodically) the spheres to find a better configuration
print("\n Translating the spheres...")
print("Achieved distmin : ")
print(sphereManip_3D.random_search(20000))
print(" ... spheres translated\n")
### Prints the result
sphereManip_3D.printCSV("Translated_spheres_3D.csv")
sphereManip_3D.printDump("Translated_spheres_3D.dump")

 

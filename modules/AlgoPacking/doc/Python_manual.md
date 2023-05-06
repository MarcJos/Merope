[[_TOC_]]

**Remark** : All the Python commands correspond 1 to 1 to the C++ commands. Consider reading the doxygen documentation for more details.

**Remark** : All the classes and functions can be subscripted `_3D` or `_2D` depending on the dimension of the considered space.

# Simplified mode
The functions below return a list of spheres
- `throwSpheres_3D(typeAlgo, nameShape, L, seed, desiredRPhi, tabPhases, mindist)` : return a list of spheres randomly thrown.
- `throwSpheres_3D(typeAlgo, nameShape, L, seed, tabRadii, tabPhases, mindist)` same as above, but we specify the radii list instead of the volume fractions.
- `fillMaxRSA_3D(nameShape, L, NbSpheres, seed, mindist)` :
Tries to fill a box with a given number of spheres, such that the resulting distribution looks like the RSA algorithm. (Hence, it adjusts the raidii of the desired spheres, these being close to the ones producing the theoretical volume fraction.)
It goes as follows :
	- fix the `radius` of spheres such that the volume fraction of `NbSpheres` of spheres of radius `radius + 0.5 mindist` inside the shape is equal to the theoretical volume fraction of the RSA packing
	- while the desired number of spheres is not achieved : Throw the RSA algorithm, stopping when  `NbSpheres` spheres are placed or at full packing. If `NbSpheres` of spheres are placed, return them. Else, lower the radius by `radius *= theta` (`theta = 0.95` but should be verified in [Interface.ixx](../src/Interface.ixx)). Repeat it until `NbSpheres` of spheres are placed.

:warning: `fillMaxRSA_3D` is practical when one wants a precise number of spheres, but does not make much sense mathematically.

Parameters :
- `typeAlgo` : `TypeAlgo` object
- `L` : lengths of the shape
- `nameShape`: `NameShape` object, defining the global shape in which the spheres are thrown
- `seed` : integer, that initializes the random generator
- `desiredRPhi` : vector of vectors of size 2 `[radius,volume fraction]`
- `tabPhases` : table of phases, each corresponding to a pair `[radius,volume fraction]` in `desiredRPhi`
- `mindist` : minimal distance between two spheres (if applicable)
- `tabRadii` : list of radii


# The pipeline

Some examples are available in [`tests`](../tests).
Here, we describe precisely each step and each feature.

# Parametrize

- The *order* in parametrization matters.
- Once fixed, arguments should not be changed.
- An object "algorithm" can be launched only once.

## Choose the algorithm

This is a constructor.
Choose simultaneously the dimension of the space (2D or 3D) and the algorithm among the RSA, WP, and Boolean algorithms.

`algo = rsa.AlgoRSA_3D()`  
`algo = rsa.AlgoRSA_2D()`  
`algo = rsa.AlgoWP_3D()`  
`algo = rsa.AlgoWP_2D()`  
`algo = rsa.AlgoBool_3D()`  
`algo = rsa.AlgoBool_2D()`

An *enum class* represent the type of the algorithm : `TypeAlgo`. Possibilities : `RSA`, `WP`, `BOOL`.

## Exclusion distances [optional]

### Exclusion distance between spheres

`algo.setExclusionDistance(d_0)`  
For an exclusion distance $d_0$, the algorithms *guarantee* that : the minimal distance between two sphere is larger than $d_0$. It is possible to ask for negative exclusion distance. The result will make sense as prescribing a maximum interpenetration distance.

**Default:** $d_0=0$.

:warning: If $d_0 < - r/2$ (notice the **minus** sign!), where $r$ is the minimal radius of the spheres, this will lead to undefined behaviors.

### Exclusion distance between the spheres and the boundary

`algo.setBoundaryExclusionDistance(d_1)`
For an exclusion distance $d_1$, the algorithms *guarantee* that : if the surrounding domain is not the torus, the minimal distance between any sphere and the boundary is larger than $d_1$.

**Default:** $d_1=0.5*d_0$, where $d_0$ is the exclusion distance between two spheres.


## Set the surrounding shape

`algo.setBigShape(L, shape)`  
Parameters:
- `L` is a 3 or 2-dimensional vector of double. It encodes the size of the surrounding box (a cuboid with origin $0$) in each direction.
- `shape` is a string, that should be chosen among `Tore` (torus, a cuboid periodic in every direction), `Cube` (actually, a cuboid, for each length may be fixed), `Sphere`, `Cylinder` (only 3D).

For the Sphere and the Cylinder, the dimensions of the surrounding shape are infered from `L` as follows :
- for the `Sphere` via : L = [2 * radius, 2 * radius, 2 * radius],
- for the `Cylinder` via : L = [2 * radius, 2 * radius, height].

An *enum class* specifies the global shape in which the spheres are thrown : `NameShape`. Possibilities : `Tore` (torus = periodic cuboid), `Cube` (cuboid), `Sphere`, `Cylinder`.

## Set the radius generator

During the algorithm, spheres are randomly placed in the shape.
The (double) *radius* and the (int) *phase identifier* of the spheres are provided by the `RadiusGenerator`.
The latter can be set by two ways: 
- by giving it a fixed sequence of radii, 
- by requesting a list of radii associated with volume fractions.

**Default:** giving the phase is *optional*. By default, it is assumed to be equal to $0$.

### Via a fixed sequence of radii

Via a vector :  
`algo.setRadiusGenerator(table_of_radii, table_of_phases)` (requires len(table_of_radii) == len(table_of_phases))  
`algo.setRadiusGenerator(table_of_radii)`  (default, phase == 0)  

Via an ASCII file :  
`algo.setRadiusGenerator(nameFile)`  
Two possibilities:  
- either a radius on each line (default, phase == 0),
- or a radius and a phase on each line, separated by a blank space.

### Via a desired volume fraction for radii

`algo.setRadiusGenerator(desiredRPhi, table_of_phases)` (requires len(desiredRPhi) == len(table_of_phases))  
`algo.setRadiusGenerator(desiredRPhi)`   (default, phase == 0)  

`desiredRPhi` is a vector of vector vectors of size 2 `[radius,volume fraction]`.

*Remark:* The prescribed volume fractions will be approached by below, but are not exactly preserved.

:warning: The algorithm computes the volume fraction by summing up all the volumes of spheres. This computation is **false** if intersections are allowed (either by negative minimal distance or using the Boolean algorithm).

### Tag the spheres

The usual way is to tag the spheres with an `int` via the mechansim using `table_of_phase` described above.
However, it is possible as well to tag them with a `string` via the function `setNamePhase` and a dictionnary relating the `table_of_phase` indices with some `string`s.
For example :  
    `desiredRPhi_2 = [[5.,0.15], [4.,0.05]]`  
	`tabPhases = [5,36]`  
	`dico = {5 : "ph1", 36 : "ph2"}`    
	`algo.setRadiusGenerator(desiredRPhi_2,tabPhases)`  
    `algo.setNamePhase(dico)`


# Launch the algo

`output_message = algo.proceed(seed, method)`  

(int) `seed` is a random seed.

`method` $\in \{1,2\}$ chooses implementation details of the RSA or WP algorithm (AlgoBool is unsensitive to this parameter).
This parameter is exclusively for developpers.
For RSA, the final results are identical (but one is slower), but not for the WP algorithm.
In general, `method=1` is a good choice.

`output_message` : a map of messages given by the algorithm. Contains in particular whether the configuration is `Packed' or not.

# Get the output

## The spheres

The classes `Sphere_3D` and `Sphere_2D` have 3 properties:
- `center` : a point in $\mathbb{R}^3$,
- `radius` : a positive scalar parameter,
- `phase` : a nonnegative integer.

## Print
`algo.printDump("sortie3D.dump")`  return an .dump file (LAMMPS format), that can be read by Ovito,  
`algo.printCSV("sortie3D.csv")`  return a .csv file, with " , " as separators,  
`algo.printPos("sortie3D.pos")`  return a .pos file, that can be used in COMBS,  
`the_spheres = algo.getPlacedSpheres()` return a python array where the spheres are encoded in 5 components [center.x, center.y, center.z, radius, phase] (in 2D, center.z == 0),  
`the_phases =  algo.getPhases()` return the list of phase related to the list of placedSpheres.

## Manipulate the spheres

:construction:

:warning: Add a method setSpheres, and think about how to articulate the classes.

# Behaviors when all spheres cannot be placed
The two algorithms react differently.
- The RSA tries to place as much spheres as possible. If it cannot achieve the request, it says "Fully packed" and stops. Hence, it is interesting to ask for too much spheres (e.g. such that the total requested volume fraction is 1).
- If the WP algorithm cannot place all the desired spheres, it returns an error.

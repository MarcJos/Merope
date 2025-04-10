# Voro++

## Useful information

See [voro++ general overview](http://math.lbl.gov/voro++/about.html) and  [voro++ reference manual](http://math.lbl.gov/voro++/doc/refman/index.html).

### Display the voronoi structure with drawGnuPlot

Use the function `voro::container_poly::draw_cells_gnuplot` from voro++, see [voro++ examples](http://math.lbl.gov/voro++/examples/single_cell/) and [voro++ manual](http://math.lbl.gov/voro++/doc/refman/classvoro_1_1container__poly.html).

The python command is  
`voroInterface.drawGnuPlot("sortie_gnuplot.txt")`  

Then, use `gnuplot` and type `splot 'sortie_gnuplot.txt' with lines`.


### Python manual (through Mérope)

**Main classes and methods** :
- `merope.VoroInterface_3D` : python wrapper for building Laguerre tessellations
    - `merope.VoroInterface_3D(L, centerTessels, periodicity = [True, True, True])` : builds a tessellations from centers of tessels with Laguerre weights.
        - `L` : of type List(float), dimensions of the cube 
        - `centerTessels`: of type List(Sphere_3D), list of spheres representing centers of tessels (the Laguerre weight is equal to the square of the sphere radius)
        - `periodicity`: of type List(bool), periodicity in each direction
    - `drawGnuPlot(nameFile)` : draw a gnuplot representation of the polycrystal
    - `drawCellsPov(nameFile)` : draw a POV-Ray representation of the polycrystal
    - `printCustom(format, nameFile)` : draw a custom output of voro++ (see  [voro++ manual](https://math.lbl.gov/voro++/doc/custom.html))
    - `addWallCylinder(xcenter, ycenter, zcenter, xvector, yvector, zvector, radius)` : change container shape to cylinder (periodicity = False)
    - `computeSolids()` : compute and regularize the geometry to make coincidence between all the vertex of the cells. Return a list of solids defined as an identifier + a list of face identifiers, and a list of faces given by an identifier + a halfspace
    - `getCellCenters()` : get the position of the center of the cells as a list defined by an identifier and a point3D
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
- `merope.VoroInterface_3D` : python wrapper for building Laguerre tessellations in (periodic cubes)
    - `merope.VoroInterface_3D(L, centerTessels, periodicity = [True, True, True])` : builds a tessellations from centers of tessels with Laguerre weights.
        - `L` : of type List(float), dimensions of the cube 
        - `centerTessels`: of type List(Sphere_3D), list of spheres representing centers of tessels (the Laguerre weight is equal to the square of the sphere radius)
        - `periodicity`: of type List(bool), periodicity in each direction
    - `drawGnuPlot(nameFile)` : draw a gnuplot representation of the polycrystal
    - `drawCellsPov(nameFile)` : draw a POV-Ray representation of the polycrystal
    - `printCustom(format, nameFile)` : draw a custom output of voro++ (see  [voro++ manual](http://math.lbl.gov/voro++/doc/refman/classvoro_1_1container__poly.html))

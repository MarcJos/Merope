# Voro++

## Useful information

See [voro++ general overview](http://math.lbl.gov/voro++/about.html) and  [voro++ reference manual](http://math.lbl.gov/voro++/doc/refman/index.html).

### Display the voronoi structure with drawGnuPlot

Use the function `voro::container_poly::draw_cells_gnuplot` from voro++, see [voro++ examples](http://math.lbl.gov/voro++/examples/single_cell/) and [voro++ manual](http://math.lbl.gov/voro++/doc/refman/classvoro_1_1container__poly.html).

The python command is  
`voroInterface.drawGnuPlot("sortie_gnuplot.txt")`  

Then, use `gnuplot` and type `splot 'sortie_gnuplot.txt' with lines`.

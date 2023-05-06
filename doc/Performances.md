# Voxelation efficiency

## Comparison with tmfft/VER

We compare Mérope and tmfft/VER when reconstructing a voxelation of a Voronoi tessellation, varying :
- the number nVox of voxels,
- the number nbSpheres of seeds (=of tessels).
We observe that, when the number of spheres rescaled by the volume of a single sphere is smaller than the number of voxels, the execution time of Mérope is quite unsensitive to the number of spheres, whereas the execution time of tmfft/VER scales linearly with the number of spheres.
Using Mérope instead of tmfft, we may gain a factor from 1 to 100 in terms of time.

See [benchmark.py](studies/performances/Merope_vs_tmfft/benchmark.py)

<img src="studies/performances/Merope_vs_tmfft/Benchmark_Merope_TMFFT_2021-12-16.png" alt="drawing" width="1000"/>

<img src="studies/performances/Merope_vs_tmfft/Benchmark_Merope_TMFFT_speedup_2021-12-16.png" alt="drawing" width="1000"/>

## Comparison with Rollpy

We compare Mérope and tmfft/VER when reconstructing a voxelation of spherical inclusions inside a matrix.
In this case, Mérope is faster than rollpy, by a factor 10-100.

See [test_performance.py](studies/performances/Merope_vs_RollPy/test_performance.py)

<img src="studies/performances/Merope_vs_RollPy/Time.png" alt="drawing" width="1000"/>

<img src="studies/performances/Merope_vs_RollPy/Speed-up.png" alt="drawing" width="1000"/>


## Comparison with Neper

We compare Mérope and Neper when reconstructing a voxelation of a Laguerre tessellation, varying :
- the number nVox of voxels,
- the number nbSpheres of seeds (=of tessels).
We observe that, using Mérope instead of Neper, we may gain a factor from 1 to 100 in terms of time. (Typically more than 10.) However, Neper supports parallel computing, which is not the case yet with Mérope.

See [Calc_benchmark.py](studies/performances/Merope_vs_neper/Calc_benchmark.py). (1 processor for Neper.)

<img src="studies/performances/Merope_vs_neper/Comparaison.png" alt="drawing" width="1000"/>

<img src="studies/performances/Merope_vs_neper/SpeedUp.png" alt="drawing" width="1000"/>


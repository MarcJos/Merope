import common
import merope

# Import py4amitex
import math as m

import sac_de_billes
import merope

from py4amitex.input import (
    Input,
    Grid,
    Algorithm,
    AlgorithmParameters,
    Diffusion,
    Materials,
    Material,
    DiffusionDriving,
    Evolution,
    Loading,
    LoadingOutput,
    Output,
    Zone,
    ReferenceMaterialD,
)

from py4amitex.input.materialbuilder import buildMaterials, VoxelSpec, IndexOrdering
from py4amitex.simulation import runSimulationExternal
from py4amitex.output.amitexoutput import AmitexOutput

# 1 Define algorithm
def makeAlgoParams():
    # Set algorithm parameters
    algorithm = Algorithm(type="Basic_Scheme", convergenceAcceleration=True)
    # Define a mechanical resolution
    diffusion = Diffusion(filter="Default",  stationary = True)
    algorithmParameters = AlgorithmParameters(algorithm, diffusion=diffusion)
    return algorithmParameters

# 2 Define loading
def makeLoadings():
    loadingOutput = LoadingOutput()

    loading = Loading()
    # 1 timesteps to final time 1
    loading.setTimeDiscretizationLinear(1, 1.0)
    # Applay strain for XX,YY,ZZ components
    for i in range(3):
        loading.setEvolution(i, DiffusionDriving.Gradient, Evolution.Linear, common.direction[i])
    # Output
    loading.setOutputVtkList([1])
    loadingOutput.add(loading)
    output = Output()
    output.setVtkFluxDGradD(fluxD=1, gradD=1)
    loadingOutput.output = output
    #
    return loadingOutput

# 3 Materials
def makeMaterials(phaseGrid):    
    materials = Materials()
    materials.referenceMaterialD = ReferenceMaterialD(common.therm_ref_coeff)    
    #################
    mat = Material()
    mat.setLawK("Fourier_iso_polarization")
    mat.setNumberCoeffK(4)
    print(phaseGrid)
    for coeff in phaseGrid:
        # Empty zone, to be filled by buildMaterials
        print(coeff)
        mat.addZone(Zone(GRID_DIMS), coeffKs=[coeff, 0, 0, -coeff])
    materials.add(mat)
    # A bit more involved with a type VoxelSpec with zone specifications:
    # every one is material 0, but phase index is taken as a zone index
    # [phaseId] -> [(0,Φ=1,zoneId=phaseId)]
    voxels = [VoxelSpec([(0, 1.0, i)]) for i, p in enumerate(phaseGrid)]
    buildMaterials(materials, GRID_DIMS, voxels, IndexOrdering.C)
    #################
    return materials


def run(GRID_DIMS, L, numberProcs, structure, homogenizationRule):
    # Define all input categories
    grid = create_homogenized_grid(structure, L, GRID_DIMS, homogenizationRule)
    algorithmParameters = makeAlgoParams()
    loadingOutput = makeLoadings()
    materials = makeMaterials(grid.get_as_list())
    grid_amitex = Grid(GRID_DIMS, [L[i] / GRID_DIMS[i] for i in range(3)])
    input = Input(grid_amitex, algorithmParameters, materials, loadingOutput)
    # Define output.std in "amitex_dir"
    input.resultsDir = "amitex_dir_therm"
    # Run amitex_fftp on two processes
    runSimulationExternal(input, numberProcs)

def create_homogenized_grid(structure, L, GRID_DIMS, homogenizationRule):
    grid = common.create_grid(structure, L, GRID_DIMS, merope.vox.VoxelRule.PolyGeom)
    grid.apply_coefficients(common.therm_pure_coeffs)
    grid.convert_to_Iso_format()
    grid.apply_homogRule(homogenizationRule)
    return grid



if __name__=="__main__":    
    GRID_DIMS = [4 for _ in range(3)]
    L = [1, 1, 1]
    numberProcs = 4
    homogenizationRule = merope.HomogenizationRule.Voigt
    structure = common.create_structure(L)
    run(GRID_DIMS, L, numberProcs, structure, merope.HomogenizationRule.Voigt)
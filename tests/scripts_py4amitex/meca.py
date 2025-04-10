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
    Mechanics,
    Materials,
    Material,
    MechanicDriving,
    Evolution,
    Loading,
    LoadingOutput,
    Output,
    Zone,
    ReferenceMaterial,
)
from py4amitex.input.materialbuilder import buildMaterials, VoxelSpec, IndexOrdering
from py4amitex.simulation import runSimulationExternal
from py4amitex.output.amitexoutput import AmitexOutput

# 1 Define algorithm
def makeAlgoParams():
    # Set algorithm parameters
    algorithm = Algorithm(type="Basic_Scheme", convergenceAcceleration=True)
    # Define a mechanical resolution
    mechanics = Mechanics(filter="Default", smallPerturbations=True)
    algorithmParameters = AlgorithmParameters(algorithm, mechanics=mechanics)
    return algorithmParameters

# 2 Define loading
def makeLoadings():
    loadingOutput = LoadingOutput()

    loading = Loading()
    # 1 timesteps to final time 1
    loading.setTimeDiscretizationLinear(1, 1.0)
    # Applay strain for XX,YY,ZZ components
    for i in range(3):
        loading.setEvolution((i, i), MechanicDriving.Strain, Evolution.Linear, common.direction[i])
    for i in range(3):
        for j in range(i + 1, 3):
            loading.setEvolution((i, j), MechanicDriving.Strain, Evolution.Linear, 0.0)
    # Output
    loading.setOutputVtkList([1])
    loadingOutput.add(loading)
    output = Output()
    output.setVtkStressStrain(stress=1, strain=1)
    output.setVtkStressStrain(stress=1, strain=1)
    loadingOutput.output = output
    #
    return loadingOutput

# 3 Materials
def makeMaterials(phaseGrid, amitex_law):
    materials = Materials()
    nPhases = common.number_phases
    materials.referenceMaterial = ReferenceMaterial(common.meca_ref_coeffs(0), common.meca_ref_coeffs(1))

    for m in range(nPhases):
        mat = Material()
        mat.setLaw("elasiso")
        mat.setCoeffs(common.meca_pure_coeffs(m))
        mat.setCoeffComposites(common.meca_pure_coeffs(m))
        materials.add(mat)
    
    buildMaterials(materials, GRID_DIMS, phaseGrid, ordering=IndexOrdering.C)
    materials.composite(0).setLaw(amitex_law)
    return materials


def run(GRID_DIMS, L, numberProcs, structure, amitex_law):
    # Define all input categories
    grid = common.create_grid(structure, L, GRID_DIMS, voxelRule)
    algorithmParameters = makeAlgoParams()
    loadingOutput = makeLoadings()
    materials = makeMaterials(grid.get_as_list(), amitex_law)
    grid_amitex = Grid(GRID_DIMS, [L[i] / GRID_DIMS[i] for i in range(3)])
    input = Input(grid_amitex, algorithmParameters, materials, loadingOutput)
    # Define output.std in "amitex_dir"
    input.resultsDir = "amitex_dir_meca"
    # Run amitex_fftp on two processes
    runSimulationExternal(input, numberProcs)


if __name__=="__main__":    
    GRID_DIMS = [64 for _ in range(3)]
    L = [1, 1, 1]
    numberProcs = 4
    amitex_law = "laminate"
    voxelRule = merope.vox.VoxelRule.Laminate
    structure = common.create_structure(L)
    run(GRID_DIMS, L, numberProcs, structure, amitex_law)

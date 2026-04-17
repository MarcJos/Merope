# #!/usr/bin/env python3
# """
# scripts/generate_examples.py
# =============================
# Génère les scripts Mérope, Combs et sessions JSON de référence.
# Appelé par : make example
# """
# import os, sys
# sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

# from microgen.models import (
#     MicrostructureConfig, SimulationBox, SphereFamily, MultiSphereFamily,
#     LaguerreTessFamily, VoxelGrid, PackingAlgo, GeneratorBackend, GlobalShape,
# )

# OUT = os.path.join(os.path.dirname(__file__))
# os.makedirs(OUT, exist_ok=True)

# # ── 1. Béton — 2 familles RSA ────────────────────────────────────────────
# beton = MicrostructureConfig(
#     name="beton_VER_RSA",
#     description="Béton : 2 familles de granulats sphériques (RSA, Mérope)",
#     backend=GeneratorBackend.MEROPE,
#     box=SimulationBox(Lx=1e-3, Ly=1e-3, Lz=1e-3, shape=GlobalShape.TORE),
#     inclusions=[
#         SphereFamily(radius=5e-5, volume_fraction=0.15, phase=1, algo=PackingAlgo.RSA),
#         SphereFamily(radius=2e-5, volume_fraction=0.10, phase=2, algo=PackingAlgo.RSA),
#     ],
#     voxel_grid=VoxelGrid(nx=128, ny=128, nz=128, output_vtk="beton_VER.vtk"),
#     seed=42,
# )
# open(f"{OUT}/merope_beton.py",    "w").write(beton.to_merope_script())
# open(f"{OUT}/combs_beton.py",     "w").write(beton.to_combs_script())
# open(f"{OUT}/beton_VER.microgen", "w").write(beton.to_json())
# print(f"[beton]    frac. vol. estimée : {beton.volume_fraction_estimate():.2%}")
# print(f"           voxel : {beton.voxel_size_um()[0]:.3f} µm")

# # ── 2. Distribution multi-rayons ─────────────────────────────────────────
# multi = MicrostructureConfig(
#     name="granulat_multi_rayons",
#     description="Béton : distribution granulométrique (sac_de_billes multi-rayons)",
#     backend=GeneratorBackend.MEROPE,
#     box=SimulationBox(Lx=2e-3, Ly=2e-3, Lz=2e-3, shape=GlobalShape.TORE),
#     inclusions=[
#         MultiSphereFamily(
#             distribution=[
#                 (1e-4, 0.10, 1),   # gros granulats
#                 (5e-5, 0.08, 2),   # granulats moyens
#                 (2e-5, 0.06, 3),   # fins
#             ],
#             algo=PackingAlgo.RSA,
#         )
#     ],
#     voxel_grid=VoxelGrid(nx=128, ny=128, nz=128, output_vtk="multiradius.vtk"),
#     seed=99,
# )
# open(f"{OUT}/merope_multiradius.py", "w").write(multi.to_merope_script())
# print(f"[multi]    frac. vol. estimée : {multi.volume_fraction_estimate():.2%}")

# # ── 3. Polycristal Laguerre ──────────────────────────────────────────────
# laguerre = MicrostructureConfig.model_construct(
#     name="polycristal_acier",
#     description="Acier polycristallin — Laguerre tessellation (voro++)",
#     backend=GeneratorBackend.MEROPE,
#     box=SimulationBox(Lx=5e-4, Ly=5e-4, Lz=5e-4, shape=GlobalShape.TORE),
#     inclusions=[LaguerreTessFamily(n_grains=100)],
#     voxel_grid=VoxelGrid(nx=64, ny=64, nz=64, output_vtk="polycristal.vtk"),
#     seed=7,
# )
# open(f"{OUT}/merope_polycristal.py", "w").write(laguerre.to_merope_script())
# print(f"[laguerre] {100} grains — voxel : {laguerre.voxel_size_um()[0]:.3f} µm")

# print(f"\nScripts générés dans {OUT}/")
# print("  merope_beton.py, combs_beton.py, beton_VER.microgen")
# print("  merope_multiradius.py")
# print("  merope_polycristal.py")




#!/usr/bin/env python3
"""
scripts/generate_examples.py
=============================
Génère les scripts Mérope, Combs et sessions JSON de référence.
Appelé par : make example
Options : --beton, --multi, --laguerre, --all (défaut)
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from microgen.models import (
    MicrostructureConfig, SimulationBox, SphereFamily, MultiSphereFamily,
    LaguerreTessFamily, VoxelGrid, PackingAlgo, GeneratorBackend, GlobalShape,
)

OUT = os.path.join(os.path.dirname(__file__))
os.makedirs(OUT, exist_ok=True)


def generate_beton():
    beton = MicrostructureConfig(
        name="beton_VER_RSA",
        description="Béton : 2 familles de granulats sphériques (RSA, Mérope)",
        backend=GeneratorBackend.MEROPE,
        box=SimulationBox(Lx=1e-3, Ly=1e-3, Lz=1e-3, shape=GlobalShape.TORE),
        inclusions=[
            SphereFamily(radius=5e-5, volume_fraction=0.15, phase=1, algo=PackingAlgo.RSA),
            SphereFamily(radius=2e-5, volume_fraction=0.10, phase=2, algo=PackingAlgo.RSA),
        ],
        voxel_grid=VoxelGrid(nx=128, ny=128, nz=128, output_vtk="beton_VER.vtk"),
        seed=42,
    )
    open(f"{OUT}/merope_beton.py", "w").write(beton.to_merope_script())
    open(f"{OUT}/combs_beton.py", "w").write(beton.to_combs_script())
    open(f"{OUT}/beton_VER.microgen", "w").write(beton.to_json())
    print(f"[beton]    frac. vol. estimée : {beton.volume_fraction_estimate():.2%}")
    print(f"           voxel : {beton.voxel_size_um()[0]:.3f} µm")


def generate_multi():
    multi = MicrostructureConfig(
        name="granulat_multi_rayons",
        description="Béton : distribution granulométrique (sac_de_billes multi-rayons)",
        backend=GeneratorBackend.MEROPE,
        box=SimulationBox(Lx=2e-3, Ly=2e-3, Lz=2e-3, shape=GlobalShape.TORE),
        inclusions=[
            MultiSphereFamily(
                distribution=[
                    (1e-4, 0.10, 1),
                    (5e-5, 0.08, 2),
                    (2e-5, 0.06, 3),
                ],
                algo=PackingAlgo.RSA,
            )
        ],
        voxel_grid=VoxelGrid(nx=128, ny=128, nz=128, output_vtk="multiradius.vtk"),
        seed=99,
    )
    open(f"{OUT}/merope_multiradius.py", "w").write(multi.to_merope_script())
    print(f"[multi]    frac. vol. estimée : {multi.volume_fraction_estimate():.2%}")


def generate_laguerre():
    laguerre = MicrostructureConfig.model_construct(
        name="polycristal_acier",
        description="Acier polycristallin — Laguerre tessellation (voro++)",
        backend=GeneratorBackend.MEROPE,
        box=SimulationBox(Lx=5e-4, Ly=5e-4, Lz=5e-4, shape=GlobalShape.TORE),
        inclusions=[LaguerreTessFamily(n_grains=100)],
        voxel_grid=VoxelGrid(nx=64, ny=64, nz=64, output_vtk="polycristal.vtk"),
        seed=7,
    )
    open(f"{OUT}/merope_polycristal.py", "w").write(laguerre.to_merope_script())
    print(f"[laguerre] 100 grains — voxel : {laguerre.voxel_size_um()[0]:.3f} µm")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--beton", action="store_true", help="Générer l'exemple béton")
    parser.add_argument("--multi", action="store_true", help="Générer l'exemple multi-rayons")
    parser.add_argument("--laguerre", action="store_true", help="Générer l'exemple polycristal Laguerre")
    parser.add_argument("--all", action="store_true", help="Générer tous les exemples")
    args = parser.parse_args()

    # Si aucun argument spécifique n'est donné, générer tous les exemples
    if not (args.beton or args.multi or args.laguerre or args.all):
        args.all = True

    if args.all or args.beton:
        generate_beton()
    if args.all or args.multi:
        generate_multi()
    if args.all or args.laguerre:
        generate_laguerre()

    print(f"\nScripts générés dans {OUT}/")
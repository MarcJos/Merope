"""
tests/test_models.py
====================
Tests unitaires du modèle pydantic v2 MICROGEN.

Couvre :
  - Construction et validation des modèles
  - Génération de code Mérope / Combs
  - Sérialisation / désérialisation JSON
  - Cas invalides (validation croisée)
  - API Mérope : noms de classes, signatures correctes
"""

import json
import math
import sys
import os
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from microgen.models import (
    SimulationBox, GlobalShape, VoxelGrid,
    SphereFamily, MultiSphereFamily, LaguerreTessFamily, EllipsoidFamily,
    PackingAlgo, GeneratorBackend, MicrostructureConfig,
)


# ── Fixtures ──────────────────────────────────────────────────────────────

@pytest.fixture
def small_box():
    return SimulationBox(Lx=1e-3, Ly=1e-3, Lz=1e-3, shape=GlobalShape.TORE)

@pytest.fixture
def sphere_vf():
    return SphereFamily(radius=5e-5, volume_fraction=0.15, phase=1)

@pytest.fixture
def sphere_n():
    return SphereFamily(radius=3e-5, n_spheres=20, phase=2, algo=PackingAlgo.WP)

@pytest.fixture
def basic_config(small_box, sphere_vf):
    return MicrostructureConfig(
        name="test_config",
        box=small_box,
        inclusions=[sphere_vf],
        voxel_grid=VoxelGrid(nx=64, ny=64, nz=64),
        seed=42,
    )


# ── Tests SimulationBox ────────────────────────────────────────────────────

class TestSimulationBox:
    def test_volume(self, small_box):
        assert math.isclose(small_box.volume, 1e-9)

    def test_as_merope_L(self, small_box):
        L = small_box.as_merope_L()
        assert L == [1e-3, 1e-3, 1e-3]
        assert len(L) == 3

    def test_negative_dimension_rejected(self):
        with pytest.raises(Exception):
            SimulationBox(Lx=-1e-3, Ly=1e-3, Lz=1e-3)

    def test_zero_dimension_rejected(self):
        with pytest.raises(Exception):
            SimulationBox(Lx=0, Ly=1e-3, Lz=1e-3)

    def test_shape_tore_default(self):
        box = SimulationBox()
        assert box.shape == GlobalShape.TORE


# ── Tests SphereFamily ─────────────────────────────────────────────────────

class TestSphereFamily:
    def test_volume_fraction_mode(self, sphere_vf, small_box):
        assert sphere_vf.volume_fraction == 0.15
        assert sphere_vf.n_spheres is None
        # n estimé > 0
        assert sphere_vf.estimated_n(small_box) > 0

    def test_n_spheres_mode(self, sphere_n, small_box):
        assert sphere_n.n_spheres == 20
        assert sphere_n.volume_fraction is None
        assert sphere_n.estimated_n(small_box) == 20

    def test_mutual_exclusion(self):
        with pytest.raises(Exception, match="mutuellement exclusifs"):
            SphereFamily(radius=1e-5, volume_fraction=0.1, n_spheres=10, phase=1)

    def test_no_sizing_rejected(self):
        with pytest.raises(Exception, match="volume_fraction OU n_spheres"):
            SphereFamily(radius=1e-5, phase=1)

    def test_inclusion_volume(self):
        s = SphereFamily(radius=1e-4, volume_fraction=0.1, phase=1)
        expected = (4/3) * math.pi * (1e-4)**3
        assert math.isclose(s.inclusion_volume(), expected, rel_tol=1e-9)

    def test_algo_wp(self, sphere_n):
        assert sphere_n.algo == PackingAlgo.WP


# ── Tests MultiSphereFamily ────────────────────────────────────────────────

class TestMultiSphereFamily:
    def test_valid_distribution(self):
        m = MultiSphereFamily(
            distribution=[(5e-5, 0.15, 1), (2e-5, 0.10, 2)]
        )
        desired, phases = m.to_merope_args()
        assert desired == [[5e-5, 0.15], [2e-5, 0.10]]
        assert phases == [1, 2]

    def test_total_vf_lt_1(self):
        with pytest.raises(Exception, match="Fraction volumique totale"):
            MultiSphereFamily(
                distribution=[(5e-5, 0.6, 1), (2e-5, 0.5, 2)]
            )

    def test_min_two_families(self):
        with pytest.raises(Exception):
            MultiSphereFamily(distribution=[(5e-5, 0.1, 1)])


# ── Tests VoxelGrid ────────────────────────────────────────────────────────

class TestVoxelGrid:
    def test_default_power2(self):
        g = VoxelGrid()
        assert g.nx == 64
        # 64 est une puissance de 2
        assert (g.nx & (g.nx - 1)) == 0

    def test_as_merope_list(self):
        g = VoxelGrid(nx=128, ny=64, nz=32)
        assert g.as_merope_list() == [128, 64, 32]

    def test_non_power2_warning(self):
        import warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            VoxelGrid(nx=100, ny=64, nz=64)
            assert any("puissance de 2" in str(warning.message) for warning in w)


# ── Tests MicrostructureConfig ─────────────────────────────────────────────

class TestMicrostructureConfig:
    def test_volume_fraction(self, basic_config, small_box, sphere_vf):
        vf = basic_config.volume_fraction_estimate()
        assert 0 < vf < 1.0

    def test_voxel_size_um(self, basic_config):
        vx, vy, vz = basic_config.voxel_size_um()
        # 1e-3 / 64 * 1e6 = 15.625 µm
        assert math.isclose(vx, 15.625, rel_tol=1e-4)

    def test_high_vf_rejected(self, small_box):
        # radius=5e-4 dans une boîte 1e-3 : V_sphere/V_box = (4/3*pi*(5e-4)^3)/(1e-9) ≈ 0.524
        # n_spheres=4 → vf ≈ 2.09 > 0.90 → rejeté
        with pytest.raises(Exception, match="Fraction volumique"):
            MicrostructureConfig(
                name="overflow",
                box=small_box,
                inclusions=[SphereFamily(radius=5e-4, n_spheres=4, phase=1)],
            )

    def test_json_roundtrip(self, basic_config):
        json_str = basic_config.to_json()
        restored = MicrostructureConfig.from_json(json_str)
        assert restored.name == basic_config.name
        assert restored.seed == basic_config.seed
        assert len(restored.inclusions) == len(basic_config.inclusions)

    def test_json_is_valid(self, basic_config):
        data = json.loads(basic_config.to_json())
        assert "name" in data
        assert "box" in data
        assert "inclusions" in data

    # ── Vérification de l'API Mérope dans le script généré ────────────────

    def test_merope_script_contains_setLength(self, basic_config):
        script = basic_config.to_merope_script()
        assert "setLength(L)" in script

    def test_merope_script_contains_setRadiusGenerator(self, basic_config):
        script = basic_config.to_merope_script()
        assert "setRadiusGenerator" in script

    def test_merope_script_contains_MultiInclusions(self, basic_config):
        script = basic_config.to_merope_script()
        assert "MultiInclusions_3D" in script

    def test_merope_script_contains_voxelation(self, basic_config):
        script = basic_config.to_merope_script()
        assert "Voxellation_3D" in script
        assert ".vtk" in script

    def test_merope_script_contains_seed(self, basic_config):
        script = basic_config.to_merope_script()
        assert "setSeed(42)" in script

    def test_merope_script_has_shebang(self, basic_config):
        script = basic_config.to_merope_script()
        assert script.startswith("#!/usr/bin/env python3")

    def test_combs_script_has_box(self, basic_config):
        script = basic_config.to_combs_script()
        assert "combs.Box(" in script

    def test_combs_script_has_generator(self, basic_config):
        script = basic_config.to_combs_script()
        assert "combs.Generator(" in script

    def test_combs_script_has_sphere_family(self, basic_config):
        script = basic_config.to_combs_script()
        assert "combs.SphereFamily(" in script

    # ── Multi-familles ────────────────────────────────────────────────────

    def test_two_sphere_families(self, small_box):
        cfg = MicrostructureConfig(
            name="two_families",
            box=small_box,
            inclusions=[
                SphereFamily(radius=5e-5, volume_fraction=0.10, phase=1),
                SphereFamily(radius=2e-5, volume_fraction=0.08, phase=2),
            ],
        )
        script = cfg.to_merope_script()
        assert "multi_0" in script
        assert "multi_1" in script

    def test_laguerre_script(self, small_box):
        # Laguerre remplit tout le volume → on désactive la validation vf via model_construct
        cfg = MicrostructureConfig.model_construct(
            name="polycristal",
            backend=GeneratorBackend.MEROPE,
            box=small_box,
            inclusions=[LaguerreTessFamily(n_grains=30)],
            voxel_grid=VoxelGrid(nx=64, ny=64, nz=64, output_vtk="polycristal.vtk"),
            seed=7,
        )
        script = cfg.to_merope_script()
        assert "LaguerreTess_3D" in script
        assert "30" in script

    # ── Sérialisation discriminée ─────────────────────────────────────────

    def test_discriminated_union_sphere(self, small_box):
        cfg = MicrostructureConfig(
            name="disc_test",
            box=small_box,
            inclusions=[SphereFamily(radius=1e-5, volume_fraction=0.05, phase=1)],
        )
        restored = MicrostructureConfig.from_json(cfg.to_json())
        assert isinstance(restored.inclusions[0], SphereFamily)

    def test_discriminated_union_laguerre(self, small_box):
        cfg = MicrostructureConfig.model_construct(
            name="laguerre_test",
            backend=GeneratorBackend.MEROPE,
            box=small_box,
            inclusions=[LaguerreTessFamily(n_grains=20)],
            voxel_grid=VoxelGrid(),
            seed=None,
        )
        restored = MicrostructureConfig.from_json(cfg.to_json())
        assert isinstance(restored.inclusions[0], LaguerreTessFamily)


# ── Tests de performance basique ──────────────────────────────────────────

class TestPerformance:
    def test_json_dump_speed(self, basic_config):
        import time
        t0 = time.perf_counter()
        for _ in range(1000):
            basic_config.to_json()
        elapsed = time.perf_counter() - t0
        # 1000 sérialisations < 1 seconde
        assert elapsed < 1.0, f"Sérialisation trop lente : {elapsed:.3f}s pour 1000 itérations"

    def test_codegen_speed(self, basic_config):
        import time
        t0 = time.perf_counter()
        for _ in range(500):
            basic_config.to_merope_script()
        elapsed = time.perf_counter() - t0
        assert elapsed < 1.0, f"Codegen trop lent : {elapsed:.3f}s pour 500 itérations"

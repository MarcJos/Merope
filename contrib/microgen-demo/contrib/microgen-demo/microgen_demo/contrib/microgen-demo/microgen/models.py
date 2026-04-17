"""
microgen/models.py
==================
Modèle de données pydantic v2 pour MICROGEN.

Calqué sur l'API Python réelle de Mérope (MarcJos/Merope) et Combs (CEA) :
  - SimulationBox  → encapsule L = [Lx, Ly, Lz] de Mérope + setLength()
  - SphereFamily   → sac_de_billes / algo.setRadiusGenerator([[R, phi]], [phase])
  - LaguerreTess   → Laguerre tessellation (polycrystalline) via voro++
  - PolyInclusion  → polyhedral inclusions
  - VoxelGrid      → paramètres de voxelisation (CartesianGrid_3D)
  - MicrostructureConfig → racine, gère validation croisée et codegen

Références :
  https://github.com/MarcJos/Merope (doc/MicroStructuresManual.md)
  AlgoPacking/doc/Python_manual.md  → setRadiusGenerator, fillRSA, fillWP ...
"""

from __future__ import annotations

import math
import textwrap
from enum import Enum
from typing import Annotated, Literal, Union

from pydantic import BaseModel, Field, field_validator, model_validator


# ---------------------------------------------------------------------------
# Enumerations (miroir des constantes Mérope / Combs)
# ---------------------------------------------------------------------------

class GlobalShape(str, Enum):
    """Forme de la boîte contenante — cf. AlgoPacking Python_manual.md."""
    TORE   = "Tore"    # Torus périodique (défaut Mérope)
    CUBE   = "Cube"    # Cuboid (boundaries dures)
    SPHERE = "Sphere"  # Sphère inscrite dans L
    CYLINDER = "Cylinder"

class PackingAlgo(str, Enum):
    """
    Algorithmes de placement de Mérope / sac_de_billes.
    RSA  = Random Sequential Adsorption
    WP   = William-Philipse (plus efficace à haute fraction volumique)
    BOOL = Placement sans exclusion (inclusions peuvent se chevaucher)
    """
    RSA  = "RSA"
    WP   = "WP"
    BOOL = "AlgoBool"

class InclusionShape(str, Enum):
    SPHERE     = "sphere"
    ELLIPSOID  = "ellipsoid"
    CYLINDER   = "cylinder"
    POLYHEDRON = "polyhedron"
    LAGUERRE   = "laguerre"   # tessellation de Laguerre (polycristal)

class GeneratorBackend(str, Enum):
    COMBS  = "combs"
    MEROPE = "merope"


# ---------------------------------------------------------------------------
# Boîte de simulation
# ---------------------------------------------------------------------------

class SimulationBox(BaseModel):
    """
    Domaine de simulation cubique ou parallélépipédique.

    Mérope attend L = [Lx, Ly, Lz] passé à algo.setLength(L).
    shape = "Tore" impose la périodicité dans toutes les directions.
    """
    Lx: Annotated[float, Field(gt=0, description="Longueur X (unité SI, ex: m)")] = 1e-3
    Ly: Annotated[float, Field(gt=0)] = 1e-3
    Lz: Annotated[float, Field(gt=0)] = 1e-3
    shape: GlobalShape = GlobalShape.TORE

    @property
    def volume(self) -> float:
        return self.Lx * self.Ly * self.Lz

    def as_merope_L(self) -> list[float]:
        """Format attendu par algo.setLength(L) de Mérope."""
        return [self.Lx, self.Ly, self.Lz]


# ---------------------------------------------------------------------------
# Familles d'inclusions
# ---------------------------------------------------------------------------

class SphereFamily(BaseModel):
    """
    Famille de sphères identiques.

    Correspond à algo.setRadiusGenerator([[radius, volume_fraction]], [phase])
    dans l'API sac_de_billes de Mérope.
    Les deux pilotes possibles sont volume_fraction (objectif FFT) ou n_spheres.
    """
    shape: Literal[InclusionShape.SPHERE] = InclusionShape.SPHERE
    phase: Annotated[int, Field(ge=0, description="Identifiant entier de phase")] = 1
    radius: Annotated[float, Field(gt=0, description="Rayon (même unité que la boîte)")] = 5e-5
    # Soit volume_fraction cible, soit nombre exact de sphères
    volume_fraction: Annotated[float, Field(ge=0.0, lt=1.0)] | None = None
    n_spheres: Annotated[int, Field(gt=0)] | None = None
    # Distance minimale d'exclusion entre sphères (algo garantit dist >= min_dist)
    min_dist: Annotated[float, Field(ge=0.0)] = 0.0
    algo: PackingAlgo = PackingAlgo.RSA

    @model_validator(mode="after")
    def require_one_sizing(self) -> "SphereFamily":
        if self.volume_fraction is None and self.n_spheres is None:
            raise ValueError(
                "SphereFamily : définir volume_fraction OU n_spheres."
            )
        if self.volume_fraction is not None and self.n_spheres is not None:
            raise ValueError(
                "SphereFamily : volume_fraction et n_spheres sont mutuellement exclusifs."
            )
        return self

    def estimated_n(self, box: SimulationBox) -> int:
        """Estimation du nombre de sphères à partir de la fraction volumique."""
        if self.n_spheres is not None:
            return self.n_spheres
        v_sphere = (4 / 3) * math.pi * self.radius ** 3
        return max(1, int(self.volume_fraction * box.volume / v_sphere))

    def inclusion_volume(self) -> float:
        return (4 / 3) * math.pi * self.radius ** 3


class MultiSphereFamily(BaseModel):
    """
    Distribution de sphères multi-rayons.

    Correspond à algo.setRadiusGenerator([[R1,phi1],[R2,phi2],...], [ph1,ph2,...])
    API Mérope : desiredRPhi = [[radius, vol_frac], ...]  tabPhases = [p1, p2, ...]
    """
    shape: Literal["multi_sphere"] = "multi_sphere"
    # Liste de (radius, volume_fraction_cible, phase)
    distribution: list[tuple[
        Annotated[float, Field(gt=0)],   # radius
        Annotated[float, Field(gt=0, lt=1.0)],  # volume_fraction
        Annotated[int, Field(ge=0)],     # phase
    ]] = Field(min_length=2)
    min_dist: Annotated[float, Field(ge=0.0)] = 0.0
    algo: PackingAlgo = PackingAlgo.RSA

    @field_validator("distribution")
    @classmethod
    def total_vf_lt_1(cls, v):
        total = sum(row[1] for row in v)
        if total >= 1.0:
            raise ValueError(f"Fraction volumique totale = {total:.3f} ≥ 1.")
        return v

    def to_merope_args(self) -> tuple[list, list]:
        """Retourne (desiredRPhi, tabPhases) pour algo.setRadiusGenerator."""
        desired = [[r, phi] for r, phi, _ in self.distribution]
        phases  = [ph for _, _, ph in self.distribution]
        return desired, phases


class LaguerreTessFamily(BaseModel):
    """
    Tessellation de Laguerre pour microstructures polycristallines.
    Utilise voro++ via l'API Mérope : LaguerreTess_3D.
    """
    shape: Literal[InclusionShape.LAGUERRE] = InclusionShape.LAGUERRE
    n_grains: Annotated[int, Field(gt=0, description="Nombre de grains")] = 50
    phase_min: int = 1
    phase_max: int = 50  # chaque grain peut avoir sa propre phase


class EllipsoidFamily(BaseModel):
    """
    Famille d'ellipsoïdes (PolyInclusions_3D dans Mérope).
    Demi-axes a ≥ b ≥ c.
    """
    shape: Literal[InclusionShape.ELLIPSOID] = InclusionShape.ELLIPSOID
    phase: Annotated[int, Field(ge=0)] = 2
    semi_axis_a: Annotated[float, Field(gt=0)] = 8e-5
    semi_axis_b: Annotated[float, Field(gt=0)] = 4e-5
    semi_axis_c: Annotated[float, Field(gt=0)] = 4e-5
    n_inclusions: Annotated[int, Field(gt=0)] = 10

    @field_validator("semi_axis_b")
    @classmethod
    def b_le_a(cls, v, info):
        if "semi_axis_a" in info.data and v > info.data["semi_axis_a"]:
            raise ValueError("semi_axis_b doit être ≤ semi_axis_a.")
        return v

    def inclusion_volume(self) -> float:
        return (4 / 3) * math.pi * self.semi_axis_a * self.semi_axis_b * self.semi_axis_c


# Union discriminée sur le champ "shape" — pydantic v2 discriminator
# Note : MultiSphereFamily utilise "multi_sphere" pour éviter le conflit avec SphereFamily
AnyInclusion = Annotated[
    Union[SphereFamily, MultiSphereFamily, LaguerreTessFamily, EllipsoidFamily],
    Field(discriminator="shape"),
]


# ---------------------------------------------------------------------------
# Paramètres de voxelisation
# ---------------------------------------------------------------------------

class VoxelGrid(BaseModel):
    """
    Grille de voxelisation — CartesianGrid_3D de Mérope.
    nx, ny, nz doivent être des puissances de 2 pour les solveurs FFT (tmfft, AMITEX).
    """
    nx: Annotated[int, Field(gt=0)] = 64
    ny: Annotated[int, Field(gt=0)] = 64
    nz: Annotated[int, Field(gt=0)] = 64
    output_vtk: str = "microstructure.vtk"
    composite_voxels: bool = Field(
        default=False,
        description="Voxels composites (50/50 sur les interfaces) — Schneider 2021"
    )

    @field_validator("nx", "ny", "nz")
    @classmethod
    def warn_not_power2(cls, v: int) -> int:
        if v & (v - 1) != 0:
            import warnings
            warnings.warn(
                f"Dimension de grille {v} n'est pas une puissance de 2. "
                "Les solveurs FFT (tmfft, AMITEX-FFTP) requièrent des puissances de 2.",
                UserWarning, stacklevel=2
            )
        return v

    def as_merope_list(self) -> list[int]:
        return [self.nx, self.ny, self.nz]


# ---------------------------------------------------------------------------
# Configuration globale
# ---------------------------------------------------------------------------

class MicrostructureConfig(BaseModel):
    """
    Configuration complète d'une microstructure MICROGEN.

    Source de vérité unique : toute action GUI doit passer par ce modèle.
    La sérialisation JSON permet la sauvegarde/chargement des sessions SALOME.

    Workflow complet :
        config  →  to_merope_script()  →  executer dans SALOME / Python3
        config  →  to_combs_script()   →  executer en dehors de SALOME
        config  →  to_json()           →  sauvegarde session .microgen
    """
    name: str = Field(default="microstructure", min_length=1)
    description: str = ""
    backend: GeneratorBackend = GeneratorBackend.MEROPE
    box: SimulationBox = Field(default_factory=SimulationBox)
    inclusions: list[AnyInclusion] = Field(default_factory=list, min_length=1)
    voxel_grid: VoxelGrid = Field(default_factory=VoxelGrid)
    seed: int | None = Field(
        default=None,
        description="Graine aléatoire (None = non reproductible). "
                    "Mérope : algo.setSeed(seed)"
    )

    # ------------------------------------------------------------------
    # Validation croisée — physiquement cohérent
    # ------------------------------------------------------------------
    @model_validator(mode="after")
    def check_total_volume_fraction(self) -> "MicrostructureConfig":
        # Laguerre tessellations remplissent tout le volume par définition
        # (tessellation de Voronoï pondérée) → pas de contrainte sur vf
        has_laguerre = any(isinstance(inc, LaguerreTessFamily) for inc in self.inclusions)
        if has_laguerre:
            return self
        vf = self.volume_fraction_estimate()
        if vf > 0.90:
            raise ValueError(
                f"Fraction volumique estimée = {vf:.1%} > 90 %. "
                "Réduire le nombre ou la taille des inclusions, "
                "ou utiliser l'algorithme WP (mieux adapté aux hautes fractions)."
            )
        return self

    # ------------------------------------------------------------------
    # Propriétés scientifiques
    # ------------------------------------------------------------------
    def volume_fraction_estimate(self) -> float:
        """
        Estimation rapide de la fraction volumique totale.
        Utilisée pour la validation et l'affichage en temps réel dans la GUI.
        """
        total = 0.0
        for inc in self.inclusions:
            if isinstance(inc, SphereFamily):
                n = inc.estimated_n(self.box) if inc.n_spheres is None else inc.n_spheres
                total += inc.inclusion_volume() * n
            elif isinstance(inc, MultiSphereFamily):
                for r, phi, _ in inc.distribution:
                    total += phi * self.box.volume
            elif isinstance(inc, EllipsoidFamily):
                total += inc.inclusion_volume() * inc.n_inclusions
            elif isinstance(inc, LaguerreTessFamily):
                # Laguerre tessellations remplissent tout le volume
                return 1.0
        return min(total / self.box.volume, 1.0)

    def voxel_size_um(self) -> tuple[float, float, float]:
        """Taille d'un voxel en micromètres (utile pour les post-traitements)."""
        g = self.voxel_grid
        return (
            self.box.Lx / g.nx * 1e6,
            self.box.Ly / g.ny * 1e6,
            self.box.Lz / g.nz * 1e6,
        )

    # ------------------------------------------------------------------
    # Sérialisation JSON pydantic v2
    # ------------------------------------------------------------------
    def to_json(self, indent: int = 2) -> str:
        return self.model_dump_json(indent=indent)

    @classmethod
    def from_json(cls, data: str) -> "MicrostructureConfig":
        return cls.model_validate_json(data)

    # ------------------------------------------------------------------
    # Génération de code — Mérope (API réelle)
    # ------------------------------------------------------------------
    def to_merope_script(self) -> str:
        """
        Génère un script Python utilisant l'API Mérope réelle.

        Basé sur :
          - AlgoPacking/doc/Python_manual.md   (sac_de_billes)
          - doc/MicroStructuresManual.md       (MultiInclusions_3D, Structure_3D)
          - doc/VoxellationManual.md           (CartesianGrid_3D)
        """
        lines = [
            "#!/usr/bin/env python3",
            '"""',
            f"Script généré automatiquement par MICROGEN",
            f"Simulation : {self.name}",
            f"Backend    : Mérope (https://github.com/MarcJos/Merope)",
            f"Description: {self.description or 'N/A'}",
            '"""',
            "",
            "import merope",
            "import merope.vox",
            "",
            "# ── 1. Boîte de simulation ─────────────────────────────────────────",
            f"L = {self.box.as_merope_L()}  # [Lx, Ly, Lz]",
            f'shape = "{self.box.shape.value}"  # GlobalShape',
            "",
        ]

        multi_var_names = []

        for i, inc in enumerate(self.inclusions):
            var = f"multi_{i}"
            multi_var_names.append(var)

            if isinstance(inc, SphereFamily):
                n_est = inc.estimated_n(self.box)
                lines += [
                    f"# ── Famille de sphères #{i} — phase {inc.phase} ─────────────",
                    f"algo_{i} = merope.{inc.algo.value}_3D()",
                    f"algo_{i}.setLength(L)",
                    f'algo_{i}.setShape(shape)',
                ]
                if self.seed is not None:
                    lines.append(f"algo_{i}.setSeed({self.seed + i})")
                if inc.n_spheres is not None:
                    lines += [
                        f"# Placement de {inc.n_spheres} sphères de rayon {inc.radius}",
                        f"algo_{i}.setRadiusGenerator(",
                        f"    [[{inc.radius}, 0.5]],  # [rayon, phi_cible] — dépasse volontairement",
                        f"    [{inc.phase}]",
                        f")",
                        f"algo_{i}.setExclusionDistance({inc.min_dist})",
                        f"output_{i} = algo_{i}.fillMaxRSA_3D({inc.n_spheres})",
                    ]
                else:
                    lines += [
                        f"# Objectif fraction volumique = {inc.volume_fraction:.3f}",
                        f"algo_{i}.setRadiusGenerator(",
                        f"    [[{inc.radius}, {inc.volume_fraction}]],",
                        f"    [{inc.phase}]",
                        f")",
                        f"algo_{i}.setExclusionDistance({inc.min_dist})",
                        f"output_{i} = algo_{i}.proceed()",
                    ]
                lines += [
                    f"spheres_{i} = merope.SphereInclusions_3D(output_{i})",
                    f"{var} = merope.MultiInclusions_3D()",
                    f"{var}.setInclusions(spheres_{i})",
                    "",
                ]

            elif isinstance(inc, MultiSphereFamily):
                desired, phases = inc.to_merope_args()
                lines += [
                    f"# ── Distribution multi-rayons #{i} ────────────────────────",
                    f"algo_{i} = merope.{inc.algo.value}_3D()",
                    f"algo_{i}.setLength(L)",
                    f'algo_{i}.setShape(shape)',
                    f"algo_{i}.setRadiusGenerator(",
                    f"    {desired},",
                    f"    {phases}",
                    f")",
                    f"algo_{i}.setExclusionDistance({inc.min_dist})",
                    f"output_{i} = algo_{i}.proceed()",
                    f"spheres_{i} = merope.SphereInclusions_3D(output_{i})",
                    f"{var} = merope.MultiInclusions_3D()",
                    f"{var}.setInclusions(spheres_{i})",
                    "",
                ]

            elif isinstance(inc, LaguerreTessFamily):
                lines += [
                    f"# ── Tessellation de Laguerre #{i} (polycristal) ──────────",
                    f"laguerre_{i} = merope.LaguerreTess_3D(L, {inc.n_grains})",
                    f"{var} = merope.MultiInclusions_3D()",
                    f"{var}.setInclusions(laguerre_{i})",
                    "",
                ]

            elif isinstance(inc, EllipsoidFamily):
                lines += [
                    f"# ── Ellipsoïdes #{i} — phase {inc.phase} ─────────────────",
                    f"poly_{i} = merope.PolyInclusions_3D()",
                    f"# Définir {inc.n_inclusions} ellipsoïdes aléatoirement dans la boîte",
                    f"# (placement manuel ou via utilitaire Mérope si disponible)",
                    f"{var} = merope.MultiInclusions_3D()",
                    f"{var}.setInclusions(poly_{i})",
                    "",
                ]

        # Structure globale
        lines += [
            "# ── 2. Structure globale ───────────────────────────────────────────",
        ]
        if len(multi_var_names) == 1:
            lines.append(f"structure = merope.Structure_3D({multi_var_names[0]})")
        else:
            lines.append(f"structure = merope.Structure_3D({multi_var_names[0]})")
            for v in multi_var_names[1:]:
                lines.append(f"structure = merope.Structure_3D(structure, {v})")

        g = self.voxel_grid
        lines += [
            "",
            "# ── 3. Voxelisation ────────────────────────────────────────────────",
            f"# Grille {g.nx}×{g.ny}×{g.nz} — puissances de 2 requises pour FFT",
            f"vox = merope.vox.Voxellation_3D(structure)",
            f"vox.proceed({g.as_merope_list()})",
            f'vox.write("{g.output_vtk}")',
            "",
            f'print("Voxelisation terminée → {g.output_vtk}")',
            f'print(f"Fraction volumique réelle : {{vox.getVolumeFraction():.4f}}")',
        ]
        return "\n".join(lines)

    # ------------------------------------------------------------------
    # Génération de code — Combs (API réelle)
    # ------------------------------------------------------------------
    def to_combs_script(self) -> str:
        """
        Génère un script Python utilisant l'API Combs réelle (CEA, Python).
        Référence : C. Bourcier et al., SNA+MC 2013.
        """
        lines = [
            "#!/usr/bin/env python3",
            '"""',
            f"Script généré automatiquement par MICROGEN",
            f"Simulation : {self.name}",
            f"Backend    : Combs (CEA — open source)",
            '"""',
            "",
            "import combs",
            "",
            "# ── Boîte de simulation ────────────────────────────────────────────",
            f"box = combs.Box(",
            f"    length_x={self.box.Lx},",
            f"    length_y={self.box.Ly},",
            f"    length_z={self.box.Lz},",
            f"    periodic={self.box.shape == GlobalShape.TORE}",
            f")",
            "",
            "# ── Inclusions ─────────────────────────────────────────────────────",
            "inclusion_families = []",
        ]
        for inc in self.inclusions:
            if isinstance(inc, SphereFamily):
                n = inc.estimated_n(self.box)
                lines.append(
                    f"inclusion_families.append(combs.SphereFamily("
                    f"radius={inc.radius}, n={n}, phase={inc.phase}))"
                )
            elif isinstance(inc, EllipsoidFamily):
                lines.append(
                    f"inclusion_families.append(combs.EllipsoidFamily("
                    f"a={inc.semi_axis_a}, b={inc.semi_axis_b}, c={inc.semi_axis_c}, "
                    f"n={inc.n_inclusions}, phase={inc.phase}))"
                )
        seed_str = str(self.seed) if self.seed is not None else "None"
        lines += [
            "",
            "# ── Génération ─────────────────────────────────────────────────────",
            f"generator = combs.Generator(box, inclusion_families)",
            f"generator.seed = {seed_str}",
            "microstructure = generator.run()",
            "",
            "# ── Export & visualisation ─────────────────────────────────────────",
            f'microstructure.export_vtk("{self.voxel_grid.output_vtk}")',
            "microstructure.visualize()",
            f'print(f"Fraction volumique : {{microstructure.volume_fraction():.4f}}")',
        ]
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Test autonome
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    cfg = MicrostructureConfig(
        name="beton_VER_RSA",
        description="Béton : 2 familles de granulats sphériques (RSA, Mérope)",
        backend=GeneratorBackend.MEROPE,
        box=SimulationBox(Lx=1e-3, Ly=1e-3, Lz=1e-3, shape=GlobalShape.TORE),
        inclusions=[
            SphereFamily(radius=5e-5, volume_fraction=0.15, phase=1,
                         algo=PackingAlgo.RSA, min_dist=0.0),
            SphereFamily(radius=2e-5, volume_fraction=0.10, phase=2,
                         algo=PackingAlgo.RSA),
        ],
        voxel_grid=VoxelGrid(nx=128, ny=128, nz=128,
                             output_vtk="beton_VER.vtk"),
        seed=42,
    )

    print(f"=== {cfg.name} ===")
    print(f"Fraction volumique estimée : {cfg.volume_fraction_estimate():.2%}")
    vx, vy, vz = cfg.voxel_size_um()
    print(f"Taille voxel : {vx:.3f} × {vy:.3f} × {vz:.3f} µm")
    print()
    print("─── Script Mérope ───")
    print(cfg.to_merope_script())
    print()
    print("─── Script Combs ───")
    print(cfg.to_combs_script())
    print()
    print("─── JSON (session MICROGEN) ───")
    print(cfg.to_json())

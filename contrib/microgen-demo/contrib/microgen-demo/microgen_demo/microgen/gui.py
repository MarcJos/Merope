"""
microgen/gui.py
===============
Interface graphique MICROGEN — PyQt5 / PySide6 compatible.

Architecture :
  ┌─────────────────────────────────────────────────────┐
  │  MicrogenMainWindow                                 │
  │  ├── ConfigPanel (gauche)                           │
  │  │   ├── BoxConfigWidget     → SimulationBox        │
  │  │   ├── InclusionListWidget → [AnyInclusion]       │
  │  │   ├── VoxelGridWidget     → VoxelGrid            │
  │  │   └── GenerationPanel     → backend, seed        │
  │  └── PreviewPanel (droite)                          │
  │      ├── Tab JSON            ← model.to_json()      │
  │      ├── Tab Script Mérope   ← model.to_merope_script()  │
  │      ├── Tab Script Combs    ← model.to_combs_script()   │
  │      └── StatusBar (vf, voxel size)                 │
  └─────────────────────────────────────────────────────┘

Le modèle pydantic est la SEULE source de vérité.
La GUI ne fait que construire/afficher le modèle.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Abstraction PyQt5 / PySide6 — objectif explicite du stage CEA
# ---------------------------------------------------------------------------
try:
    from PyQt5.QtWidgets import (
        QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
        QFormLayout, QGroupBox, QLabel, QLineEdit, QDoubleSpinBox,
        QSpinBox, QComboBox, QPushButton, QTextEdit, QTabWidget,
        QMessageBox, QCheckBox, QSplitter, QScrollArea, QFrame,
        QSizePolicy, QToolButton, QStackedWidget,
    )
    from PyQt5.QtCore import Qt, pyqtSignal as Signal, QTimer
    from PyQt5.QtGui import QFont, QColor, QPalette, QSyntaxHighlighter, QTextCharFormat
    QT_BACKEND = "PyQt5"
except ImportError:
    from PySide6.QtWidgets import (
        QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
        QFormLayout, QGroupBox, QLabel, QLineEdit, QDoubleSpinBox,
        QSpinBox, QComboBox, QPushButton, QTextEdit, QTabWidget,
        QMessageBox, QCheckBox, QSplitter, QScrollArea, QFrame,
        QSizePolicy, QToolButton, QStackedWidget,
    )
    from PySide6.QtCore import Qt, Signal, QTimer
    from PySide6.QtGui import QFont, QColor, QPalette, QSyntaxHighlighter, QTextCharFormat
    QT_BACKEND = "PySide6"

import sys
import re

from microgen.models import (
    MicrostructureConfig, SimulationBox, VoxelGrid,
    SphereFamily, MultiSphereFamily, EllipsoidFamily, LaguerreTessFamily,
    PackingAlgo, GeneratorBackend, GlobalShape, InclusionShape,
)


# ---------------------------------------------------------------------------
# Surbrillance Python minimale pour les onglets de script
# ---------------------------------------------------------------------------
class PythonHighlighter(QSyntaxHighlighter):
    KEYWORDS = r'\b(import|from|def|class|if|else|elif|for|while|return|True|False|None|and|or|not|in|is)\b'
    STRINGS   = r'(\"\"\".*?\"\"\"|\'\'\'.*?\'\'\'|\"[^\"]*\"|\'[^\']*\')'
    COMMENTS  = r'#[^\n]*'
    NUMBERS   = r'\b\d+\.?\d*([eE][+-]?\d+)?\b'

    def __init__(self, parent):
        super().__init__(parent)
        def fmt(color, bold=False):
            f = QTextCharFormat()
            f.setForeground(QColor(color))
            if bold:
                f.setFontWeight(700)
            return f
        self._rules = [
            (re.compile(self.KEYWORDS), fmt("#569cd6", bold=True)),
            (re.compile(self.STRINGS,  re.DOTALL), fmt("#ce9178")),
            (re.compile(self.COMMENTS), fmt("#6a9955")),
            (re.compile(self.NUMBERS), fmt("#b5cea8")),
        ]

    def highlightBlock(self, text):
        for pattern, fmt in self._rules:
            for m in pattern.finditer(text):
                self.setFormat(m.start(), m.end() - m.start(), fmt)


# ---------------------------------------------------------------------------
# Widget — Boîte de simulation
# ---------------------------------------------------------------------------
class BoxConfigWidget(QGroupBox):
    changed = Signal()

    def __init__(self, parent=None):
        super().__init__("📦  Boîte de simulation (VER)", parent)
        layout = QFormLayout()

        self.Lx = self._spin(1e-6, 1.0, 1e-3)
        self.Ly = self._spin(1e-6, 1.0, 1e-3)
        self.Lz = self._spin(1e-6, 1.0, 1e-3)

        self.shape = QComboBox()
        self.shape.addItems([s.value for s in GlobalShape])
        self.shape.setToolTip(
            "Tore = conditions périodiques (défaut Mérope)\n"
            "Cube = boundaries dures"
        )

        layout.addRow("Lx (m) :", self.Lx)
        layout.addRow("Ly (m) :", self.Ly)
        layout.addRow("Lz (m) :", self.Lz)
        layout.addRow("Forme :", self.shape)
        self.setLayout(layout)

        for w in (self.Lx, self.Ly, self.Lz):
            w.valueChanged.connect(self.changed)
        self.shape.currentIndexChanged.connect(self.changed)

    @staticmethod
    def _spin(lo, hi, default) -> QDoubleSpinBox:
        s = QDoubleSpinBox()
        s.setRange(lo, hi)
        s.setValue(default)
        s.setDecimals(8)
        s.setSingleStep(1e-4)
        s.setStepType(QDoubleSpinBox.AdaptiveDecimalStepType)
        return s

    def to_model(self) -> SimulationBox:
        return SimulationBox(
            Lx=self.Lx.value(), Ly=self.Ly.value(), Lz=self.Lz.value(),
            shape=GlobalShape(self.shape.currentText()),
        )


# ---------------------------------------------------------------------------
# Widget — Famille de sphères
# ---------------------------------------------------------------------------
class SphereFamilyWidget(QGroupBox):
    removed = Signal(object)
    changed = Signal()

    def __init__(self, index: int, parent=None):
        super().__init__(f"🔵  Sphères — famille #{index + 1}", parent)
        layout = QFormLayout()

        self.radius    = self._spin(1e-8, 1.0, 5e-5)
        self.phase     = QSpinBox(); self.phase.setRange(0, 999); self.phase.setValue(1)
        self.algo      = QComboBox(); self.algo.addItems([a.value for a in PackingAlgo])
        self.min_dist  = self._spin(0, 1.0, 0.0)
        self.min_dist.setToolTip("Distance minimale d'exclusion entre sphères (algo garantie)")

        # Sizing : volume_fraction OU n_spheres
        self.use_vf    = QCheckBox("Par fraction volumique")
        self.use_vf.setChecked(True)
        self.vf        = self._spin(0.001, 0.89, 0.15)
        self.n_spheres = QSpinBox(); self.n_spheres.setRange(1, 100000); self.n_spheres.setValue(50)
        self.n_spheres.setEnabled(False)

        self.use_vf.stateChanged.connect(self._toggle_sizing)
        self.use_vf.stateChanged.connect(self.changed)

        btn_del = QPushButton("✕")
        btn_del.setFixedWidth(28)
        btn_del.setToolTip("Supprimer cette famille")
        btn_del.clicked.connect(lambda: self.removed.emit(self))

        layout.addRow("Rayon (m) :", self.radius)
        layout.addRow("Phase :", self.phase)
        layout.addRow("Algorithme :", self.algo)
        layout.addRow("Excl. dist. :", self.min_dist)
        layout.addRow("", self.use_vf)
        layout.addRow("Frac. vol. :", self.vf)
        layout.addRow("Nb sphères :", self.n_spheres)
        layout.addRow("", btn_del)
        self.setLayout(layout)

        for w in (self.radius, self.vf, self.min_dist):
            w.valueChanged.connect(self.changed)
        for w in (self.phase, self.n_spheres):
            w.valueChanged.connect(self.changed)
        self.algo.currentIndexChanged.connect(self.changed)

    def _toggle_sizing(self, state):
        checked = bool(state)
        self.vf.setEnabled(checked)
        self.n_spheres.setEnabled(not checked)

    @staticmethod
    def _spin(lo, hi, default) -> QDoubleSpinBox:
        s = QDoubleSpinBox()
        s.setRange(lo, hi)
        s.setValue(default)
        s.setDecimals(8)
        s.setStepType(QDoubleSpinBox.AdaptiveDecimalStepType)
        return s

    def to_model(self) -> SphereFamily:
        kwargs = dict(
            radius=self.radius.value(),
            phase=self.phase.value(),
            algo=PackingAlgo(self.algo.currentText()),
            min_dist=self.min_dist.value(),
        )
        if self.use_vf.isChecked():
            kwargs["volume_fraction"] = self.vf.value()
        else:
            kwargs["n_spheres"] = self.n_spheres.value()
        return SphereFamily(**kwargs)


# ---------------------------------------------------------------------------
# Widget — Tessellation de Laguerre
# ---------------------------------------------------------------------------
class LaguerreWidget(QGroupBox):
    removed = Signal(object)
    changed = Signal()

    def __init__(self, index: int, parent=None):
        super().__init__(f"🟩  Laguerre Tess. — famille #{index + 1}", parent)
        layout = QFormLayout()
        self.n_grains = QSpinBox()
        self.n_grains.setRange(2, 10000)
        self.n_grains.setValue(50)
        self.n_grains.setToolTip("Nombre de grains du polycristal (voro++)")
        btn_del = QPushButton("✕"); btn_del.setFixedWidth(28)
        btn_del.clicked.connect(lambda: self.removed.emit(self))
        layout.addRow("Nb grains :", self.n_grains)
        layout.addRow("", btn_del)
        self.setLayout(layout)
        self.n_grains.valueChanged.connect(self.changed)

    def to_model(self) -> LaguerreTessFamily:
        return LaguerreTessFamily(n_grains=self.n_grains.value())


# ---------------------------------------------------------------------------
# Widget — Grille de voxelisation
# ---------------------------------------------------------------------------
class VoxelGridWidget(QGroupBox):
    changed = Signal()

    def __init__(self, parent=None):
        super().__init__("🔲  Voxelisation (CartesianGrid_3D)", parent)
        layout = QFormLayout()

        self.nx = self._spin2(16, 512, 64)
        self.ny = self._spin2(16, 512, 64)
        self.nz = self._spin2(16, 512, 64)
        self.output = QLineEdit("microstructure.vtk")
        self.composite = QCheckBox("Voxels composites (interfaces 50/50)")
        self.composite.setToolTip("Schneider 2021 — mixage de phase sur les bords")

        layout.addRow("nx :", self.nx)
        layout.addRow("ny :", self.ny)
        layout.addRow("nz :", self.nz)
        layout.addRow("Fichier VTK :", self.output)
        layout.addRow("", self.composite)
        self.setLayout(layout)

        for w in (self.nx, self.ny, self.nz):
            w.valueChanged.connect(self.changed)

    @staticmethod
    def _spin2(lo, hi, default) -> QSpinBox:
        s = QSpinBox(); s.setRange(lo, hi); s.setValue(default)
        s.setToolTip("Puissance de 2 recommandée pour solveurs FFT (tmfft, AMITEX)")
        return s

    def to_model(self) -> VoxelGrid:
        return VoxelGrid(
            nx=self.nx.value(), ny=self.ny.value(), nz=self.nz.value(),
            output_vtk=self.output.text() or "microstructure.vtk",
            composite_voxels=self.composite.isChecked(),
        )


# ---------------------------------------------------------------------------
# Fenêtre principale
# ---------------------------------------------------------------------------
class MicrogenMainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"MICROGEN — pydantic v2 + Mérope  [{QT_BACKEND}]")
        self.setMinimumSize(1100, 700)
        self._inclusion_widgets: list = []

        # Timer pour debounce (éviter de régénérer à chaque frappe)
        self._debounce = QTimer()
        self._debounce.setSingleShot(True)
        self._debounce.setInterval(300)
        self._debounce.timeout.connect(self._update_preview)

        # ── Layout principal ──────────────────────────────────────────────
        splitter = QSplitter(Qt.Horizontal)

        # ── Panneau gauche ────────────────────────────────────────────────
        left_scroll = QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_scroll.setMinimumWidth(380)
        left_scroll.setMaximumWidth(480)

        left_inner = QWidget()
        left_layout = QVBoxLayout(left_inner)
        left_layout.setSpacing(6)

        # Nom & backend
        meta_box = QGroupBox("ℹ️  Métadonnées")
        meta_form = QFormLayout()
        self.name_edit = QLineEdit("beton_VER_demo")
        self.backend   = QComboBox()
        self.backend.addItems([b.value for b in GeneratorBackend])
        self.seed_spin = QSpinBox()
        self.seed_spin.setRange(-1, 999999)
        self.seed_spin.setValue(42)
        self.seed_spin.setSpecialValueText("aléatoire")
        self.seed_spin.setMinimum(-1)
        meta_form.addRow("Nom :", self.name_edit)
        meta_form.addRow("Backend :", self.backend)
        meta_form.addRow("Seed :", self.seed_spin)
        meta_box.setLayout(meta_form)
        left_layout.addWidget(meta_box)

        # Boîte
        self.box_widget = BoxConfigWidget()
        self.box_widget.changed.connect(self._debounce.start)
        left_layout.addWidget(self.box_widget)

        # Inclusions
        inc_header = QGroupBox("🔶  Inclusions")
        inc_header_layout = QVBoxLayout()
        btn_add_sphere   = QPushButton("+ Sphères (RSA/WP)")
        btn_add_laguerre = QPushButton("+ Laguerre (polycristal)")
        btn_add_sphere.clicked.connect(self._add_sphere)
        btn_add_laguerre.clicked.connect(self._add_laguerre)
        hbtn = QHBoxLayout()
        hbtn.addWidget(btn_add_sphere)
        hbtn.addWidget(btn_add_laguerre)
        inc_header_layout.addLayout(hbtn)
        self._inc_container_layout = QVBoxLayout()
        inc_header_layout.addLayout(self._inc_container_layout)
        inc_header.setLayout(inc_header_layout)
        left_layout.addWidget(inc_header)

        # Voxelisation
        self.vox_widget = VoxelGridWidget()
        self.vox_widget.changed.connect(self._debounce.start)
        left_layout.addWidget(self.vox_widget)

        # Bouton générer
        btn_generate = QPushButton("🚀  Valider & Générer la configuration")
        btn_generate.setFixedHeight(38)
        btn_generate.setStyleSheet(
            "QPushButton { background:#005a9e; color:white; border-radius:4px; font-weight:bold; }"
            "QPushButton:hover { background:#0078d4; }"
        )
        btn_generate.clicked.connect(self._generate)
        left_layout.addWidget(btn_generate)
        left_layout.addStretch()

        left_scroll.setWidget(left_inner)

        # ── Panneau droit ─────────────────────────────────────────────────
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)

        self.status_label = QLabel("")
        self.status_label.setAlignment(Qt.AlignCenter)
        self.status_label.setFixedHeight(28)
        right_layout.addWidget(self.status_label)

        tabs = QTabWidget()
        mono = QFont("Courier New", 9)

        self.json_output   = self._make_editor(mono)
        self.merope_output = self._make_editor(mono)
        self.combs_output  = self._make_editor(mono)

        PythonHighlighter(self.merope_output.document())
        PythonHighlighter(self.combs_output.document())

        tabs.addTab(self.json_output,   "JSON  (session MICROGEN)")
        tabs.addTab(self.merope_output, "Script  Mérope")
        tabs.addTab(self.combs_output,  "Script  Combs")
        right_layout.addWidget(tabs)

        # ── Assemblage ────────────────────────────────────────────────────
        splitter.addWidget(left_scroll)
        splitter.addWidget(right_widget)
        splitter.setSizes([420, 680])
        self.setCentralWidget(splitter)

        # Connexions métadonnées
        self.name_edit.textChanged.connect(self._debounce.start)
        self.backend.currentIndexChanged.connect(self._debounce.start)
        self.seed_spin.valueChanged.connect(self._debounce.start)

        # Famille de sphères par défaut
        self._add_sphere()
        self._update_preview()

    # ------------------------------------------------------------------
    @staticmethod
    def _make_editor(font) -> QTextEdit:
        e = QTextEdit()
        e.setReadOnly(True)
        e.setFont(font)
        e.setStyleSheet("background:#1e1e1e; color:#d4d4d4;")
        return e

    # ------------------------------------------------------------------
    def _add_sphere(self):
        w = SphereFamilyWidget(len(self._inclusion_widgets))
        self._register_inclusion(w)

    def _add_laguerre(self):
        w = LaguerreWidget(len(self._inclusion_widgets))
        self._register_inclusion(w)

    def _register_inclusion(self, w):
        w.removed.connect(self._remove_inclusion)
        w.changed.connect(self._debounce.start)
        self._inclusion_widgets.append(w)
        self._inc_container_layout.addWidget(w)
        self._debounce.start()

    def _remove_inclusion(self, widget):
        self._inclusion_widgets.remove(widget)
        self._inc_container_layout.removeWidget(widget)
        widget.deleteLater()
        self._debounce.start()

    # ------------------------------------------------------------------
    def _build_config(self) -> MicrostructureConfig | None:
        try:
            inclusions = []
            for w in self._inclusion_widgets:
                inclusions.append(w.to_model())
            seed = self.seed_spin.value() if self.seed_spin.value() >= 0 else None
            return MicrostructureConfig(
                name=self.name_edit.text() or "microstructure",
                backend=GeneratorBackend(self.backend.currentText()),
                box=self.box_widget.to_model(),
                inclusions=inclusions,
                voxel_grid=self.vox_widget.to_model(),
                seed=seed,
            )
        except Exception as e:
            self.status_label.setText(f'<span style="color:#f44;">⚠ {e}</span>')
            return None

    def _update_preview(self):
        cfg = self._build_config()
        if cfg is None:
            return
        self.json_output.setPlainText(cfg.to_json())
        self.merope_output.setPlainText(cfg.to_merope_script())
        self.combs_output.setPlainText(cfg.to_combs_script())
        vf = cfg.volume_fraction_estimate()
        vx, vy, vz = cfg.voxel_size_um()
        color = "lime" if vf < 0.5 else "orange" if vf < 0.75 else "#f44"
        self.status_label.setText(
            f'<span style="color:{color}">Frac. vol. ≈ {vf:.1%}</span>'
            f'  |  Voxel : {vx:.2f} × {vy:.2f} × {vz:.2f} µm'
            f'  |  Backend : {cfg.backend.value}'
        )

    def _generate(self):
        cfg = self._build_config()
        if cfg is None:
            return
        try:
            validated = MicrostructureConfig.model_validate(cfg.model_dump())
            self._update_preview()
            vf = validated.volume_fraction_estimate()
            QMessageBox.information(
                self, "✔ Configuration validée",
                f"Simulation : {validated.name}\n"
                f"Backend    : {validated.backend.value}\n"
                f"Frac. vol. : {vf:.2%}\n"
                f"Grille     : {validated.voxel_grid.nx}³\n\n"
                "Les scripts Mérope et Combs sont disponibles dans les onglets →"
            )
        except Exception as e:
            QMessageBox.critical(self, "Erreur pydantic", str(e))


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import os, sys
    # Ajouter le répertoire parent au PYTHONPATH pour l'import relatif
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    win = MicrogenMainWindow()
    win.show()
    sys.exit(app.exec() if QT_BACKEND == "PySide6" else app.exec_())

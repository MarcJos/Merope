PYTHON     ?= python3
PIP        ?= pip3
PYTEST     ?= pytest
SRC        := microgen
TESTS      := tests
SCRIPTS    := scripts

.PHONY: all install install-qt5 install-qt6 install-dev \
        test test-fast coverage demo gui example \
        lint typecheck format clean help

all: install test

install:
	@echo "→ Installation des dépendances..."
	$(PIP) install "pydantic>=2.0" pytest pytest-cov -q
	@echo "  Pour la GUI : make install-qt5  ou  make install-qt6"

install-qt5:
	$(PIP) install PyQt5

install-qt6:
	$(PIP) install PySide6

install-dev:
	$(PIP) install "pydantic>=2.0" pytest pytest-cov ruff mypy

test:
	@echo "→ Tests unitaires (37 tests)..."
	PYTHONPATH=. $(PYTEST) $(TESTS)/ -v --tb=short

test-fast:
	PYTHONPATH=. $(PYTEST) $(TESTS)/ -v --tb=short -x

coverage:
	PYTHONPATH=. $(PYTEST) $(TESTS)/ --cov=$(SRC) --cov-report=html --cov-report=term-missing
	@echo "   HTML → htmlcov/index.html"

demo:
	PYTHONPATH=. $(PYTHON) $(SRC)/models.py

gui:
	PYTHONPATH=. $(PYTHON) $(SRC)/gui.py

example:
	@echo "→ Génération des scripts Mérope/Combs..."
	PYTHONPATH=. $(PYTHON) $(SCRIPTS)/generate_examples.py

lint:
	@command -v ruff >/dev/null 2>&1 && ruff check $(SRC)/ $(TESTS)/ || echo "ruff non installé"

typecheck:
	@command -v mypy >/dev/null 2>&1 && PYTHONPATH=. mypy $(SRC)/models.py --ignore-missing-imports || echo "mypy non installé"

format:
	@command -v ruff >/dev/null 2>&1 && ruff format $(SRC)/ $(TESTS)/ || echo "ruff non installé"

clean:
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete 2>/dev/null || true
	rm -rf .pytest_cache htmlcov .coverage .mypy_cache .ruff_cache
	rm -f $(SCRIPTS)/merope_*.py $(SCRIPTS)/combs_*.py $(SCRIPTS)/*.microgen
	@echo "Nettoyage terminé."

help:
	@echo ""
	@echo "  microgen-demo — Contribution stage CEA Saclay"
	@echo "  make install     pydantic + pytest"
	@echo "  make install-qt5 PyQt5 (GUI)"
	@echo "  make install-qt6 PySide6 (GUI)"
	@echo "  make test        37 tests unitaires"
	@echo "  make coverage    Rapport HTML"
	@echo "  make demo        CLI (JSON + scripts)"
	@echo "  make gui         Interface graphique"
	@echo "  make example     Scripts béton + polycristal"
	@echo "  make lint        ruff"
	@echo "  make typecheck   mypy"
	@echo "  make clean       supprime artefacts"
	@echo ""
	
	
example-laguerre:
	@echo "→ Génération du script Mérope (polycristal Laguerre)..."
	PYTHONPATH=. $(PYTHON) $(SCRIPTS)/generate_examples.py --laguerre

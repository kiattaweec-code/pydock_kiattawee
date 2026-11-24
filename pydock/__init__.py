"""pydock -- a flexible-ligand molecular docking engine.

A from-scratch, GOLD-shaped docking program: genetic-algorithm conformational
search over a torsion tree, rotatable receptor hydroxyls, and two
interchangeable scoring functions (empirical ChemScore/PLP-style, and an
AutoDock-style force field).

Quick start
-----------
>>> from pydock import DockingConfig, run_docking
>>> cfg = DockingConfig(receptor="protein.pdb", ligands=["ligand.sdf"],
...                     reference_ligand="crystal.sdf", scoring="empirical")
>>> results = run_docking(cfg)
>>> results[0].solutions[0].score
"""
from .dock import DockingConfig, LigandResult, run_docking, write_poses
from .ga import Docker, GAConfig, Solution, cluster_solutions
from .ligand import Ligand, load_ligand, prepare_ligand
from .receptor import Receptor, prepare_receptor, site_from_ligand, site_from_residues
from .scoring import Scorer, build_scorer

__version__ = "0.1.1"
__all__ = [
    "DockingConfig", "LigandResult", "run_docking", "write_poses",
    "Docker", "GAConfig", "Solution", "cluster_solutions",
    "Ligand", "load_ligand", "prepare_ligand",
    "Receptor", "prepare_receptor", "site_from_ligand", "site_from_residues",
    "Scorer", "build_scorer",
]

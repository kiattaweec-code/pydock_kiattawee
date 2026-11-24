"""Command-line interface.

Two ways in: flags, or a GOLD-style configuration file.

    pydock -r protein.pdb -l ligand.sdf --ref crystal.sdf
    pydock -c dock.conf

A configuration file is `key = value` lines, `#` for comments, and `ligand`
may be repeated:

    receptor         = protein.pdb
    ligand           = actives.sdf
    reference_ligand = crystal.sdf
    scoring          = empirical
    runs             = 10
    poses            = 5
    output           = results
"""
from __future__ import annotations

import argparse
import glob
import sys

from .dock import DockingConfig, run_docking
from .ga import GAConfig

_ALIASES = {
    "protein": "receptor", "rec": "receptor",
    "ligands": "ligand", "lig": "ligand",
    "ref": "reference_ligand", "reference": "reference_ligand",
    "output": "output_dir", "out": "output_dir",
    "scoring_function": "scoring", "fitness": "scoring",
    "runs": "n_runs", "ga_runs": "n_runs",
    "poses": "n_poses", "n_solutions": "n_poses",
    "population": "pop_size", "popsize": "pop_size",
    "islands": "n_islands", "niches": "n_islands",
    "evals": "n_evals", "evaluations": "n_evals",
    "grid_spacing": "spacing",
}

_GA_KEYS = {"pop_size", "n_islands", "n_evals", "niche_size", "niche_rmsd",
            "p_crossover", "p_mutate", "migrate_every", "tournament",
            "ls_frequency", "ls_steps", "ls_rho", "local_opt",
            "local_opt_maxiter", "early_stop", "seed", "autoscale",
            "use_cavity", "cavity_jump"}

_FLOATS = {"spacing", "padding", "radius", "cluster_rmsd", "niche_rmsd",
           "p_crossover", "p_mutate", "ls_frequency", "ls_rho", "cavity_jump",
           "cavity_spacing", "cavity_probe", "rot_h_cutoff"}
_INTS = {"n_runs", "n_poses", "pop_size", "n_islands", "n_evals", "niche_size",
         "migrate_every", "tournament", "ls_steps", "local_opt_maxiter",
         "early_stop", "seed", "cavity_min_psp"}
_BOOLS = {"keep_waters", "rotate_receptor_oh", "rotate_ligand_polar_h",
          "local_opt", "autoscale", "use_cavity", "verbose", "detect_cavity"}


def _coerce(key: str, value: str):
    if key in _BOOLS:
        return value.strip().lower() in ("1", "true", "yes", "on")
    if key in _INTS:
        return int(float(value))
    if key in _FLOATS:
        return float(value)
    return value.strip()


def read_config(path: str) -> dict:
    out: dict = {"ligand": []}
    with open(path) as fh:
        for raw in fh:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            if "=" not in line:
                parts = line.split(None, 1)
                if len(parts) != 2:
                    continue
                key, value = parts
            else:
                key, value = line.split("=", 1)
            key = _ALIASES.get(key.strip().lower(), key.strip().lower())
            if key == "ligand":
                out["ligand"].append(value.strip())
            elif key == "centre" or key == "center":
                out["centre"] = tuple(float(x) for x in
                                      value.replace(",", " ").split())
            elif key == "site_residues":
                res = []
                for tok in value.replace(",", " ").split():
                    chain, num = (tok.split(":") if ":" in tok
                                  else (tok[0], tok[1:]))
                    res.append((chain, int(num)))
                out["site_residues"] = res
            else:
                out[key] = _coerce(key, value)
    return out


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pydock",
        description="Flexible-ligand docking with a genetic algorithm.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    p.add_argument("-c", "--config", help="configuration file")
    p.add_argument("-r", "--receptor", help="receptor PDB")
    p.add_argument("-l", "--ligand", action="append", default=[],
                   help="ligand file; repeatable, globs allowed")
    p.add_argument("--ref", "--reference", dest="reference_ligand",
                   help="reference ligand, defines the site and gives an RMSD")
    p.add_argument("--centre", "--center", nargs=3, type=float,
                   metavar=("X", "Y", "Z"), help="site centre")
    p.add_argument("--radius", type=float, help="site radius with --centre")
    p.add_argument("--site-residues", nargs="+", metavar="CHAIN:RESNUM",
                   help="define the site from these residues")
    p.add_argument("-s", "--scoring", default="empirical",
                   choices=["empirical", "chemscore", "chemplp", "plp",
                            "forcefield", "ff", "autodock"],
                   help="scoring function (default: empirical)")
    p.add_argument("-n", "--runs", dest="n_runs", type=int, default=10,
                   help="independent GA runs per ligand (default: 10)")
    p.add_argument("-p", "--poses", dest="n_poses", type=int, default=5,
                   help="poses written per ligand (default: 5)")
    p.add_argument("-o", "--output", dest="output_dir", default="pydock_out")
    p.add_argument("--spacing", type=float, default=0.4, help="grid spacing")
    p.add_argument("--padding", type=float, default=4.0,
                   help="box padding around the reference ligand")
    p.add_argument("--evals", dest="n_evals", type=int,
                   help="score evaluations per run (disables autoscaling)")
    p.add_argument("--population", dest="pop_size", type=int, default=100)
    p.add_argument("--islands", dest="n_islands", type=int, default=5)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--keep-waters", action="store_true")
    p.add_argument("--no-cavity", action="store_true",
                   help="search the whole box instead of the detected cavity")
    p.add_argument("--no-receptor-oh", action="store_true",
                   help="keep receptor hydroxyls fixed")
    p.add_argument("-q", "--quiet", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    fields: dict = {}
    if args.config:
        fields = read_config(args.config)
        fields["ligands"] = fields.pop("ligand", [])

    if args.receptor:
        fields["receptor"] = args.receptor
    if args.ligand:
        fields["ligands"] = args.ligand
    if args.reference_ligand:
        fields["reference_ligand"] = args.reference_ligand
    if args.centre:
        fields["centre"] = tuple(args.centre)
    if args.radius is not None:
        fields["radius"] = args.radius
    if args.site_residues:
        fields["site_residues"] = [
            (t.split(":")[0], int(t.split(":")[1])) for t in args.site_residues]

    defaults = build_parser().parse_args([])
    for key in ("scoring", "n_runs", "n_poses", "output_dir", "spacing",
                "padding"):
        if getattr(args, key) != getattr(defaults, key) or key not in fields:
            fields[key] = getattr(args, key)
    if args.keep_waters:
        fields["keep_waters"] = True
    if args.no_cavity:
        fields["detect_cavity"] = False
    if args.no_receptor_oh:
        fields["rotate_receptor_oh"] = False
    fields["verbose"] = not args.quiet

    if not fields.get("receptor"):
        print("error: no receptor given (-r or `receptor =` in a config file)",
              file=sys.stderr)
        return 2
    expanded: list[str] = []
    for pattern in fields.get("ligands", []):
        hits = sorted(glob.glob(pattern))
        expanded.extend(hits or [pattern])
    fields["ligands"] = expanded
    if not expanded:
        print("error: no ligands given (-l or `ligand =` in a config file)",
              file=sys.stderr)
        return 2

    ga_fields = {k: fields.pop(k) for k in list(fields) if k in _GA_KEYS}
    for key in ("pop_size", "n_islands", "seed"):
        ga_fields.setdefault(key, getattr(args, key))
    if args.n_evals:
        ga_fields["n_evals"] = args.n_evals
    if "n_evals" in ga_fields:
        ga_fields["autoscale"] = False
        ga_fields.setdefault("early_stop", ga_fields["n_evals"])
    fields["ga"] = GAConfig(**ga_fields)

    cfg = DockingConfig(**fields)
    results = run_docking(cfg)

    if cfg.verbose and len(results) > 1:
        print("\n" + "=" * 64)
        print(f"{'ligand':<28} {'score':>10} {'modes':>6} {'rmsd':>7}")
        for r in sorted(results, key=lambda x: x.solutions[0].score):
            rms = "-" if r.rmsd_to_reference is None else f"{r.rmsd_to_reference:.2f}"
            print(f"{r.name[:28]:<28} {r.solutions[0].score:10.3f} "
                  f"{len(r.clusters):6d} {rms:>7}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

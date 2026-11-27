# pydock

A flexible-ligand molecular docking engine written from scratch — GOLD-shaped,
but no GOLD code, no CSD data, and no licence restrictions.

Genetic-algorithm conformational search over a torsion tree, cavity detection,
rotatable receptor hydroxyls as search variables, and two interchangeable
scoring functions.

```bash
pip install -e .
pydock -r protein.pdb -l ligand.sdf --ref crystal.sdf -s empirical -n 10
```

### As a desktop app

Download an installer for Windows or macOS from the releases page, or build
one yourself — `packaging/README.md` has the details. It opens in a window,
runs entirely on your machine, and needs no Python installed.

```
pydock                        # open the window
pydock --no-open --port 8000  # headless, for a shared workstation
pydock --selftest             # diagnose an installation
```

### From source

Or run it in a browser — upload, watch the search, inspect poses in 3D:

```bash
pip install fastapi "uvicorn[standard]" python-multipart
uvicorn pydock_web.app:app --port 8000      # then open http://127.0.0.1:8000
```

See [`pydock_web/README.md`](pydock_web/README.md) for the HTTP API.

## What it does

| Capability | Status |
|---|---|
| Flexible ligand (torsion tree, any number of rotors) | yes |
| Rigid-body placement (translation + orientation) | yes |
| Rotatable receptor OH / SH / NH₃⁺ as genes | yes, GOLD-style |
| Automatic polar-hydrogen placement on the receptor | yes |
| LIGSITE-style cavity detection | yes |
| Island-model GA with niching and migration | yes |
| Lamarckian local search (Solis–Wets) | yes |
| Two scoring functions, switchable at runtime | yes |
| Precomputed interaction grids | yes |
| Multiple binding modes with pose clustering | yes |
| Multi-record library files (SDF / SMILES / MOL2) | yes |
| Whole-protein pocket detection and ranking | yes |
| Arbitrary box (centre + x/y/z size), not just a cube | yes |
| Screening ranking with ligand efficiency, as CSV | yes |
| Symmetry-corrected RMSD | yes |
| Metal coordination term | yes |
| Ranked SDF output with per-term breakdown | yes |
| Structure preparation in the browser (split a complex) | yes |
| Ligand extraction to MOL2 with Sybyl atom types | yes |

## Preparing a structure

A PDB entry is not a receptor. The bound ligand is still in it, and the engine
keeps every hetero atom it is given, so a ligand left in the file fills the
pocket it was meant to reveal. Crystallographic hydrogens cause a quieter
problem: the engine places polar hydrogens itself and only treats hydroxyls as
rotatable when the file arrives without any, so a protonated receptor silently
loses that part of the search.

The **Prepare** page does the splitting. Load a complex, click the ligand, and
it writes the protein without it alongside the ligand as MOL2 and SDF, then
hands both to the docking page without a download. It also finds the waters,
buffers, cryoprotectants and metals, and removes them by default -- except
metals, which are usually part of the site.

One thing is worth understanding. A PDB file records no bond orders and no
protonation, so they have to come from somewhere:

- **Paste the ligand's SMILES** (from RCSB or PubChem) and the bond orders are
  transferred onto the crystal coordinates. The result is exactly the compound
  you name, with its formal charges, and the coordinates are untouched. The
  SMILES is checked against the residue's heavy-atom count, so the wrong
  compound is refused rather than quietly accepted.
- **Leave it empty** and they are worked out from the geometry. This succeeds
  on small, mostly-saturated ligands and fails on drug-sized ones: for a
  43-heavy-atom ergoline no net charge gives a valid structure. When it fails
  the file is still written, with every bond single and a warning saying so --
  good for placing a box or measuring an RMSD, not for docking.

RDKit cannot write MOL2, so `pydock/mol2.py` does, with Sybyl atom and bond
types (`C.ar`, `N.am`, `O.co2`, `C.cat`, `S.o2`, `am` bonds, and so on).
`packaging/test_mol2.py` writes 17 molecules chosen to exercise the awkward
cases, reads each one back with RDKit and compares.

## Scoring functions

Selected with `-s`. Lower is better in both cases.

**`empirical`** — ChemScore-style free-energy terms (Eldridge et al., *JCAMD*
1997): directional hydrogen bonds, metal coordination, lipophilic contact
surface, and a rotor entropy penalty, plus polarity-resolved steric repulsion.
This is the closer analogue of GOLD's ChemScore. `-s chemplp` adds a
piecewise-linear shape term on top.

**`forcefield`** — AutoDock4-style (Huey et al., *JCC* 2007): Lennard-Jones
12-6 van der Waals, screened electrostatics with the Mehler–Solmajer
distance-dependent dielectric, volume-based desolvation, and a directional
12-10 hydrogen bond.

Both evaluate hydrogen bonds **pairwise rather than from a grid**, which is
what keeps full donor *and* acceptor directionality and lets receptor hydroxyls
reorient during the search. Everything type-dependent and non-directional is
gridded and read by trilinear interpolation.

## Validation

Redocking ergotamine into the 5-HT receptor site it came from — 43 heavy
atoms, 7 heavy-atom rotors, 20 rotatable receptor hydroxyls in range, 35 search
variables in total. The starting conformer is generated from scratch by ETKDG,
so the search gets no information from the crystal geometry.

| Scoring function | Top-ranked RMSD | Best pose found | Runs on the right mode |
|---|---|---|---|
| `empirical` | **1.18 Å** | 1.16 Å | 5 / 8 |
| `forcefield` | 2.29 Å | 1.87 Å | 4 / 8 |

The parameterisation floor is **0.64 Å**: the rigid fragments keep their ETKDG
internal geometry, so no amount of searching could do better with this input
conformer. 1.18 Å is close to that floor.

Two details worth noting. The empirical function scores its own best pose at
−15.16 against −14.31 for the crystal pose, i.e. it is optimising a real
minimum rather than drifting. And every run that landed in the wrong binding
mode (9.6–10.0 Å) scored worse than every run that found the right one — the
ranking separates correct from incorrect, which is the property that actually
matters for prospective work.

Reproduce with:

```bash
python make_testcase.py     # builds test_data/ from a bundled complex
python validate.py 250000 8 # 8 runs x 250k evaluations, both functions
python rmsd_floor.py        # the 0.64 A floor
```

## Usage

### Command line

```bash
# site from a known ligand, which also reports an RMSD
pydock -r protein.pdb -l ligand.sdf --ref crystal.sdf

# site from coordinates
pydock -r protein.pdb -l library.sdf --centre 12.4 -3.1 28.7 --radius 11

# a whole library in one file: every record is docked and ranked
pydock -r protein.pdb -l library.sdf --ref crystal.sdf -n 10

# site from residues
pydock -r protein.pdb -l lig.sdf --site-residues A:129 A:201 A:330

# force field scoring, more effort, several ligand files
pydock -r protein.pdb -l 'actives/*.sdf' -s forcefield -n 20 --evals 400000
```

### Configuration file

GOLD-style `key = value`; `ligand` may be repeated.

```
receptor         = protein.pdb
ligand           = actives.sdf
ligand           = decoys.sdf
reference_ligand = crystal.sdf
scoring          = empirical
runs             = 10
poses            = 5
population       = 100
islands          = 5
evals            = 250000
output           = results
```

```bash
pydock -c dock.conf
```

### Python

```python
from pydock import DockingConfig, GAConfig, run_docking

cfg = DockingConfig(
    receptor="protein.pdb",
    ligands=["ligand.sdf"],
    reference_ligand="crystal.sdf",
    scoring="empirical",
    n_runs=10,
    ga=GAConfig(pop_size=100, n_islands=5, n_evals=250_000, autoscale=False),
)
result = run_docking(cfg)[0]
best = result.solutions[0]
print(best.score, best.terms)          # total and per-term breakdown
print(len(result.clusters))            # distinct binding modes
```

Lower-level pieces are usable on their own — score an existing pose without
searching, for instance:

```python
from pydock import prepare_receptor, load_ligand, build_scorer, site_from_ligand

rec = prepare_receptor("protein.pdb")
lig = load_ligand("pose.sdf")
centre, half = site_from_ligand(lig.coords0[lig.heavy_idx])
scorer = build_scorer("empirical", rec, lig, centre, half)
print(scorer.score(lig.coords0, breakdown=True), scorer.terms())
```

## Output

One SDF per ligand, poses ranked best first, each tagged with the total score,
the GA run that produced it, and every individual term:

```
pydock_rank = 1          pydock_hbond = 2.0
pydock_score = -15.1553  pydock_lipo = 183.3341
pydock_scoring = empirical   pydock_clash = 0.0267
pydock_run = 3           pydock_internal = 0.0
```

## Architecture

```
pydock/
  geometry.py   quaternions, axis rotations, torsion-tree transform, RMSD
  params.py     every parameter table and coefficient, in one place
  ligand.py     RDKit typing, rotatable bonds, torsion tree, symmetry maps
  receptor.py   PDB parsing, typing, polar-H placement, rotatable hydroxyls
  cavity.py     LIGSITE cavity detection
  grid.py       grid construction and trilinear interpolation
  scoring.py    the two scoring functions and the pairwise H-bond kernels
  ga.py         island GA, niching, Solis-Wets local search, clustering
  dock.py       orchestration, pose output
  cli.py        flags and config files
```

Chromosome layout, the same decomposition GOLD uses:

```
[ tx ty tz | qw qx qy qz | torsion_1..T | receptorOH_1..R ]
```

Torsions are emitted in breadth-first order from the largest rigid fragment, so
a child rotation's axis is already in its final position when it is applied.

### Design decisions worth knowing about

**Steric repulsion is split by the polarity of both partners.** Two polar atoms
at 2.6 Å are hydrogen bonded; two apolar atoms at 2.6 Å are clashing. A single
type-independent clash map would penalise every hydrogen bond it found, so
there are four clash maps instead.

**Intramolecular strain is repulsion-only.** An attractive intramolecular well
rewards the ligand for folding onto itself, which is not what a docking search
should be optimising.

**Cavity detection matters more than any other single feature.** On the
validation case it cut the search volume to 12% of the box and the receptor
hydroxyl genes from 35 to 20; RMSD went from 9.8 Å to 2.2 Å at identical
effort, and the run got 4× faster.

**No slicing in njit loops.** Taking `coords[i]` inside a Numba loop allocates
an array view per call. Rewriting the hydrogen-bond kernel with scalar indexing
made the whole score 2.7× faster with bit-identical output.

## Honest limitations

Things a production docking program has and this does not:

- **Parameters are published values, not fitted ones.** The coefficients come
  from the ChemScore and AutoDock4 papers. They were fitted to *those*
  functions as those authors implemented them, not to this code. Expect the
  absolute numbers to be unreliable as binding affinities; the useful signal is
  the ranking.
- **No CSD-derived statistics**, so no ASP-style potential and no torsion
  preference library. GOLD's torsional and hydrogen-bond geometry priors come
  from the Cambridge Structural Database, which is licensed data.
- **Receptor charges are approximate.** A compact Amber-flavoured table, not a
  real force field assignment. This is the likeliest reason the `forcefield`
  function underperforms the empirical one here.
- **Protonation is assumed, not predicted.** Histidine defaults to the neutral
  ND1-H tautomer; aspartate, glutamate, lysine and arginine are taken as
  charged. For anything pH-sensitive, protonate the receptor yourself first
  (PDB2PQR, `reduce`) — explicit hydrogens in the input file are used as given.
- **Protein backbone and side-chain heavy atoms are rigid.** Only hydroxyl
  rotation is sampled. No induced fit.
- **Waters are discarded by default** (`--keep-waters` keeps them, but they are
  then fixed, not toggled on and off the way GOLD can).
- **No covalent docking, no constraints, no pharmacophore restraints,
  no ensemble docking.**
- **Speed.** Roughly 25–35 s per GA run for a 7-rotor ligand on one core. About
  80 µs per score evaluation; the remaining overhead is Python in the GA step
  itself. Fusing a GA generation into a single Numba kernel is the obvious next
  optimisation and should be worth several times.
- **Validated on one complex.** One redocking success is a working proof, not
  a performance claim. A real assessment needs the Astex Diverse Set or
  CASF-2016, scored as a distribution.

## Next steps, roughly in order of value

1. **Benchmark properly** — Astex Diverse Set (85 complexes), report the
   fraction under 2 Å. That number is what makes the engine citable.
2. **Hydrogen-bond fitting points.** GOLD's biggest search advantage: it
   precomputes where a donor or acceptor *could* sit and maps ligand polar
   atoms onto those points instead of searching translation blindly.
3. **Fit the coefficients.** With a training set and the terms already
   separated, a linear fit against measured affinities is straightforward and
   turns borrowed parameters into real ones.
4. **Fuse the GA step into Numba**, and parallelise runs across cores — the
   runs are independent, so this is close to free.
5. **Side-chain flexibility** for a chosen set of residues, as extra genes.
6. **ML rescoring** of the surviving poses, trained on PDBbind, which is where
   most of the remaining accuracy is.

## Requirements

Python ≥ 3.10, NumPy, Numba, RDKit. Everything runs on CPU.

SciPy is deliberately *not* a dependency: it was there for one derivative-free
minimiser, and `pydock.optimize` now provides the same adaptive Nelder-Mead
(Gao & Han 2012) in about seventy lines. It reproduces SciPy's results on the
validation case exactly — same poses, same scores, run for run — and takes
109 MB off the packaged application.

## Licence

The code here is original. The parameter tables are published values cited in
`params.py`. Nothing is derived from GOLD, the CCDC, or the CSD.

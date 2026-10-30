# PISTAK — π-stack analysis from CIF files

*[Italiano](README.md) · English*

Stand-alone Python application based on the algorithm of the "π-stack analyser from CIF files" tool on
moleculesfriends.eu. It runs the same analyses as the website (interplanar distance, longitudinal and lateral slip,
angle between planes, contacts, overlap) and adds: extra parameters, stack columns, a ring-pair table (Janiak style),
processing of whole folders, synchronisable 3D and projection views of the parameters, and a plugin system.

**User manual**: [pistak/manual/en.md](pistak/manual/en.md) (Italian: [pistak/manual/it.md](pistak/manual/it.md)),
with a guided example, an explanation of every view and every quantity, the thresholds, the command line and the
library. It is also inside the application: *Help ▸ Manual* (F1). This README is the guide for developers.

## Starting

```bash
./start.py                       # GUI (with the demo crystal)
./start.py examples/herringbone_anthracene.cif
./start.py analyze a.cif --csv out.csv --json out.json     # terminal only (--block NAME for files with several structures)
./start.py batch folder/ --csv summary.csv                  # many files in parallel
./start.py compare a.cif b.cif c.cif --ref 1 --csv comparison.csv   # compare the stacking motif
./start.py stats folder/ --kind stack --corr --csv pairs.csv --summary description.csv
./start.py qc a.cif --program psi4 --method sapt -o inputs/ # stacked dimers for ORCA / Psi4 / Gaussian
./start.py check a.cif                                       # quick checks (formula, density, H, missed symmetry)
./start.py geom a.cif --csv geom.csv                         # bonds, angles, torsions, vdW contacts with s.u.
./start.py periodic a.cif --program qe -o calc/              # the cell for VASP / Quantum ESPRESSO / CP2K / CRYSTAL
./start.py session.pistak                                    # reopen a saved session
./start.py params                # available parameters (including those from plugins)
./start.py config                # where personal settings, last session and plugins live
```

The first time, `start.py` creates `.venv/` with a Python 3.11–3.13 (if your default one is newer it looks for a
suitable one, or use `./start.py --python /path/to/python3.12`) and installs `requirements.txt` into it. Starting
again is immediate; if you change the requirements it reinstalls by itself. Your own extra packages: put them in
`requirements-user.txt`. `./start.py --dev` adds pytest, `./start.py --reinstall` redoes the installation.

## How it is built

```
pistak/
  core/      pure computation (numpy/scipy), no GUI: usable as a library
    cif.py          CIF reading (data_ blocks, cell, symmetry, atoms, disorder rule)
    spacegroups.py  symmetry operations from the Hall or H-M symbol (when the CIF does not list them)
    structure.py    cell → symmetrised atoms → bonds → whole molecules
    cores.py        rings and π fragments of each molecule (planar ring systems), long axis from the atoms
    pairs.py        pairs of neighbouring π fragments, symmetry equivalence, classification (stack / parallel / tilted)
    stacks.py       stack columns (rank 0–3) from the pair graph
    ringpairs.py    ring–ring parameters (Cg···Cg, α, β, γ)
    config.py       thresholds (Limits) and preferences: pistak/default.toml + the user's user.toml
    errors.py       errors with a stable code (CifError, StructureError, InputError)
    stages.py       analysis stages and the thresholds each one depends on (partial re-analysis with cache)
    analysis.py     the whole chain, `reanalyse`, the PairContext and the per-molecule summary (MoleculeContact)
    uncertainty.py  parameter uncertainties with a Monte Carlo over the s.u. of the CIF
    crystal.py      cell content (with all disordered components), density, Z′, refinement values
    validation.py   quick checks (pre-checkCIF): A/B/C/G alerts with stable codes; ADDSYM with optional spglib
    intramolecular.py  bonds, angles, torsions and van der Waals contacts of the independent molecules, symmetry
                       codes (C2#1), s.u. propagated from the coordinates
    diagnosis.py    why an analysis has no pairs (saturated rings, no rings, distant cores…), for the empty views
    elements.py     atomic numbers, masses, radii, colours
  params/    registry of the measurable parameters (Registry; builtin = website, extra = phase 2, plugins) and
             registry of whole-structure quantities (structure.py: StructureRegistry, `@structure_parameter`)
  io/        CSV/JSON/XYZ export with fixed keys; batch (parallel), stats (statistics), compare (comparison);
             qc = dimer input for ORCA, Psi4, Gaussian; periodic = cell for VASP, Quantum ESPRESSO, CP2K, CRYSTAL
  gui/       window (PySide6): session.py = shared state and synchronisation; views/ = the views;
             limits_panel.py = live threshold sliders;
             collection_window.py = comparison and statistics; session_file.py = .pistak format
  errors.py  error messages in the current language
  labels.py  user-facing texts built from the results (pair relation, column type)
  i18n.py    translations (i18n_en.py: the English dictionary)
  cli.py     command-line commands
  default.toml  default thresholds and preferences
  manual/    the manual (it.md, en.md) and its figures (img/<language>/), also shown in Help ▸ Manual
plugins/     your own parameters (loaded at start-up)
tools/       manual_figures.py: regenerates the manual figures from the real window (no screen needed)
examples/    synthetic test CIFs (not experimental)
tests/       pytest, including parity with the original algorithm of the website
```

Dependency rule (checked by `tests/test_structure.py`): `core` imports only itself, no translations; `params` uses
`core` (and `i18n` for the parameter labels); `io` uses `core` and `params` and writes files with fixed English
headers, identical in every language; `gui` and `cli` use everything and translate with `labels`, `errors`, `i18n`.

**Personal settings.** They live outside the project, in the user's folder (`./start.py config` shows it):
macOS `~/Library/Application Support/PISTAK/`, Windows `%APPDATA%\PISTAK\`, Linux `~/.config/pistak/`, or the folder
given by `PISTAK_CONFIG_DIR`. It holds `user.toml` (thresholds and preferences that override `pistak/default.toml`),
`last.pistak` (the last session) and `plugins/`. Earlier versions used `config/` inside the project: on first start
`user.toml` and `last.pistak` are copied from there to the new folder.

## Windows and synchronisation

Every view (table, two-parameter scatter, histogram, 3D, projection, detail) is a panel that you can drag, stack or
detach into a window of its own; *View ▸ New view* opens further copies. At the top left of every view there is the 🔗
selector:

* **Group A/B/C**: the views of the same group follow the same selected pair and the same filter. Click on a point of
  the scatter → the table highlights the row, the 3D and the projection show that pair.
* **Not synchronised**: the view keeps a selection of its own (useful to compare two pairs side by side).

In the scatter, "Box selection" filters the pairs inside the rectangle; in the histogram a click on a bar filters that
interval; the filter applies to the group (the table shows only those pairs).

## Comparison and statistics across several structures

*Analysis ▸ Comparison and statistics…* (Ctrl+Shift+C) opens a window where you drag CIF files and folders
(analysed in parallel). All the pairs of all the files end up in the same table:

* **Structures**: for each one the closest π-stack pair, the differences (Δ) from the reference file (right click ▸
  "Set as reference") and the *motif distance*: distance in Å, in the (interplanar, L slip, W slip) space, between the
  two most similar stacking motifs. Double click: opens that file in the main window.
* **Distributions**: histogram of any parameter, stacked by pair type or by file.
* **Scatter**: two parameters of your choice, colour by type, file or a third parameter; a click on a point and "Open in
  the main window" goes straight to that pair.
* **Correlations**: Pearson matrix between the ticked parameters (only pairs with finite values; NaN if a parameter is
  constant).
* **Summary**: mean, standard deviation, minimum, quartiles and maximum; plus one row per file.

The "Pairs" filter (π-stack / parallel / tilted / all) and "Weight by equivalent pairs" apply to all tabs.
From the command line: `compare` and `stats` (see above). Everything can be exported to CSV.

## Sessions

*File ▸ Save session…* writes a `.pistak` file (readable JSON) with: the CIF (embedded, so the session opens even if you
move the file), the thresholds, the selected pair and the filter of each group, the views with their settings (chosen
axes, table columns, 3D options and camera…), the layout and geometry of the panels, and the list of files of the
comparison window. It is reopened from *File ▸ Open a session…*, by dragging the file onto the window, or with
`./start.py session.pistak`. On closing, the last session is saved automatically as `last.pistak` in the personal
settings folder (*File ▸ Reopen last session*). A new view that you write saves its settings by implementing
`extra_state()` and `apply_state()`; unknown keys from a newer file are ignored.

## Language and window sizes

* **Italiano / English:** *Language* menu (the window rebuilds itself keeping session, selection and panels), or
  `./start.py --lang en`, `PISTAK_LANG=en`, or `language = "en"` in `user.toml` (`[gui]` section; `auto` = system
  language). Texts are written in Italian inside `tr("...")` and translated by `pistak/i18n_en.py`;
  `python -m pistak.i18n_check` lists those without a translation. Plugins translate their own labels with
  `pistak.i18n.add_translations("en", {...})`.
* **Light / dark theme:** *View ▸ Theme* (automatic like the system, light, dark), or `background` in `user.toml`
  (`[gui]` section).
* **Windows:** the main window never exceeds the screen (94% of the usable area), can be shrunk down to about 640×420
  and enlarged at will; the view toolbars wrap instead of widening the window; sessions saved on a larger screen are
  brought back inside the screen.

## Extending

* **A new parameter**: copy `plugins/example_param.py` (into `plugins/` of the project or of the personal settings
  folder), change the function, *Analysis ▸ Reload plugins*. It appears everywhere (table, axes, export). Available
  context: see the docstring of the file and `core/analysis.py::PairContext`. From code you can also compute a different
  set of parameters with a registry of your own: `reg = Registry()`, `@reg.parameter(...)`,
  `compute_all(analysis, reg)`; parameter errors end up in `analysis.param_errors`.
* **A whole-structure quantity**: `@structure_parameter("key", "Label", unit=…, group=…)` on a function that receives
  the `Analysis` (also in plugins). It appears in Detail ▸ Structure, in `analyze`, in the JSON (`summary.structure`)
  and, if you add it to `STRUCTURE_KEYS` in `io/batch.py`, in the `batch` summary.
* **A new quick check**: `@check("name")` in `core/validation.py` returns a list of `Alert(code, level, values)`; the
  user-facing text goes in `ALERT_TEXT` of `labels.py` (and the translation in `i18n_en.py`).
* **A new view**: inherit from `gui/views/base.py::PairView` (synchronisation is already done) and add it to
  `VIEW_TYPES` in `gui/views/__init__.py`.
* **Thresholds**: *Analysis ▸ Thresholds…* or `user.toml` in the personal settings folder (`[limits]`).
* **Another 3D renderer**: the scene is described by `gui/scene.py` (numpy only); `views/view3d.py` draws it with
  PyVista.

## Analysis conventions

All of them are chosen in *Analysis ▸ Thresholds…*, in `user.toml` (`[limits]`) or from the command line
(`--pair-mode`, `--disorder`, `--overlap-measure`, `--min-overlap`, `--max-slip`). *Website values* brings them back to
those of the website (`Limits.site()`).

* **Pairs** (`pair_mode`). `canonical` (default): each contact between two molecules appears once; two contacts are
  the same pair if an operation of the space group maps one onto the other (not because the rounded values
  coincide); the long axis L is the bisector of the long axes of the two molecules, so longitudinal and lateral slip
  do not depend on which one is called A; "Equivalent pairs" counts the contacts per cell. The neighbouring cells
  explored widen by themselves as much as needed (small or very oblique cells), with `neighbour_shells` as a
  minimum. `site`: like the website — each contact twice (a→b and b→a, with slips that may differ), equivalence by
  rounded values, always ±`neighbour_shells` cells.
* **π fragments** (`fragments`, `min_fragment_atoms`). `all` (default): every system of fused planar rings of a
  molecule is a fragment (π1 = the largest, the usual π core; the others with at least `min_fragment_atoms` atoms).
  Pairs are between fragments: for two molecules in contact, π1–π1 and the other fragment combinations closer than
  `contact_cutoff` + 1 Å are measured; with several fragments in the same molecule the intramolecular pairs too. The
  relation in the table says which ones ("molecule 2 [1 0 0] · π2–π1", "intramolecular · π1–π3"); the parameters
  `fragment_a`/`fragment_b` give them as numbers. Symmetry equivalence also between fragments. `largest`: only π1
  (the website). The **per-molecule summary** (*Detail ▸ Contacts*, `molecule_contacts` in the JSON) gathers, for each
  contact between molecules, its fragment pairs, how many stacks, the shortest contact and the best pair (double click
  selects it).
* **Long axis** (`long_axis_atoms`). Empty: the direction of maximum extension of the fragment; for nearly isotropic
  cores (parameter *Core isotropy* close to 1: coronene, porphyrins) it is ill-defined, and so is the split of the
  slip into L and W. With two CIF labels ("C1,C9") the axis goes from the first atom to the second, projected onto the
  plane, for the fragments that contain them (or for π1). If a label appears more than once in the molecule (molecule
  on a symmetry element) the farthest pair is taken: "C1,C1" is fine.
* **Signed slips**: `long_slip_signed`, `lat_slip_signed` (components of s from A towards B along L and W) and
  `slip_direction` (−180…180°). The sign has physical meaning if the axis is defined by atoms and the molecule does not
  sit on an inversion centre (which reverses the axis); swapping A and B changes the sign.
* **Uncertainties (σ)**: *Analysis ▸ Uncertainties (Monte Carlo)…*, `analyze --mc 200` (`--mc-seed`, `--mc-no-cell`) or
  `pistak.api.uncertainty(analysis, 200)`. From the s.u. of the CIF (the digits in parentheses of the coordinates and
  the cell) N perturbed geometries are generated: the random displacement of each atom of the asymmetric unit is passed
  to its symmetry copies through the operation that generates them, the planes and axes of the fragments are
  recomputed and all pairs are measured again. The σ of every parameter of the registry (plugins included) is the
  standard deviation over the samples; in addition, for every pair, the fraction of samples in which the type
  (stack/parallel/tilted) stays the same. The Detail shows "value ± σ", the CSV adds the `<key>_su` and
  `kind_agreement` columns, the JSON `su` per pair. Approximations: independent coordinates (no covariance), special
  positions not constrained, atoms without s.u. kept fixed. A large σ on the L/W slips of a nearly isotropic core
  (phenyl, coronene) says that the long axis is not defined.
* **Dimer for quantum-chemistry calculations**: *Export ▸ Selected pair for quantum chemistry…* or
  `./start.py qc file.cif` (`--pair 1,3-5` or `--kind stack|parallel|tilted|all`, `-o` file or folder; without `-o` and
  with a single pair the input goes to the terminal). The monomers are the **two whole molecules** of the pair in the
  crystal geometry (B already translated into the right cell), not just the π fragments: for a contact between
  fragments the energy is still that of the whole dimer. X-ray hydrogens are too short: with normalisation (default,
  `--raw-h` to remove it) every X–H is lengthened along the bond to the mean neutron distances (Allen 1987/2010: C–H
  1.083, N–H 1.009, O–H 0.983, B–H 1.19, S–H 1.338 Å); heavy atoms do not move.
  - **DFT-D** (`--method dft-d`, default B3LYP-D3(BJ)/def2-TZVP; `--functional`, `--dispersion d3bj|d3zero|d4|none`,
    `--basis`): interaction energy with counterpoise (Boys–Bernardi). ORCA: three jobs in one file (`$new_job`):
    dimer, A with ghost B (`El :`), B with ghost A; E_int = E(AB) − E(A^AB) − E(B^AB). Psi4:
    `energy('b3lyp-d3bj', bsse_type='cp')`, which does everything by itself. Gaussian: `Counterpoise=2` with
    `(Fragment=n)` on the atoms (no D4: `EmpiricalDispersion=GD3BJ|GD3`). `--no-cp` gives the dimer calculation only.
  - **SAPT** (`--method sapt`, Psi4 only; `--sapt-level sapt0|ssapt0|sapt2|sapt2+|sapt2+(3)|sapt2+3`, default basis
    jun-cc-pVDZ, the usual choice for SAPT0): energy split into electrostatics, exchange, induction and dispersion.
  - Charges and multiplicities of the monomers: `--charges 0,0`, `--mult 1,1` (the dimer has the sum of the charges
    and, for high-spin coupling, 2S = 2S_A + 2S_B). `--nprocs`, `--mem` (GB) go into the program headers.
  - The file starts with a comment saying where it comes from (CIF, block, pair, symmetry, distance, hydrogens) and is
    all ASCII. Intramolecular pairs (between fragments of the same molecule) are not exported: they are not a dimer.
    The keywords follow the current versions (ORCA 5–6, Psi4 1.7+, Gaussian 16): check them against yours.
* **Staged analysis and live thresholds**: the chain is split into stages (`core/stages.py`) — CIF reading
  (`disorder`), structure (bonds, cell), π fragments (rings, long axis), contacts (`contact_cutoff`,
  `neighbour_shells`, `pair_mode`: the measurement of all pairs, the expensive part) and pairs (the stack criterion:
  `max_dihedral`, `max_interplanar`, `min_overlap`, `max_slip`, `overlap_measure`). Changing the thresholds redoes
  only from the first stage touched (`pistak.api.reanalyse(analysis, limits)`; the window does it by itself): with
  those of the stack criterion the pairs already measured are reclassified without redoing the geometry (a few ms even
  with tens of pairs), and the registry parameters declared `geometric=True` (all those of `params/extra.py`) are
  reused. *Analysis ▸ Live thresholds* (Ctrl+Shift+L) opens the sliders: table, scatter (with the thresholds as dashed
  lines), histogram and columns follow the slider, selection and filters stay on the same pairs; the maximum contact
  measures again and applies on release, in a separate thread. *Save as default* writes `user.toml`. The Monte Carlo
  uncertainties have to be recomputed after a change of thresholds. A plugin that depends only on the geometry of the
  pair can declare it with `@parameter(..., geometric=True)`.
* **Disordered atoms** (`disorder`): `major` (default) keeps, for each assembly (`_atom_site_disorder_assembly`), the
  group with the highest mean occupancy; `first` group 1 (the website); `all` all of them; a number, that group.
* **Stack criterion**: besides the angle between planes and the interplanar distance, the overlap of the projected
  cores must exceed `min_overlap`, measured as the fraction of atoms inside the outline of the other core (`atoms`,
  the website) or as common area / area of the smaller core (`area`); with `max_slip` > 0 pairs that are too slipped
  are also discarded.
* **Files with several structures**: every `data_` block with atoms is a structure. The window asks which one to open,
  the comparison and `batch`/`stats`/`compare` analyse them all (named "file.cif [block]"), `analyze` takes the first
  or the one given by `--block`.
* **Missing symmetry operations**: if the CIF does not list them they are derived from the Hall symbol
  (`_space_group_name_Hall`) or the H-M symbol (table of the most common groups in `core/spacegroups.py`; with the
  optional `gemmi` package installed, all of them). The JSON summary says where they come from (`symops_source`).

## Parity with the website

`tests/test_parity.py` compares the port, with `Limits.site()`, with the numbers produced by the original TypeScript
sources of the website on 6 synthetic crystals (all three pair types). To regenerate the reference values:
`SITE_DIR=/path/to/MYSITE node tests/parity/dump_site.js`.

Intentional differences even in `site` mode: the limit of atoms in the cell is 6000 (2500 in the browser); the
directions of the axes (L, W) are fixed in a reproducible way; in addition there are the extra parameters and columns.

## Development

`./start.py --dev --setup-only` installs pytest and ruff into `.venv/`. The project is also an installable package
(`pyproject.toml`, `pistak` command): `pip install -e ".[gui,dev]"`.

```bash
.venv/bin/python -m pytest -m "not gui"     # computation only, a few seconds
.venv/bin/python -m pytest                  # everything, including the windows (offscreen)
.venv/bin/ruff check .                      # static check (rules in pyproject.toml)
.venv/bin/python tools/manual_figures.py    # regenerate the manual figures (after interface changes)
```

**Manual.** `pistak/manual/it.md` and `en.md` must be kept aligned: `tests/test_manual.py` checks that they have the
same chapters and the same figures, that internal links lead to a heading, and that the library example really runs.
The figures are regenerated with the script above (only the CIFs in `examples/`, 256-colour PNG).

Markers: `gui` (windows), `parity` (comparison with the website), `snapshot` (snapshots). **Snapshots**:
`tests/snapshots/site/*.json` and `tests/snapshots/canonical/*.json` contain all the results on the examples and on the
real COD structures in `tests/data/cod/` (public domain) in the two modes (molecules, cores, pairs with all the
parameters, ring pairs, columns, structure quantities, alerts, van der Waals contacts); any change of a number makes
`pytest -m snapshot` fail. After an intended change: `pytest -m snapshot --update-snapshots` and check the diff. For
real structures that cannot be distributed (CSD): `PISTAK_SNAPSHOT_DIR=/folder/cif pytest -m snapshot --update-snapshots`
writes to `tests/snapshots/local/<mode>/`, which stays out of the repository.

## Notes for exFAT/FAT volumes (macOS)

macOS creates `._name` files next to copied files. `start.py` removes them from `.venv` after every installation (they
break matplotlib) and the plugin/CIF loaders ignore them. If you prefer, move the project to an APFS disk.

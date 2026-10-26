# PISTAK — π-stack analysis from CIF files

Stand-alone Python application based on the algorithm of the "π-stack analyser" on moleculesfriends.eu.
Same analyses as the website plus extra parameters, stack columns, ring-pair tables, folder processing,
comparison and statistics across structures, synchronised 3D / projection / parameter views, sessions and plugins.

```bash
./start.py                       # GUI (creates .venv and installs packages the first time)
./start.py --lang en             # English interface (also PISTAK_LANG=en, or the Language menu)
./start.py analyze a.cif --csv out.csv
./start.py batch folder/ --csv summary.csv
./start.py compare a.cif b.cif --csv comparison.csv
./start.py stats folder/ --kind stack --corr
./start.py session.pistak        # reopen a saved session
./start.py config                # where personal settings, last session and plugins live
```

Exported files (CSV/JSON/XYZ) always use fixed English column names, whatever the interface language.
By default each intermolecular contact is listed once and contacts related by a space-group operation are grouped
(`pair_mode = "canonical"`); `--pair-mode site` (or *Website values* in the thresholds dialog) reproduces the website.
Molecules with several planar ring systems get one π fragment each (`fragments = "all"`): pairs are fragment–fragment
(π1–π2, also intramolecular), grouped per molecule contact in *Detail ▸ Contacts* and in the JSON. The long axis can be
set from two atom labels (`long_axis_atoms = "C1,C9"`), and signed slips are available as parameters.
Parameter uncertainties (σ) come from a Monte Carlo over the standard uncertainties in the CIF (*Analysis ▸
Uncertainties*, `analyze --mc 200`, `pistak.api.uncertainty`), symmetry-consistent, for every registry parameter.
Disorder policy (`--disorder major|first|all|N`), stack criterion (`--overlap-measure`, `--min-overlap`, `--max-slip`),
multi-block CIFs (`--block`) and missing symmetry operations (derived from the Hall or H-M symbol) are described in `README.md`.
Personal settings (`user.toml`), the last session and a `plugins/` folder live in the user's config folder
(macOS `~/Library/Application Support/PISTAK/`, Windows `%APPDATA%\\PISTAK\\`, Linux `~/.config/pistak/`, or `PISTAK_CONFIG_DIR`).

Windows fit the screen (never larger than 94% of the usable area), can be resized freely, and toolbars wrap.
The full documentation is in `README.md` (Italian); code comments are in Italian, UI texts are translated via
`pistak/i18n_en.py` (`python -m pistak.i18n_check` lists missing translations).

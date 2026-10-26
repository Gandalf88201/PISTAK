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
Personal settings (`user.toml`), the last session and a `plugins/` folder live in the user's config folder
(macOS `~/Library/Application Support/PISTAK/`, Windows `%APPDATA%\\PISTAK\\`, Linux `~/.config/pistak/`, or `PISTAK_CONFIG_DIR`).

Windows fit the screen (never larger than 94% of the usable area), can be resized freely, and toolbars wrap.
The full documentation is in `README.md` (Italian); code comments are in Italian, UI texts are translated via
`pistak/i18n_en.py` (`python -m pistak.i18n_check` lists missing translations).

# Review log — English summary

`validation/REGISTRO.md` is the full review log of PISTAK, in Italian: one row per problem found or check closed,
with the reference used and the effect on the results. This file summarizes every entry in English and gives a
glossary of the Italian file names, so that the entries cited in the paper can be checked without reading Italian.
Where the two differ, `REGISTRO.md` is authoritative.

Status words: *corrected* (code changed), *verified* (checked, no change needed), *option* (new setting, default
unchanged), *deferred* (left to the author), *open*.

## Glossary of Italian file and folder names

| Name | Meaning |
|---|---|
| `REGISTRO.md` | review log (register of checks) |
| `check_struttura.py` | structure checks: cell, formulae, oblique cell and origin invariance, rings, RDKit, radii |
| `check_incertezze.py` | uncertainty checks: Monte Carlo against linear propagation, convergence, `_geom_` tables |
| `check_sensibilita.py` | threshold sensitivity |
| `check_riproducibilita.py` | reproducibility (1 *vs* 4 processes, Monte Carlo seeds) |
| `check_simmetria.py` | symmetry (Hall and Hermann–Mauguin symbols against gemmi) |
| `check_cif.py` | CIF reading (sites, s.u., elements, disorder) against gemmi |
| `check_registro.py` | small checks of individual register entries |
| `check_export.py` | QC and periodic export files |
| `check_snapshot.py` | stored reference results recomputed with NumPy |
| `check_platon.py` | ring–ring parameters against PLATON |
| `rapporto_platon.py` | report of the PLATON comparison |
| `scarica_cod.py` | downloads the 200 COD CIFs of `cod_sample.txt` |
| `cod_cache/`, `benchmark_cache/` | downloaded CIFs |
| `platon/`, `platon_polimeri/` | PLATON listings for molecular crystals and for coordination polymers |
| `ambiente.txt` | exact package versions used for the paper |
| `dichiarazione_ia.md` | draft statement on the use of AI |
| `LEGGIMI.md` | read-me |
| `confronto.txt` | comparison output |
| Italian words in outputs | *coppie* pairs; *anelli* rings; *abbinate* matched; *scarto* difference; *mediana* median; *uguale* equal; *diversi* different; *stack* π-stacked; *predefinito* default; *strutture* structures; *nuove/perse* new/lost; *soglia* threshold; *valore* value; *cambiano tipo* change type |

## Entries

| # | Area | Summary | Status |
|---|---|---|---|
| 1 | parameters | `beta_b` was labelled γ (numbers right, key wrong) | corrected |
| 2 | parameters | `mean_atom_height` equalled the interplanar distance (order of mean and absolute value) | corrected |
| 3 | rings | the one-cycle-per-bond ring search lost the central ring of coronene in some atom orderings; now all shortest cycles per bond | corrected |
| 4 | geometry | `unit()` of a null vector: shown to be unreachable | verified |
| 5 | pairs | the cap of ±8 cells in the neighbour search was silent; now alert `shellsTruncated` | corrected |
| 6 | pairs | contact cut-off 4.0 Å made the 4.0–4.2 Å range of the stack criterion unreachable; default now 4.5 Å | corrected |
| 7 | QC export | X–H distances for P and Si without a source; C–H by hybridization from Allen & Bruno (2010) | corrected |
| 8 | repository | licence file missing; MIT added | corrected |
| 9 | tests | the synthetic-crystal generator uses PISTAK's cell vectors: shown not to be circular | verified |
| 10 | CIF parser | infinite loop on non-ASCII white space | corrected |
| 11 | space groups | setting suffixes (`:1`, `:2`, `:R`) were dropped | corrected |
| 12 | CIF parser | element from labels (`BR2` read as boron, centroids and Q peaks kept as atoms, `dum` flags ignored) | corrected |
| 13 | quick checks | listed operations inconsistent with the space-group symbol: now alert `symopsMismatch` | corrected |
| 14 | disorder | negative disorder groups (1/−1 used as alternatives; groups on special positions) | corrected |
| 15 | cell contents | positions occupied by two different CIF atoms counted once | corrected |
| 16 | radii | covalent and van der Waals radii checked against their sources (C, Lu/Hf, H) | corrected |
| 17 | structure | reconstructed formulae against `_chemical_formula_moiety` | corrected |
| 18 | geometry tables | torsions through three collinear atoms (undefined) were reported | corrected |
| 19 | geometry tables | measurements with equal labels were taken as symmetry-equivalent | corrected |
| 20 | quick checks | false missing-H alerts on isonitriles, cyanides, metal carbonyls | corrected |
| 21 | quick checks | false ADDSYM alerts on disordered structures | corrected |
| 22 | documentation | A/B/C/G alert levels are PISTAK's own, not checkCIF alerts | corrected |
| 23 | quick checks | C≡X and sp-angle thresholds now cited | verified |
| 24 | QC export | dimer spin multiplicity is now an option | corrected |
| 25 | QC export | ORCA `%maxcore` margin | corrected |
| 26 | QC/periodic export | provenance and assumptions written in the input files | corrected |
| 27 | QC/periodic export | running the inputs with the target programs | deferred |
| 28 | tests | tolerances tightened to 10⁻⁴ Å and 10⁻³° | corrected |
| 29 | tests | stored snapshots recomputed with NumPy | verified |
| 30 | validation | invariance to origin shift and oblique cell | verified |
| 31 | configuration | user settings file silently changed analyses; tests and paper scripts isolated | corrected |
| 32 | tests | coverage measured (90 % overall, 97 % without the GUI); property-based tests not written | partly open |
| 33 | tests | test names translated to English | corrected |
| 34 | quick checks | spglib failures were reported as "spglib missing" | corrected |
| 35 | Monte Carlo | NaN in some samples gave a σ over the remaining samples silently | corrected |
| 36 | radii | 48 elements had no van der Waals radius; Alvarez (2013) for d and f blocks | corrected |
| 37 | GUI | licence notice of code adapted from a Qt example | corrected |
| 38 | documentation | dependency licences completed | corrected |
| 39 | environment | exact versions recorded in `ambiente.txt` | corrected |
| 40 | code quality | ruff clean | verified |
| 41 | paper | AI statement drafted | approved as draft |
| 42 | repository | continuous integration added (Linux, Python 3.11–3.13) | corrected |
| 43 | rings | rings compared with RDKit | verified |
| 44 | validation | threshold sensitivity | verified |
| 45 | repository | version for the paper | closed |
| 46 | ordering | ordering of symmetry-equivalent results differed between macOS and Linux | corrected |
| 47 | extended networks | repeat unit of coordination polymers split a ring between cells; ring pairs of polymers checked against PLATON | corrected |
| 48 | parameters | contact parameters of intramolecular pairs measured the molecule against itself | corrected |
| 49 | classification | two fragments bonded to each other were accepted as π-stacked (chelate rings sharing a metal) | corrected |
| 50 | classification | default overlap criterion (one atom inside the hull) accepted minimal overlaps; now area of the union of ring polygons > 0.10 | corrected |
| 51 | rings | η⁵/η⁶ rings were hidden by M–C–C triangles; two-pass ring search | corrected |
| 52 | fragments | `fragment_merge = "conjugated"` for coplanar non-fused systems | option |
| 53 | QC export | warnings for disorder, solvent masks; r²SCAN-3c added | partly corrected |
| 54 | Monte Carlo | percentile intervals; alert `noSU` | corrected |
| 55 | pair frame | undefined long axis (isotropy (λ₂/λ₃)^½ > 0.9) and twist limit (75°); new `twist`, `planes_cross`, `local_interplanar` | corrected |
| 56 | fragments | linear skeletons rejected as fallback fragments | corrected |
| 57 | radii | Br van der Waals radius from Mantina et al. (2009) | corrected |
| 58 | classification | ferrocene Cp rings were an intramolecular π-stack; fragments bonded to a common atom are now *bridged* | corrected |
| 59 | fragments | no π fragments in structures without C–C or C–H bonds (`pi_scope`) | corrected |
| 60 | structure | positions within 0.3 Å merged only for the same site or repeated ordered atoms | corrected |
| 61 | Monte Carlo | displacements of atoms on special positions projected onto the site symmetry | corrected |
| 62 | Monte Carlo | per-fragment alert `fragmentNoSU` | corrected |
| 63 | overlap | self-overlapping ring sets (convex hull; found on COD 2100449, whose cause was entry 65, kept as a safeguard), macrocycle triangles over rings, `planes_cross` for coplanar fragments | corrected |
| 64 | fragments | duplicated rings in the conjugated merge (H₂Pc) | corrected |
| 65 | structure | after entry 60, the same atom written twice at half occupancy (COD 2100449) was no longer merged; positions coinciding within 0.05 Å are merged whatever their occupancy | corrected |

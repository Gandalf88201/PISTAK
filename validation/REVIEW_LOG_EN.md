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
| 66 | ordering | pairs with the same interplanar distance up to the last binary digit swapped places between machines (COD 2008099 on Linux, Python 3.12); distances rounded to 10⁻⁶ Å in the sort key, deterministic tie-break | corrected |
| 67 | fragments | `pi_scope = "organic"` also searches structures with no C–C or C–H bond that have a planar ring of non-metal atoms with B, C or N (melamine, cyanuric acid, borazines); carbides, carbonates, cyanides, silicate-type rings still skipped. COD sample: only 7206203 (melaminium salt) enters, +6 pairs, no new π-stack | corrected |
| 68 | Monte Carlo | special positions: one deviate per free parameter with its own s.u., constraints give the dependent coordinates (the averaging projector gave σ/√2 along x = y); no σ of the 102 π-stacked pairs of the COD sample changes | corrected |
| 69 | definitions | isotropy and linearity test: the code has always used (λ₂/λ₃)^½ (since 26 October 2026) with thresholds 0.9 and 0.05; only the round-2 text said λ₂/λ₃. Applied to the ratio, 0.9 would give 47 instead of 52 of 72 undefined axes, 0.05 would drop 2 fallback fragments (1 π-stack) | verified |
| 70 | fragments | the planar-skeleton fallback accepts only sp²/sp skeletons (no atom with four or more neighbours, H included); fallback fragments 55 → 10 in the COD sample, pairs 4450 → 4284, π-stacked pairs 102 → 99 (all three removed involved sp³ solvents); `check_incertezze.py converge` selects its pairs explicitly | corrected |
| 71 | Monte Carlo | cell parameters: only the free ones are drawn, constraints derived from the symmetry operations (G = RᵀGR); alert `cellConstraints` when the CIF cell is not invariant. COD sample: all 200 cells consistent (no alert); cell deviates change in 8 structures, only 7102228 has π-stacked pairs (σ changes ≤ 0.8 %); medians over the 99 π-stacked pairs unchanged | corrected |
| 72 | reproducibility | compatibility-mode snapshot of the carbide COD 1538535 (R-3m): 6 pairs on macOS, 7 on Linux (contacts exactly on lattice translations, a tie decided by the last binary digit); that snapshot is excluded, the default-mode snapshot is kept and identical on both platforms; root cause still open | open (worked around) |

## Round 4 of the paper review (10–11 November 2026): items 73–87 (v1.3.0, not yet released)

- **73** Fallback long axis for nearly isotropic fragments (single rings): long axis of the whole molecule projected
  onto the fragment plane (`axis_fallback = "molecule"`, default), of the fragment with its substituents, or none.
  The origin of the axis is written with every pair (`long_axis_a`, `long_axis_b`).
- **74** Slips in the frame of each fragment (`long_slip_a`, `lat_slip_a`, `long_slip_b`, `lat_slip_b`): the pitch
  and roll distances of Curtis et al. (2004); defined whatever the twist between the two long axes.
- **75** `max_axis_twist` lowered from 75° to 30°: beyond it the bisector frame describes neither molecule.
- **76** *Clamped* pairs: two fragments of one molecule closer than `min_interplanar` = 3.2 Å are flagged and are not
  π-stacked (COD 7239588). Effect on the COD sample: 99 → 98 π-stacked pairs.
- **77** Complete symmetry classes: a class missed members when symmetry exchanged fragments of equal size, so its
  multiplicity was too low and its representative changed with `contact_cutoff`. Every class now holds its whole
  orbit and has a canonical representative. The number of pairs does not change.
- **78** Contacts per molecule (`per_molecule`), independent of Z.
- **79** Monte Carlo: optional covariance matrix of the coordinates (`--mc-covariance`); the reported value is the
  nominal one; results are reproducible from the seed.
- **80** Dimer export: keywords stated for ORCA 6; warnings for impossible charge/multiplicity, metals and salts.
  The exported inputs have not been run.
- **81** Output terms aligned with the paper (π-stacked, angle); no "interplanar distance" printed for pairs without
  projected overlap; JSON format version 5.
- **82** Windows and macOS added to continuous integration (to be confirmed after the first push).

- **83** Planar molecules without rings (nitrate, guanidinium) are π fragments only on request
  (`ringless_fragments = "planar"`; default `"none"`): with 500 structures the rule classed two nitrate–arene
  contacts (COD 7032025, 7104713) as π-stacked. Author's decision.
- **84** Two elements on one site (mixed occupancy, COD 2108505): angles and torsions through both atoms of the
  site are degenerate and no longer listed.
- **85** COD sample enlarged from 200 to 500 structures (the first 500 of the same seeded draw; the previous 200 are
  its first 200). The Monte Carlo over the sample runs in blocks of 100 structures.
- **86** (11 November 2026) Symmetry classes: when no molecule had more than one fragment, an operation mapping a
  molecule onto a copy of itself was assumed to map its fragment onto the fragment. This fails for a molecule on a
  symmetry element whose fragment is one of two equivalent rings of which only one is a fragment (COD 7044129), and
  for disorder removed asymmetrically (COD 4326386): non-equivalent pairs were merged into one class. The fragment
  maps are now always used. Found by the independent multiplicity check (2 of 364 structures of discrete
  molecules). COD sample: 10 875 → 10 883 pairs, no π-stacked pair changes.
- **87** (12 November 2026) The table of the graphical interface prints a dash, like the terminal, in place of the
  "interplanar distance" of pairs that are not π-stacked and have no projected overlap; one rule for both. The value
  stays in the scatter plot, CSV and JSON.
- **88** (12 November 2026, round 5, N1) Extended networks: no result depends any more on how the repeat unit is
  cut. The *clamped* flag and the pair categories now test membership of the same bonded network (a translation in
  the lattice of periodicities of the network), not of the same repeat unit; a bridging atom is looked up in the bond
  graph of the infinite network, not only in the unit; the nodes of the stacking graph are the discrete molecules
  and the π fragments of networks; and the symmetry images of fragments that share an atom are found whichever
  fragment the atom is listed in. A random re-cut of the repeat units of the 98 network structures of the sample
  (three cuts each; `paper/scripts/recut_invariance.py`) changed 22 of them with the previous rules (12 cuts with a
  different motif rank) and none now. Test: `tests/test_round5.py`.

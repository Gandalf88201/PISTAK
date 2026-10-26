"""Comandi da terminale. Senza argomenti (o con un file .cif) apre la finestra grafica.

    python -m pistak                       finestra grafica
    python -m pistak struttura.cif         finestra grafica con quel file
    python -m pistak analyze a.cif [--block NOME] [--csv out.csv] [--json out.json] [--xyz 1 out.xyz]
    python -m pistak batch cartella/ --csv riepilogo.csv
    python -m pistak compare a.cif b.cif c.cif [--ref 1] [--csv confronto.csv]
    python -m pistak stats cartella/ [--kind stack] [--csv coppie.csv] [--summary descrizione.csv] [--corr]
    python -m pistak gui sessione.pistak   riapre una sessione salvata dalla finestra
    python -m pistak params                elenco dei parametri disponibili (con quelli dei plugin)
    python -m pistak demo > demo.cif       scrive il cristallo dimostrativo
    python -m pistak config                dove stanno impostazioni personali, ultima sessione e plugin
    python -m pistak qc a.cif --pair 1 --program psi4 --method sapt -o dimero.in
                                           il dimero di una coppia come input per ORCA, Psi4 o Gaussian
    python -m pistak --lang en ...         interfaccia in inglese (anche PISTAK_LANG=en)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pistak import i18n
from pistak.core.config import Limits, load_limits
from pistak.errors import ERROR_TEXT, describe_error  # noqa: F401  (ERROR_TEXT: per chi lo importava da qui)
from pistak.i18n import tf, tr
from pistak.labels import relation


def _limits(args) -> Limits:
    lim = load_limits()
    over = {k: getattr(args, k) for k in ("max_dihedral", "max_interplanar", "contact_cutoff", "pair_mode", "disorder",
                                          "overlap_measure", "min_overlap", "max_slip") if getattr(args, k, None) is not None}
    return Limits.from_dict(over, lim) if over else lim


def _add_limit_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--max-dihedral", type=float, help=tr("angolo massimo tra i piani (°)"))
    p.add_argument("--max-interplanar", type=float, help=tr("distanza interplanare massima (Å)"))
    p.add_argument("--contact-cutoff", type=float, help=tr("contatto massimo tra atomi pesanti (Å)"))
    p.add_argument("--pair-mode", choices=["canonical", "site"], help=tr("canonical: una coppia per contatto (predefinito); site: come il sito"))
    p.add_argument("--disorder", help=tr("atomi disordinati: major | first | all | numero del gruppo"))
    p.add_argument("--overlap-measure", choices=["atoms", "area"], help=tr("come si misura la sovrapposizione per il criterio di stack"))
    p.add_argument("--min-overlap", type=float, help=tr("sovrapposizione minima per uno stack (0–1)"))
    p.add_argument("--max-slip", type=float, help=tr("slittamento massimo per uno stack (Å, 0 = nessun limite)"))


def cmd_analyze(args) -> int:
    from pistak import analyse_file
    from pistak.core.analysis import read_cif_file
    from pistak.core.cif import block_names
    from pistak.io import export

    try:
        names = block_names(read_cif_file(args.file))
        a = analyse_file(args.file, _limits(args), args.block)
    except Exception as e:  # noqa: BLE001
        print(tf("errore: {msg}", msg=describe_error(e)), file=sys.stderr)
        return 1
    if len(names) > 1:
        print(tf("il file contiene {n} strutture ({names}); analizzo {block} (scegli con --block)",
                 n=len(names), names=", ".join(names), block=a.cif.name), file=sys.stderr)
    s = export.summary_dict(a)
    c = s["cell"]
    print(tf("{name}  gruppo {sg}  a={a:.3f} b={b:.3f} c={c:.3f}  α={al:.2f} β={be:.2f} γ={ga:.2f}  atomi in cella: {n}",
             name=a.name, sg=s["space_group"] or "?", a=c["a"], b=c["b"], c=c["c"], al=c["alpha"], be=c["beta"], ga=c["gamma"],
             n=s["cell_atoms"]))
    for m in s["molecules"]:
        frags = "" if m["pi_fragments"] <= 1 else ", " + tf("{n} frammenti π per molecola", n=m["pi_fragments"])
        print(tf("  {count} × {formula}  ({core} con nucleo π", count=m["count"], formula=m["formula"], core=m["with_pi_core"])
              + frags + ")")
    if not a.pairs:
        print(tr("nessuna coppia di molecole vicine con nuclei π"))
    else:
        print(f"{'#':>3} {tr('tipo'):<8} {tr('interpl.'):>8} {'L-slip':>7} {'W-slip':>7} {tr('diedro'):>7} {'Cg-Cg':>7} "
              f"{tr('contatto'):>8}  {tr('relazione')}")
        for p in a.pairs[: args.rows]:
            print(f"{p.index + 1:>3} {p.kind:<8} {p.interplanar:8.2f} {p.long_slip:7.2f} {p.lat_slip:7.2f} {p.dihedral:7.1f} "
                  f"{p.centroid_distance:7.2f} {p.contact:8.2f}  {relation(p)}")
        if len(a.pairs) > args.rows:
            print(tf("… altre {n} coppie (usa --rows o --csv)", n=len(a.pairs) - args.rows))
    if args.mc:
        _print_uncertainty(a, args)
    if args.csv:
        print(tr("scritto"), export.write_csv(a, args.csv))
    if args.json:
        print(tr("scritto"), export.write_json(a, args.json))
    if args.xyz:
        n, out = int(args.xyz[0]), args.xyz[1]
        if not 1 <= n <= len(a.pairs):
            print(tf("errore: la coppia {n} non esiste", n=n), file=sys.stderr)
            return 1
        print(tr("scritto"), export.write_pair_xyz(a, a.pairs[n - 1], out))
    return 0


def _print_uncertainty(a, args) -> None:
    from pistak.api import uncertainty

    prog = lambda i, n: print(f"\r[MC {i}/{n}]", end="", flush=True, file=sys.stderr)  # noqa: E731
    u = uncertainty(a, args.mc, args.mc_seed, include_cell=not args.mc_no_cell, progress=prog)
    print(file=sys.stderr)
    if not u.has_su:
        print(tr("attenzione: il CIF non ha incertezze standard (s.u.): tutte le σ sono 0"), file=sys.stderr)
    print("\n" + tf("incertezze (σ, Monte Carlo su {n} campioni dalle s.u. del CIF):", n=u.samples))
    print(f"{'#':>3} {tr('interpl.'):>15} {'L-slip':>15} {'W-slip':>15} {tr('diedro'):>13} {tr('stabilità del tipo'):>20}")
    for p in a.pairs[: args.rows]:
        s = u.sigma[p.index]
        print(f"{p.index + 1:>3} {p.interplanar:7.3f}±{s['interplanar']:<7.3f} {p.long_slip:7.3f}±{s['long_slip']:<7.3f} "
              f"{p.lat_slip:7.3f}±{s['lat_slip']:<7.3f} {p.dihedral:6.2f}±{s['dihedral']:<6.2f} "
              f"{100 * u.kind_agreement[p.index]:>19.0f}%")


def _entries(args):
    from pistak.io.batch import find_cifs, run_batch

    files: list[Path] = []
    for item in args.inputs:
        files += find_cifs(item)
    if not files:
        print(tr("nessun file .cif trovato"), file=sys.stderr)
        return None
    prog = lambda i, n, name: print(f"[{i}/{n}] {name}", flush=True, file=sys.stderr)  # noqa: E731
    return run_batch(files, _limits(args), workers=getattr(args, "workers", 0), progress=prog)


def cmd_batch(args) -> int:
    from pistak.io.batch import write_summary_csv

    args.inputs = [args.folder]
    results = _entries(args)
    if results is None:
        return 1
    out = write_summary_csv(results, args.csv)
    ok = sum(r.error is None for r in results)
    print(tf("{ok}/{n} file analizzati; riepilogo in {out}", ok=ok, n=len(results), out=out))
    return 0 if ok else 1


def cmd_compare(args) -> int:
    from pistak.io.compare import COMPARE_KEYS, compare, write_csv

    results = _entries(args)
    if results is None:
        return 1
    ref = min(max(args.ref, 1), len(results)) - 1
    rows = compare(results, ref)
    keys = ["interplanar", "long_slip", "lat_slip", "dihedral"]
    print(tf("riferimento: {name}", name=results[ref].name))
    print(f"{'file':<34} {'stack':>5} {'col':>3} " + " ".join(f"{k[:9]:>9}" for k in keys) + "   " + tr("Δ rispetto al rif.") + "            "
          + tr("dist.motivo"))
    for r in rows:
        if r.error:
            print(f"{r.name:<34} {tr('ERRORE')}: {r.error}")
            continue
        b = " ".join(f"{r.best[k]:9.2f}" if k in r.best else f"{'–':>9}" for k in keys)
        d = " ".join(f"{r.delta[k]:+6.2f}" if k in r.delta else f"{'–':>6}" for k in keys[:3])
        mm = "–" if r.motif_distance != r.motif_distance else f"{r.motif_distance:.2f} Å"
        print(f"{r.name:<34} {r.n_stack:>5} {r.n_columns:>3} {b}   {d}   {mm:>10}")
    if args.csv:
        print(tr("scritto"), write_csv(rows, args.csv, COMPARE_KEYS))
    return 0


def cmd_stats(args) -> int:
    from pistak.io.stats import PairTable, correlation, describe, write_describe_csv

    results = _entries(args)
    if results is None:
        return 1
    table = PairTable.from_entries([r for r in results if not r.error])
    if not len(table):
        print(tr("nessuna coppia trovata"), file=sys.stderr)
        return 1
    kind = None if args.kind == "all" else args.kind
    keys = args.keys.split(",") if args.keys else ["interplanar", "long_slip", "lat_slip", "dihedral", "contact", "overlap_area"]
    rows = describe(table, keys, kind, weighted=not args.unweighted)
    print(tf("{n} coppie in {f} file; tipo: {kind}", n=len(table), f=len({*table.files}), kind=args.kind))
    heads = [tr("media"), tr("dev.st."), "min", "q25", tr("mediana"), "q75", "max"]
    print(f"{tr('parametro'):<20}{'n':>6}" + "".join(f"{h:>9}" for h in heads))
    for r in rows:
        if not r["n"]:
            print(f"{r['key']:<20}{0:>6}")
            continue
        print(f"{r['key']:<20}{r['n']:>6}{r['mean']:9.3f}{r['std']:9.3f}{r['min']:9.3f}{r['q25']:9.3f}{r['median']:9.3f}{r['q75']:9.3f}{r['max']:9.3f}")
    if args.corr:
        c = correlation(table, keys, kind)
        print("\n" + tr("correlazioni (Pearson):"))
        print(" " * 20 + "".join(f"{k[:8]:>9}" for k in keys))
        for k, row in zip(keys, c):
            print(f"{k:<20}" + "".join(f"{'–':>9}" if v != v else f"{v:9.2f}" for v in row))
    if args.csv:
        print(tr("scritto"), table.write_csv(args.csv))
    if args.summary:
        print(tr("scritto"), write_describe_csv(rows, args.summary))
    return 0


def cmd_params(_args) -> int:
    from pistak.api import ensure_plugins
    from pistak.params import load_errors, specs

    ensure_plugins()
    group = None
    for s in sorted(specs(), key=lambda s: (s.group, s.key)):
        if s.group != group:
            group = s.group
            print(f"\n[{s.group_name}]")
        unit = f" ({s.unit})" if s.unit else ""
        print(f"  {s.key:<20} {s.name}{unit}  <{s.source}>")
    for name, msg in load_errors.items():
        print("\n" + tf("plugin {name} non caricato:", name=name) + f"\n{msg}", file=sys.stderr)
    return 0


def _pair_list(text: str) -> list[int]:
    """'1,3-5' → [1, 3, 4, 5]."""
    out: list[int] = []
    for part in filter(None, (p.strip() for p in text.split(","))):
        lo, _, hi = part.partition("-")
        out += list(range(int(lo), int(hi or lo) + 1))
    return out


def cmd_qc(args) -> int:
    from pistak import analyse_file
    from pistak.io import qc

    try:
        a = analyse_file(args.file, _limits(args), args.block)
        charges = tuple(int(x) for x in args.charges.split(","))
        mults = tuple(int(x) for x in args.mult.split(","))
        if len(charges) != 2 or len(mults) != 2:
            raise ValueError("--charges/--mult: due valori, A,B")
        options = qc.QCOptions(args.program, args.method, args.functional, args.dispersion, args.basis or "",
                               args.sapt_level, not args.no_cp, not args.raw_h, charges, mults, args.nprocs, args.mem)
        options.check()
        if args.pair:
            pairs = [a.pairs[n - 1] for n in _pair_list(args.pair) if 1 <= n <= len(a.pairs)]
        else:
            pairs = qc.exportable(a.pairs, None if args.kind == "all" else args.kind)
        if not pairs:
            print(tr("nessuna coppia da esportare"), file=sys.stderr)
            return 1
        if len(pairs) == 1 and not args.output:
            sys.stdout.write(qc.dimer_input(a, pairs[0], options))
            return 0
        out = Path(args.output or ".")
        # cartella se ci sono più coppie, se esiste già come cartella o se il nome non ha estensione (o finisce con /)
        if len(pairs) > 1 or out.is_dir() or not out.suffix or str(args.output).endswith(("/", "\\")):
            out.mkdir(parents=True, exist_ok=True)
        for p in pairs:
            print(tr("scritto"), qc.write_dimer_input(a, p, out, options))
    except Exception as e:  # noqa: BLE001
        print(tf("errore: {msg}", msg=describe_error(e)), file=sys.stderr)
        return 1
    return 0


def cmd_config(_args) -> int:
    from pistak.core.config import DEFAULTS_FILE, user_config_dir, user_config_file
    from pistak.params import plugin_dirs

    f = user_config_file()
    print(tf("cartella delle impostazioni: {path}", path=user_config_dir()))
    print(tf("impostazioni personali:      {path}", path=f) + ("" if f.exists() else "  " + tr("(non ancora creato)")))
    print(tf("valori predefiniti:          {path}", path=DEFAULTS_FILE))
    print(tr("cartelle dei plugin:"))
    for d in plugin_dirs():
        print(f"  {d}" + ("" if d.is_dir() else "  " + tr("(non esiste)")))
    return 0


def cmd_demo(_args) -> int:
    from pistak.core.demo import demo_cif

    sys.stdout.write(demo_cif())
    return 0


def cmd_gui(args) -> int:
    from pistak.gui.app import run

    return run(args.file, screenshot=getattr(args, "screenshot", None))


def _extract_lang(argv: list[str]) -> tuple[list[str], str | None]:
    """Toglie `--lang X` (o `--lang=X`) dagli argomenti: la lingua serve prima di costruire l'help."""
    out: list[str] = []
    lang = None
    it = iter(argv)
    for a in it:
        if a == "--lang":
            lang = next(it, None)
        elif a.startswith("--lang="):
            lang = a.split("=", 1)[1]
        else:
            out.append(a)
    return out, lang


def main(argv: list[str] | None = None) -> int:
    argv, lang = _extract_lang(list(sys.argv[1:] if argv is None else argv))
    if lang:
        i18n.set_language(lang)
    commands = {"analyze", "batch", "compare", "stats", "params", "demo", "config", "qc", "gui"}
    if not argv or argv[0] not in commands | {"-h", "--help"}:
        argv = ["gui", *argv]  # `pistak file.cif` apre la GUI
    ap = argparse.ArgumentParser(prog="pistak", description=tr("Analisi dei π-stack da file CIF."),
                                 epilog=tr("Lingua: --lang it|en, oppure la variabile PISTAK_LANG."))
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("analyze", help=tr("analizza un file CIF e stampa/esporta le coppie"))
    p.add_argument("file")
    p.add_argument("--rows", type=int, default=15, help=tr("righe da stampare (15)"))
    p.add_argument("--csv")
    p.add_argument("--json")
    p.add_argument("--xyz", nargs=2, metavar=("N", "FILE"), help=tr("scrive la coppia N in formato XYZ"))
    p.add_argument("--block", help=tr("blocco data_ da analizzare (predefinito: il primo con atomi)"))
    p.add_argument("--mc", type=int, default=0, metavar="N",
                   help=tr("incertezze dei parametri con un Monte Carlo di N campioni dalle s.u. del CIF (es. 200)"))
    p.add_argument("--mc-seed", type=int, default=0, help=tr("seme del generatore casuale del Monte Carlo"))
    p.add_argument("--mc-no-cell", action="store_true", help=tr("non perturbare i parametri di cella"))
    _add_limit_args(p)
    p.set_defaults(fn=cmd_analyze)

    p = sub.add_parser("batch", help=tr("analizza tutti i .cif di una cartella"))
    p.add_argument("folder")
    p.add_argument("--csv", default="riepilogo.csv")
    p.add_argument("--workers", type=int, default=0, help=tr("processi in parallelo (0 = automatico)"))
    _add_limit_args(p)
    p.set_defaults(fn=cmd_batch)

    p = sub.add_parser("compare", help=tr("confronta il motivo di stack di più strutture"))
    p.add_argument("inputs", nargs="+", help=tr("file .cif o cartelle"))
    p.add_argument("--ref", type=int, default=1, help=tr("struttura di riferimento (numero, 1 = la prima)"))
    p.add_argument("--csv")
    p.add_argument("--workers", type=int, default=0)
    _add_limit_args(p)
    p.set_defaults(fn=cmd_compare)

    p = sub.add_parser("stats", help=tr("statistiche sulle coppie di molti file"))
    p.add_argument("inputs", nargs="+", help=tr("file .cif o cartelle"))
    p.add_argument("--kind", choices=["all", "stack", "parallel", "tilted"], default="stack")
    p.add_argument("--keys", help=tr("parametri separati da virgola (predefiniti: i principali)"))
    p.add_argument("--unweighted", action="store_true", help=tr("ogni coppia conta una volta, senza le equivalenti per simmetria"))
    p.add_argument("--corr", action="store_true", help=tr("stampa la matrice di correlazione"))
    p.add_argument("--csv", help=tr("scrive tutte le coppie di tutti i file"))
    p.add_argument("--summary", help=tr("scrive la descrizione statistica"))
    p.add_argument("--workers", type=int, default=0)
    _add_limit_args(p)
    p.set_defaults(fn=cmd_stats)

    sub.add_parser("params", help=tr("elenca i parametri disponibili")).set_defaults(fn=cmd_params)
    sub.add_parser("demo", help=tr("scrive il CIF dimostrativo su stdout")).set_defaults(fn=cmd_demo)
    sub.add_parser("config", help=tr("mostra dove stanno impostazioni personali e plugin")).set_defaults(fn=cmd_config)

    p = sub.add_parser("qc", help=tr("il dimero di una coppia come input per ORCA, Psi4 o Gaussian"))
    p.add_argument("file")
    p.add_argument("--block")
    p.add_argument("--pair", help=tr("coppie da esportare, numeri della tabella (es. 1 oppure 1,3-5)"))
    p.add_argument("--kind", choices=["stack", "parallel", "tilted", "all"], default="stack",
                   help=tr("senza --pair: esporta tutte le coppie di questo tipo (stack)"))
    p.add_argument("--program", choices=["orca", "psi4", "gaussian"], default="orca")
    p.add_argument("--method", choices=["dft-d", "sapt"], default="dft-d", help=tr("dft-d con counterpoise, oppure sapt (solo psi4)"))
    p.add_argument("--functional", default="B3LYP")
    p.add_argument("--dispersion", choices=["d3bj", "d3zero", "d4", "none"], default="d3bj")
    p.add_argument("--basis", help=tr("base (predefinita: def2-TZVP per dft-d, jun-cc-pVDZ per sapt)"))
    p.add_argument("--sapt-level", default="sapt0", help=tr("livello SAPT di Psi4: sapt0, ssapt0, sapt2, sapt2+, sapt2+(3)…"))
    p.add_argument("--no-cp", action="store_true", help=tr("senza correzione counterpoise"))
    p.add_argument("--raw-h", action="store_true", help=tr("lascia le distanze X–H del CIF (non normalizzate)"))
    p.add_argument("--charges", default="0,0", help=tr("cariche dei monomeri A,B"))
    p.add_argument("--mult", default="1,1", help=tr("molteplicità dei monomeri A,B"))
    p.add_argument("--nprocs", type=int, default=4)
    p.add_argument("--mem", type=int, default=4, help=tr("memoria in GB"))
    p.add_argument("-o", "--output", help=tr("file (una coppia) o cartella (più coppie); senza: su stdout"))
    _add_limit_args(p)
    p.set_defaults(fn=cmd_qc)

    p = sub.add_parser("gui", help=tr("apre la finestra grafica"))
    p.add_argument("file", nargs="?")
    p.add_argument("--screenshot", metavar="PNG", help=tr("salva una cattura della finestra ed esce (per prove automatiche)"))
    p.set_defaults(fn=cmd_gui)

    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())

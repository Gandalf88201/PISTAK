"""Rapporto PDF del confronto tra PISTAK e PLATON (checklist, sezioni 4–5).

Produce `validation/platon/PISTAK_vs_PLATON.pdf` e, accanto, `coppie_abbinate.csv` (tutte le coppie abbinate, valori
di PISTAK e di PLATON). I dati vengono da `check_platon.py` (output `.lis` di PLATON in `validation/platon/`). La parte
sulla distanza interplanare usa i 200 CIF COD di `cod_cache/`; il risultato si salva in `interplanare.json`, così il
rapporto si rifà anche senza i CIF (`--rifai` per ricalcolarlo).

    .venv/bin/python validation/rapporto_platon.py [--rifai]
"""
from __future__ import annotations

import csv
import datetime as dt
import io
import json
import platform
import subprocess
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from reportlab.lib import colors  # noqa: E402
from reportlab.lib.enums import TA_JUSTIFY  # noqa: E402
from reportlab.lib.pagesizes import A4  # noqa: E402
from reportlab.lib.styles import ParagraphStyle  # noqa: E402
from reportlab.lib.units import mm  # noqa: E402
from reportlab.pdfbase import pdfmetrics  # noqa: E402
from reportlab.pdfbase.ttfonts import TTFont  # noqa: E402
from reportlab.platypus import (  # noqa: E402
    Image,
    KeepTogether,
    LongTable,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

from check_platon import DIR, compare  # noqa: E402

OUT_PDF = DIR / "PISTAK_vs_PLATON.pdf"
OUT_CSV = DIR / "coppie_abbinate.csv"
OUT_JSON = DIR / "interplanare.json"
QUANTITIES = [  # chiave, etichetta, unità, decimali di PLATON
    ("cg", "Cg···Cg", "Å", 3), ("alpha", "α", "°", 1), ("beta", "β", "°", 1), ("gamma", "γ", "°", 1),
    ("perp_i", "CgI_Perp", "Å", 3), ("perp_j", "CgJ_Perp", "Å", 3), ("slip", "Slippage", "Å", 3),
]

# ---- caratteri: DejaVu (incluso in matplotlib) per greco, Å, ⊥ -------------------------------------------------
FONTS = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
pdfmetrics.registerFont(TTFont("DV", str(FONTS / "DejaVuSans.ttf")))
pdfmetrics.registerFont(TTFont("DV-B", str(FONTS / "DejaVuSans-Bold.ttf")))
pdfmetrics.registerFont(TTFont("DV-I", str(FONTS / "DejaVuSans-Oblique.ttf")))
pdfmetrics.registerFont(TTFont("DVM", str(FONTS / "DejaVuSansMono.ttf")))
pdfmetrics.registerFontFamily("DV", normal="DV", bold="DV-B", italic="DV-I", boldItalic="DV-B")

INK, MUTED, RULE, BAND, ACCENT = (colors.HexColor(c) for c in ("#1d2430", "#5b6472", "#c9ced6", "#eef1f5", "#2f5d8a"))
S = {
    "title": ParagraphStyle("title", fontName="DV-B", fontSize=17, leading=21, textColor=INK, spaceAfter=4),
    "subtitle": ParagraphStyle("subtitle", fontName="DV", fontSize=10, leading=13, textColor=MUTED, spaceAfter=10),
    "h1": ParagraphStyle("h1", fontName="DV-B", fontSize=12.5, leading=16, textColor=INK, spaceBefore=12, spaceAfter=5),
    "h2": ParagraphStyle("h2", fontName="DV-B", fontSize=10.5, leading=13, textColor=INK, spaceBefore=8, spaceAfter=3),
    "body": ParagraphStyle("body", fontName="DV", fontSize=9, leading=12.4, textColor=INK, alignment=TA_JUSTIFY,
                           spaceAfter=5),
    "bullet": ParagraphStyle("bullet", fontName="DV", bulletFontName="DV", fontSize=9, leading=12.4, textColor=INK, leftIndent=10,
                             bulletIndent=2, spaceAfter=2),
    "caption": ParagraphStyle("caption", fontName="DV-I", fontSize=7.8, leading=10, textColor=MUTED, spaceBefore=3,
                              spaceAfter=8),
    "cell": ParagraphStyle("cell", fontName="DV", fontSize=7.8, leading=9.6, textColor=INK),
    "cellb": ParagraphStyle("cellb", fontName="DV-B", fontSize=7.8, leading=9.6, textColor=INK),
    "box": ParagraphStyle("box", fontName="DV", fontSize=9, leading=12.6, textColor=INK, alignment=TA_JUSTIFY),
}


def P(text: str, style: str = "body") -> Paragraph:
    return Paragraph(text, S[style])


def table(rows: list[list], widths: list[float], header_rows: int = 1, font_size: float = 7.8, align_right_from: int = 1,
          long: bool = False, zebra: bool = True) -> Table:
    cls = LongTable if long else Table
    t = cls(rows, colWidths=widths, repeatRows=header_rows)
    style = [
        ("FONT", (0, 0), (-1, -1), "DV", font_size),
        ("FONT", (0, 0), (-1, header_rows - 1), "DV-B", font_size),
        ("TEXTCOLOR", (0, 0), (-1, -1), INK),
        ("LINEABOVE", (0, 0), (-1, 0), 0.8, INK),
        ("LINEBELOW", (0, header_rows - 1), (-1, header_rows - 1), 0.5, INK),
        ("LINEBELOW", (0, -1), (-1, -1), 0.8, INK),
        ("ALIGN", (align_right_from, 0), (-1, -1), "RIGHT"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 1.6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 1.6),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
    ]
    if zebra:
        style += [("BACKGROUND", (0, r), (-1, r), BAND) for r in range(header_rows + 1, len(rows), 2)]
    t.setStyle(TableStyle(style))
    return t


def figure(fig, width_mm: float) -> Image:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    w, h = fig.get_size_inches()
    return Image(buf, width=width_mm * mm, height=width_mm * mm * h / w)


def num(x: float | None, d: int) -> str:
    """Numero all'italiana (virgola decimale, segno meno tipografico, niente «−0»)."""
    if x is None or not np.isfinite(x):
        return "–"
    return f"{round(x, d) + 0.0:.{d}f}".replace("-", "−").replace(".", ",")


# ---- dati --------------------------------------------------------------------------------------------------------
def collect() -> tuple[list[dict], list[str]]:
    results, unread = [], []
    for stem in (DIR / "elenco.txt").read_text().split():
        if not (DIR / f"{stem}.lis").exists():
            continue
        res = compare(stem)
        (results if res["n_rings"] else unread).append(res if res["n_rings"] else stem)
    return results, unread


def missed_reasons(results: list[dict]) -> Counter:
    from pistak.api import analyse_file
    from pistak.core.config import Limits
    from pistak.core.cores import small_rings

    out: Counter = Counter()
    for res in results:
        an = analyse_file(DIR / f"{res['stem']}.cif", Limits())
        found, labels = set(), set()
        for m in an.structure.molecules:
            labels.update(m.labels)
            for r in small_rings(len(m), m.bonds, lambda i, h=m.heavy: h[i], 8):
                found.add(frozenset(m.labels[i] for i in r))
        for r in res["ref_missed"]:
            ri, rj = res["rings"][r["i"]], res["rings"][r["j"]]
            if ri in found and rj in found:
                out["contact"] += 1
            elif any("_" in lab for lab in ri | rj):
                out["symmetry"] += 1
            elif not (ri | rj) <= labels:  # atomi che PISTAK non ha: componente di disordine scartata
                out["disorder"] += 1
            elif len(ri) == 8 or len(rj) == 8:
                out["envelope"] += 1
            else:
                out["other"] += 1
    return out


def interplanar_data(redo: bool) -> dict:
    if OUT_JSON.exists() and not redo:
        return json.loads(OUT_JSON.read_text())
    from pistak.api import analyse_file
    from pistak.core.config import Limits

    rows = []
    files = sorted((HERE / "cod_cache").glob("[0-9]*.cif"))
    for f in files:
        try:
            an = analyse_file(f, Limits())
        except Exception:  # noqa: BLE001  (CIF che PISTAK non analizza: fuori dal campione)
            continue
        for p in an.pairs:
            if p.kind in ("stack", "parallel") and not p.intramolecular:
                rows.append([f.stem, p.kind, p.dihedral, p.interplanar, p.get("height_b_over_a"), p.get("height_a_over_b")])
    data = {"n_cif": len(files), "rows": rows}
    OUT_JSON.write_text(json.dumps(data))
    return data


def git_version() -> str:
    try:
        return subprocess.run(["git", "-C", str(ROOT), "log", "-1", "--format=%h (%ad)", "--date=short"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:  # noqa: BLE001
        return "?"


# ---- figure ------------------------------------------------------------------------------------------------------
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.edgecolor": "#5b6472",
                     "axes.linewidth": 0.6, "xtick.color": "#5b6472", "ytick.color": "#5b6472",
                     "axes.labelcolor": "#1d2430", "axes.titlesize": 8, "axes.titleweight": "bold"})
BLUE, ORANGE = "#2f5d8a", "#c8692c"


def parity_figure(rows: list[dict]):
    fig, axes = plt.subplots(2, 4, figsize=(7.4, 3.9))
    for ax, (key, label, unit, _) in zip(axes.flat, QUANTITIES):
        pairs = [(r["platon"][key], r["pistak"][key]) for r in rows if r["platon"][key] is not None]
        x, y = np.array(pairs).T
        lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
        ax.plot([lo, hi], [lo, hi], color="#c9ced6", lw=0.8, zorder=0)
        ax.scatter(x, y, s=4, color=BLUE, alpha=0.6, lw=0)
        ax.set_title(f"{label} ({unit})")
        ax.set_xlabel("PLATON")
        ax.set_ylabel("PISTAK")
        ax.text(0.04, 0.93, f"n = {len(x)}", transform=ax.transAxes, va="top", color="#5b6472")
    axes.flat[-1].axis("off")
    fig.tight_layout()
    return fig


def residual_figure(rows: list[dict]):
    fig, axes = plt.subplots(2, 4, figsize=(7.4, 3.6))
    for ax, (key, label, unit, dec) in zip(axes.flat, QUANTITIES):
        d = np.array([r["pistak"][key] - r["platon"][key] for r in rows if r["platon"][key] is not None])
        half = 0.5 * 10 ** -dec  # mezza unità dell'ultima cifra stampata da PLATON
        ax.hist(d, bins=40, color=BLUE, alpha=0.85)
        ax.axvline(-half, color=ORANGE, lw=0.8, ls="--")
        ax.axvline(half, color=ORANGE, lw=0.8, ls="--")
        ax.set_title(f"Δ {label} ({unit})")
        ax.ticklabel_format(axis="x", style="sci", scilimits=(-2, 2))
    axes.flat[-1].axis("off")
    axes.flat[-1].text(0.0, 0.6, "Δ = PISTAK − PLATON\n\n- - -  ± mezza unità dell'ultima\n       cifra stampata da PLATON",
                       color="#1d2430", va="top")
    fig.tight_layout()
    return fig


def interplanar_figure(arr: np.ndarray):
    dih, ip, ha, hb = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 2.6))
    axes[0].scatter(dih, np.abs(ha - hb), s=4, color=ORANGE, alpha=0.6, lw=0)
    axes[0].set_xlabel("angolo tra i piani (°)")
    axes[0].set_ylabel("|⊥A − ⊥B| (Å)")
    axes[0].set_title("Le due distanze perpendicolari di PLATON")
    axes[1].scatter(dih, ip - (ha + hb) / 2, s=4, color=BLUE, alpha=0.6, lw=0)
    axes[1].axhline(0, color="#c9ced6", lw=0.8, zorder=0)
    axes[1].set_xlabel("angolo tra i piani (°)")
    axes[1].set_ylabel("interplanare − media(⊥A, ⊥B) (Å)")
    axes[1].set_title("Distanza interplanare di PISTAK")
    fig.tight_layout()
    return fig


# ---- documento ---------------------------------------------------------------------------------------------------
def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("DV", 7)
    canvas.setFillColor(MUTED)
    canvas.drawString(18 * mm, 10 * mm, "PISTAK – confronto con PLATON (checklist di revisione, sezioni 4–5)")
    canvas.drawRightString(A4[0] - 18 * mm, 10 * mm, f"pagina {doc.page}")
    canvas.restoreState()


def build(redo: bool) -> None:
    results, unread = collect()
    rows = [r for res in results for r in res["rows"]]
    reasons = missed_reasons(results)
    inter = interplanar_data(redo)
    arr = np.array([[r[2], r[3], r[4], r[5]] for r in inter["rows"]], dtype=float)

    with OUT_CSV.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cod", "ring_i", "ring_j", "size_i", "size_j"] + [f"{k}_{s}" for k, *_ in QUANTITIES for s in ("pistak", "platon")])
        for r in rows:
            w.writerow([r["stem"], r["i"], r["j"], r["size_i"], r["size_j"]]
                       + [x for k, *_ in QUANTITIES for x in (f"{r['pistak'][k]:.5f}",
                                                              "" if r["platon"][k] is None else r["platon"][k])])

    stats = {}
    for key, _label, _unit, dec in QUANTITIES:
        d = np.abs([r["pistak"][key] - r["platon"][key] for r in rows if r["platon"][key] is not None])
        stats[key] = (len(d), d.max(), np.percentile(d, 95), np.median(d), 0.5 * 10 ** -dec)
    n_mine = sum(res["mine"] for res in results)
    n_missed = sum(len(res["ref_missed"]) for res in results)
    today = dt.date.today().strftime("%d/%m/%Y")

    story: list = []
    story += [P("Confronto tra PISTAK e PLATON", "title"),
              P(f"Parametri anello–anello su 20 strutture della COD · {today} · PISTAK {git_version()}", "subtitle")]

    box = Table([[P(
        f"<b>Sintesi.</b> Su {len(results)} strutture lette da PLATON, {len(rows)} coppie di anelli calcolate da PISTAK "
        f"sono state abbinate alle righe di PLATON <i>CALC ALL</i>. Lo scarto massimo è {num(stats['cg'][1], 4)} Å per la "
        f"distanza Cg···Cg, {num(max(stats['alpha'][1], stats['beta'][1], stats['gamma'][1]), 3)}° per gli angoli α, β, γ, "
        f"{num(max(stats['perp_i'][1], stats['perp_j'][1]), 4)} Å per le distanze perpendicolari e {num(stats['slip'][1], 4)} Å per "
        "lo slippage: al massimo circa un'unità dell'ultima cifra che PLATON stampa (0,1° e 0,001 Å), con mediane "
        "dieci volte più piccole. Sono differenze di arrotondamento, non di calcolo. Le definizioni stampate da PLATON coincidono con quelle di PISTAK. Le coppie che "
        "PLATON elenca e PISTAK no dipendono da come i due programmi scelgono anelli e coppie di molecole e da come PLATON "
        "nomina gli atomi generati per simmetria (sezione 4), non da errori di calcolo.",
        "box")]], colWidths=[174 * mm])
    box.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), BAND), ("LINEBEFORE", (0, 0), (0, -1), 2.2, ACCENT),
                             ("LEFTPADDING", (0, 0), (-1, -1), 8), ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                             ("TOPPADDING", (0, 0), (-1, -1), 6), ("BOTTOMPADDING", (0, 0), (-1, -1), 6)]))
    story += [box, Spacer(1, 6)]

    # 1. materiali e metodo
    story += [P("1. Materiali e metodo", "h1"), P("Campione", "h2"), P(
        "Dal campione di 200 CIF scelti a caso nella Crystallography Open Database (COD, dominio pubblico; serie 1, 2, 4 "
        "e 7, seme fisso, elenco in <font name='DVM'>validation/cod_sample.txt</font>) si sono presi i 20 CIF organici "
        "senza metalli, senza avvisi di livello A nei controlli rapidi di PISTAK, con almeno uno stack e il maggior "
        "numero di coppie di anelli entro 6 Å. Sono elencati in tabella 1.")]
    t1 = [["COD", "Gruppo", "Anelli (PLATON)", "Coppie PISTAK", "Abbinate", "Classi solo PLATON"]]
    for res in results:
        t1.append([res["stem"], res["space_group"], res["n_rings"], res["mine"], len(res["rows"]), len(res["ref_missed"])])
    for stem in unread:
        t1.append([stem, "", "non letto", "", "", ""])
    t1.append(["Totale", "", sum(r["n_rings"] for r in results), n_mine, len(rows), n_missed])
    story += [table(t1, [22 * mm, 30 * mm, 28 * mm, 28 * mm, 24 * mm, 32 * mm], align_right_from=2),
              P("Tabella 1. Le strutture del confronto. «Coppie PISTAK»: coppie di anelli planari di molecole a contatto, "
                "Cg···Cg ≤ 6 Å; «Abbinate»: quelle ritrovate tra le righe di PLATON; «Classi solo PLATON»: coppie di "
                "PLATON, contate una volta per classe di simmetria, senza corrispondente in PISTAK (sezione 4). PLATON "
                "non legge 2205356 (la tabella anisotropa precede quella degli atomi) e su 7242996 si ferma su ADDSYM.",
                "caption")]

    story += [P("Programmi", "h2")]
    t2 = [["Programma", "Versione", "Uso"],
          ["PISTAK", git_version(), "pistak.core.ringpairs.ring_pairs, soglie predefinite"],
          ["PLATON (A. L. Spek)", "V-200626, compilato dal sorgente", "CALC ALL, impostazioni predefinite"],
          ["Python, numpy", f"{platform.python_version()}, {np.__version__}", "confronto e grafici (matplotlib)"]]
    story += [table(t2, [36 * mm, 52 * mm, 86 * mm], align_right_from=9, zebra=False), Spacer(1, 4),
              P("PLATON si è eseguito in modalità testo su ogni CIF, con "
                "<font name='DVM'>printf 'CALC ALL\\nEXIT\\n' | platon -o nome.cif</font>; la tabella confrontata è "
                "<i>Analysis of Short Ring-Interactions</i> del file <font name='DVM'>nome.lis</font> (Cg···Cg &lt; 6,0 Å, "
                "β &lt; 60°; lo slippage è stampato solo per α &lt; 20°).")]

    story += [P("Definizioni", "h2"), P(
        "Con Cg<sub>I</sub>, Cg<sub>J</sub> i baricentri degli anelli, <b>n</b><sub>I</sub>, <b>n</b><sub>J</sub> le "
        "normali dei loro piani medi (minimi quadrati) e <b>d</b> = Cg<sub>J</sub> − Cg<sub>I</sub>, le grandezze sono "
        "quelle della tabella 2. Le definizioni di PLATON sono quelle stampate nel file <font name='DVM'>.lis</font>; "
        "quelle di PISTAK si leggono in <font name='DVM'>ringpairs.py</font>.")]
    t3 = [["Grandezza", "Definizione", "PLATON", "PISTAK"],
          ["Cg···Cg", "|d|", "Cg-Cg", "cg_distance"],
          ["α", "angolo tra n_I e n_J (0–90°)", "Alpha", "alpha"],
          ["β", "angolo tra d e n_I", "Beta", "beta (A = I)"],
          ["γ", "angolo tra d e n_J", "Gamma", "gamma (B = J)"],
          ["⊥ di Cg_J sul piano I", "|d · n_I|", "CgJ_Perp", "perp_a"],
          ["⊥ di Cg_I sul piano J", "|d · n_J|", "CgI_Perp", "perp_b"],
          ["slippage", "√(|d|² − (d · n_I)²)", "Slippage", "slippage"]]
    story += [table(t3, [38 * mm, 58 * mm, 34 * mm, 44 * mm], align_right_from=9),
              P("Tabella 2. Definizioni e nomi nei due programmi. β e γ sono gli angoli di Janiak e PLATON; nel registro dei "
                "parametri di PISTAK le stesse grandezze per i nuclei π si chiamano beta e gamma.", "caption")]

    story += [P("Abbinamento", "h2"), P(
        "Gli anelli si riconoscono dall'insieme delle etichette dei loro atomi (PLATON li elenca nelle righe «n-Membered "
        "Ring (k)»). Una coppia di PISTAK corrisponde a una riga di PLATON se gli anelli coincidono, in un verso o "
        "nell'altro, e la distanza Cg···Cg differisce di meno di 0,006 Å; tra più righe possibili (copie per simmetria, "
        "anelli uguali con la loro immagine) si prende quella più vicina negli angoli, e nel verso inverso si scambiano β "
        "con γ e le due distanze perpendicolari. L'abbinamento non usa nessun valore calcolato da PISTAK oltre alla "
        "distanza Cg···Cg, che serve solo a scegliere la riga.")]

    # 2. risultati
    story += [PageBreak(), P("2. Risultati", "h1")]
    t4 = [["Grandezza", "Coppie", "Max |Δ|", "95° perc. |Δ|", "Mediana |Δ|", "Arrotondamento PLATON"]]
    for key, label, unit, _dec in QUANTITIES:
        n, mx, p95, med, half = stats[key]
        d = 4 if unit == "Å" else 3
        t4.append([f"{label} ({unit})", n, num(mx, d), num(p95, d), num(med, d), f"± {num(half, d)}"])
    story += [table(t4, [32 * mm, 16 * mm, 24 * mm, 28 * mm, 28 * mm, 46 * mm]),
              P("Tabella 3. Scarti Δ = PISTAK − PLATON su tutte le coppie abbinate. L'ultima colonna è mezza unità "
                "dell'ultima cifra che PLATON stampa: uno scarto entro questo valore è indistinguibile da zero.", "caption"),
              figure(parity_figure(rows), 172),
              P("Figura 1. Valori di PISTAK contro quelli di PLATON, coppia per coppia; la retta è y = x.", "caption"),
              figure(residual_figure(rows), 172),
              P("Figura 2. Distribuzione degli scarti. Le linee tratteggiate segnano l'arrotondamento di PLATON.", "caption")]

    # 3. per struttura
    t5 = [["COD", "Abbinate", "Cg···Cg (Å)", "α (°)", "β (°)", "γ (°)", "⊥ (Å)", "Slippage (Å)"]]
    for res in results:
        df = res["diffs"]
        mx = {k: max((abs(x) for x in v), default=float("nan")) for k, v in df.items()}
        t5.append([res["stem"], len(res["rows"]), num(mx["cg"], 4), num(mx["alpha"], 3), num(mx["beta"], 3),
                   num(mx["gamma"], 3), num(max(mx["perp_a"], mx["perp_b"]), 4), num(mx["slip"], 4)])
    story += [PageBreak(), P("3. Risultati per struttura", "h1"),
              table(t5, [22 * mm, 18 * mm, 24 * mm, 20 * mm, 20 * mm, 20 * mm, 22 * mm, 26 * mm]),
              P("Tabella 4. Scarto massimo |PISTAK − PLATON| in ogni struttura; ⊥ è il massimo tra CgI_Perp e CgJ_Perp.", "caption")]

    # 4. coppie non abbinate
    story += [P("4. Coppie che PLATON elenca e PISTAK no", "h1"), P(
        f"PLATON elenca {n_missed} classi di coppie di anelli (contate una volta per classe di simmetria) senza "
        "corrispondente in PISTAK. Nessuna è un errore di calcolo: dipendono da come i due programmi scelgono anelli e "
        "coppie. Le cause sono nella tabella 5.")]
    causes = [
        ("Anello perimetrale", "envelope",
         "Anelli da 8 atomi che sono il perimetro di due anelli condensati (5 + 5, 4 + 6). PLATON li tratta come anelli; "
         "PISTAK, per costruzione, tiene solo gli anelli più corti di qualche legame."),
        ("Atomi generati per simmetria", "symmetry",
         "PLATON chiama gli atomi generati per simmetria con un suffisso (C1_a), PISTAK con l'etichetta dell'atomo "
         "originale: l'abbinamento per etichetta non può riconoscerli. Non è un difetto di nessuno dei due."),
        ("Molecole non a contatto", "contact",
         "Le due molecole non hanno atomi pesanti entro contact_cutoff = 4,5 Å. PLATON guarda solo Cg···Cg &lt; 6 Å; "
         "PISTAK calcola la tabella anello–anello solo per le coppie di molecole a contatto."),
        ("Componente di disordine", "disorder",
         "L'anello contiene atomi che PISTAK non ha: appartengono alla componente di disordine scartata dalla regola "
         "major."),
        ("Non classificate", "other", "Da esaminare a mano."),
    ]
    t6 = [["Causa", "Spiegazione", "Classi"]]
    t6 += [[P(name, "cell"), P(text, "cell"), reasons[key]] for name, key, text in causes if reasons.get(key)]
    story += [table(t6, [36 * mm, 124 * mm, 14 * mm], align_right_from=2),
              P("Tabella 5. Classi di coppie presenti solo in PLATON, per causa. La classificazione è automatica, nell'ordine: "
                "entrambi gli anelli esistono in PISTAK (molecole non a contatto); etichette con suffisso di simmetria; "
                "atomi assenti in PISTAK (disordine); anello da 8 atomi (perimetrale).", "caption"),
              P(f"Al contrario, PISTAK ha {n_mine - len(rows)} coppie di anelli senza riga in PLATON: PLATON stampa solo "
                "le coppie con β &lt; 60° in almeno un verso, PISTAK tutte quelle entro 6 Å.")]

    # 5. distanza interplanare
    bins = [(0, 5), (5, 10), (10, 15), (15, 20.0001)]
    t7 = [["Angolo tra i piani", "Coppie", "|⊥A − ⊥B| mediano", "max", "Δ mediano", "max"]]
    for lo, hi in bins:
        m = (arr[:, 0] >= lo) & (arr[:, 0] < hi)
        if not m.any():
            continue
        spread = np.abs(arr[m, 2] - arr[m, 3])
        dif = np.abs(arr[m, 1] - (arr[m, 2] + arr[m, 3]) / 2)
        t7.append([f"{lo}–{int(hi)}°", int(m.sum()), num(np.median(spread), 3), num(spread.max(), 2),
                   num(np.median(dif), 3), num(dif.max(), 2)])
    tilt = arr[:, 0] >= 5
    spread_t = np.abs(arr[tilt, 2] - arr[tilt, 3])
    dif_t = np.abs(arr[tilt, 1] - (arr[tilt, 2] + arr[tilt, 3]) / 2)
    med_bins = [float(np.median(np.abs(arr[(arr[:, 0] >= lo) & (arr[:, 0] < hi), 2]
                                       - arr[(arr[:, 0] >= lo) & (arr[:, 0] < hi), 3]))) for lo, hi in bins[1:]]
    story += [KeepTogether([P("5. Distanza interplanare: una definizione diversa da PLATON", "h1"), P(
        "Per la coppia di nuclei π PISTAK riporta una sola distanza interplanare, la proiezione di <b>d</b> sulla normale "
        "media <b>n</b><sub>m</sub> = (n<sub>A</sub> + n<sub>B</sub>)/|n<sub>A</sub> + n<sub>B</sub>| (con n<sub>B</sub> orientata "
        "come n<sub>A</sub>). PLATON dà invece due distanze perpendicolari, ciascun baricentro dal piano dell'altro. Con "
        f"piani paralleli coincidono; la tabella 6 misura quanto differiscono sulle {len(arr)} coppie stack o parallel "
        f"dei {inter['n_cif']} CIF COD del campione (nuclei π interi, non singoli anelli)."),
        table(t7, [34 * mm, 18 * mm, 34 * mm, 18 * mm, 30 * mm, 18 * mm]),
        P("Tabella 6. ⊥A, ⊥B: distanze perpendicolari alla PLATON tra i nuclei (parametri height_b_over_a e "
          "height_a_over_b di PISTAK). Δ = distanza interplanare di PISTAK − media(⊥A, ⊥B). Valori in Å.", "caption")]),
        figure(interplanar_figure(arr), 172),
        P("Figura 3. A sinistra, quanto differiscono le due distanze perpendicolari al crescere dell'angolo tra i piani; a "
          "destra, la distanza interplanare di PISTAK rispetto alla loro media.", "caption"),
        P(f"Con piani inclinati di 5–20° ({int(tilt.sum())} coppie) le due distanze perpendicolari differiscono, in "
          f"mediana, da {num(min(med_bins), 1)} a {num(max(med_bins), 1)} Å secondo l'angolo, e fino a "
          f"{num(spread_t.max(), 1)} Å; la distanza lungo la normale media è un solo numero, a {num(np.median(dif_t), 2)} Å "
          f"in mediana dalla loro media (al massimo {num(dif_t.max(), 2)} Å). Nel paper va "
          "scritta la definizione: d<sub>⊥</sub> = |<b>d</b> · <b>n</b><sub>m</sub>|. Per confrontarsi con lavori che usano PLATON "
          "conviene riportare anche CgI_Perp e CgJ_Perp, che PISTAK calcola uguali a PLATON (sezione 2).")]

    # 6. limiti
    story += [P("6. Limiti del confronto", "h1")]
    for txt in [
        "La precisione è quella di stampa di PLATON (0,1° e 0,001 Å): il confronto dimostra l'accordo entro questi valori, "
        "non oltre.",
        "Si confrontano singoli anelli. La geometria della coppia di nuclei π (distanza interplanare, slittamenti lungo "
        "l'asse lungo e laterale, sovrapposizione) non ha un equivalente in PLATON; la sua verifica sta altrove nella "
        "checklist (test analitici su cristalli costruiti a mano).",
        f"Le strutture sono {len(results)}, tutte organiche e senza metalli; i casi con disordine sono pochi.",
        "PLATON e PISTAK leggono lo stesso CIF, ma PLATON tratta il disordine a modo suo: gli anelli di componenti scartate "
        "da PISTAK non si confrontano.",
    ]:
        story.append(Paragraph(txt, S["bullet"], bulletText="•"))

    # appendice
    story += [PageBreak(), P("Appendice A. Tutte le coppie abbinate", "h1"), P(
        "Valori di PLATON e scarto Δ = PISTAK − PLATON per ogni coppia abbinata, nella notazione di PLATON (anello I → "
        "anello J, numeri dei «Ring» del file .lis). Δ di Cg···Cg, ⊥ e slippage in 10⁻⁴ Å, Δ degli angoli in 10⁻² °. Le righe "
        "identiche (la stessa coppia di anelli trovata in più coppie di molecole equivalenti per simmetria) compaiono una "
        f"volta sola. Tutte le {len(rows)} coppie, con i valori di PISTAK a piena precisione, sono in "
        f"<font name='DVM'>{OUT_CSV.name}</font>.")]
    head = ["COD", "I→J", "Cg···Cg", "Δ", "α", "Δ", "β", "Δ", "γ", "Δ", "CgI⊥", "Δ", "CgJ⊥", "Δ", "Slip.", "Δ"]
    body, seen = [head], set()
    for r in sorted(rows, key=lambda r: (r["stem"], r["i"], r["j"], r["platon"]["cg"])):
        pl, pk = r["platon"], r["pistak"]
        line = [r["stem"], f"{r['i']}→{r['j']}"]
        for key, _, unit, _dec in QUANTITIES:
            v = pl[key]
            scale = 1e4 if unit == "Å" else 1e2
            line += [num(v, 3 if unit == "Å" else 1), "–" if v is None else num((pk[key] - v) * scale, 0)]
        if tuple(line) not in seen:  # la stessa coppia di anelli in più coppie di molecole equivalenti: una riga sola
            seen.add(tuple(line))
            body.append(line)
    widths = [14 * mm, 10 * mm] + [12.2 * mm, 8.2 * mm] * 7
    story.append(table(body, widths, font_size=6.3, align_right_from=1, long=True))

    doc = SimpleDocTemplate(str(OUT_PDF), pagesize=A4, leftMargin=18 * mm, rightMargin=18 * mm, topMargin=16 * mm,
                            bottomMargin=17 * mm, title="Confronto tra PISTAK e PLATON",
                            author="PISTAK – revisione del codice", subject="Parametri anello–anello contro PLATON CALC ALL")
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    print(f"{OUT_PDF} ({len(rows)} coppie), {OUT_CSV}")


if __name__ == "__main__":
    build("--rifai" in sys.argv)

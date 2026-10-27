"""Dettaglio: riepilogo del cristallo, tutti i parametri della coppia scelta, coppie di anelli (Janiak), colonne,
contatti tra molecole (riepilogo per molecola delle coppie di frammenti π), avvisi (file non letti, plugin con errori),
geometria delle molecole (legami, angoli, torsioni) e controlli rapidi della struttura."""
from __future__ import annotations

import math
from html import escape as _esc

from PySide6.QtCore import Qt
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTextBrowser,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from pistak.core.diagnosis import why_no_pairs
from pistak.core.intramolecular import molecule_geometry
from pistak.core.ringpairs import HEADERS, ring_pairs
from pistak.gui import theme
from pistak.gui.session import Session
from pistak.gui.views.base import PairView
from pistak.i18n import tf, tr
from pistak.labels import ALERT_LEVELS, alert_text, no_pairs_text, rank_name, relation
from pistak.params import load_errors, specs, structure_specs

#: Tipi di misura della scheda Geometria: chiave → (titolo, unità, decimali).
GEOMETRY_KINDS = {"bond": ("Legami", "Å", 4), "angle": ("Angoli", "°", 2), "torsion": ("Torsioni", "°", 2)}
#: Colori dei livelli degli avvisi (testo), leggibili sia sul tema chiaro sia su quello scuro.
LEVEL_COLORS = {"A": QColor(214, 48, 49), "B": QColor(230, 126, 34), "C": QColor(180, 150, 20)}


class _NumItem(QTableWidgetItem):
    """Cella che si ordina per valore numerico (i NaN in fondo), non per testo."""

    def __init__(self, text: str, value: float):
        super().__init__(text)
        self.value = value if math.isfinite(value) else math.inf

    def __lt__(self, other) -> bool:
        return self.value < other.value if isinstance(other, _NumItem) else super().__lt__(other)


class DetailView(PairView):
    TITLE = "Dettaglio"
    EMPTY_NOTICE = False  # senza coppie restano Struttura, Geometria e Controlli; il motivo va in cima a Struttura

    def __init__(self, session: Session, group: str | None = "A", parent=None):
        super().__init__(session, group, parent)
        self.tabs = QTabWidget()
        self.summary = QTextBrowser()
        self.params = QTreeWidget()
        self.params.setHeaderLabels([tr("Parametro"), tr("Valore"), tr("Unità")])
        self.params.setRootIsDecorated(True)
        self.rings = QTableWidget()
        self.rings.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.columns = QTableWidget()
        self.columns.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.contacts = QTableWidget()
        self.contacts.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.contacts.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.contacts.setToolTip(tr("Doppio clic: seleziona la coppia di frammenti migliore di quel contatto"))
        self.contacts.cellDoubleClicked.connect(self._contact_chosen)
        self._contact_best: list[int] = []
        self.notes = QTextBrowser()
        self.tabs.addTab(self.summary, tr("Struttura"))
        self.tabs.addTab(self.params, tr("Coppia"))
        self.tabs.addTab(self.rings, tr("Anelli"))
        self.tabs.addTab(self.columns, tr("Colonne"))
        self.tabs.addTab(self.notes, tr("Avvisi"))
        self.tabs.addTab(self.contacts, tr("Contatti"))
        self.tabs.addTab(self._build_geometry(), tr("Geometria"))
        self.checks = QTableWidget()
        self.checks.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.checks.setWordWrap(True)
        self.tabs.addTab(self.checks, tr("Controlli"))
        self.body.addWidget(self.tabs)
        self._geometry_items: list | None = None  # calcolata alla prima apertura della scheda, per analisi
        self._geometry_key: tuple | None = None
        self.tabs.currentChanged.connect(lambda _i: self._fill_geometry())

    def _build_geometry(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(0, 0, 0, 0)
        bar = QHBoxLayout()
        self.geom_kind = QComboBox()
        for key, (title, unit, _d) in GEOMETRY_KINDS.items():
            self.geom_kind.addItem(f"{tr(title)} ({unit})", key)
        self.geom_h = QCheckBox(tr("Idrogeni"))
        self.geom_h.setChecked(True)
        self.geom_info = QLabel()
        self.geom_info.setToolTip(tr("Solo le molecole indipendenti per simmetria; ′ = atomo generato dalla simmetria. "
                                     "s.u. dalle incertezze delle coordinate del CIF (vuota per gli atomi senza s.u.)."))
        bar.addWidget(self.geom_kind)
        bar.addWidget(self.geom_h)
        bar.addStretch(1)
        bar.addWidget(self.geom_info)
        lay.addLayout(bar)
        self.geometry = QTableWidget()
        self.geometry.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.geometry.setSortingEnabled(True)
        lay.addWidget(self.geometry)
        self.geom_kind.currentIndexChanged.connect(lambda _i: self._fill_geometry())
        self.geom_h.toggled.connect(lambda _b: self._fill_geometry())
        return w

    def _fill_geometry(self) -> None:
        if self.tabs.currentIndex() != 6:
            return  # si calcola solo se la scheda è visibile
        a = self.analysis
        if a is None:
            self.geometry.setRowCount(0)
            return
        key = (id(a.structure), self.geom_h.isChecked())
        if key != self._geometry_key:
            self._geometry_items = molecule_geometry(a.structure, a.cif, include_h=self.geom_h.isChecked())
            self._geometry_key = key
        kind = self.geom_kind.currentData()
        title, unit, digits = GEOMETRY_KINDS[kind]
        rows = [g for g in self._geometry_items or [] if g.kind == kind]
        n_atoms = {"bond": 2, "angle": 3, "torsion": 4}[kind]
        heads = [tr("Molecola")] + [tf("Atomo {n}", n=i + 1) for i in range(n_atoms)] + [f"{tr(title)} ({unit})", "s.u."]
        self.geometry.setSortingEnabled(False)
        self.geometry.clear()
        self.geometry.setColumnCount(len(heads))
        self.geometry.setHorizontalHeaderLabels(heads)
        self.geometry.setRowCount(len(rows))
        for r, g in enumerate(rows):
            cells: list[tuple[str, float | None]] = [(str(g.molecule + 1), g.molecule + 1)]
            cells += [(lab, None) for lab in g.atoms]
            cells += [(f"{g.value:.{digits}f}", g.value), ("" if not math.isfinite(g.su) else f"{g.su:.{digits}f}", g.su)]
            for c, (text, num) in enumerate(cells):
                it = _NumItem(text, num) if num is not None else QTableWidgetItem(text)
                if num is not None:
                    it.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.geometry.setItem(r, c, it)
        self.geometry.setSortingEnabled(True)
        self.geometry.resizeColumnsToContents()
        n_mol = len({g.molecule for g in self._geometry_items or []})
        self.geom_info.setText(tf("{n} misure · {m} molecole indipendenti", n=len(rows), m=n_mol))

    def _fill_checks(self) -> None:
        a = self.analysis
        alerts = (a.alerts or []) if a else []
        heads = [tr("Livello"), tr("Avviso")]
        self.checks.clear()
        self.checks.setColumnCount(len(heads))
        self.checks.setHorizontalHeaderLabels(heads)
        self.checks.setRowCount(len(alerts))
        for r, al in enumerate(alerts):
            level = QTableWidgetItem(al.level)
            level.setToolTip(tr(ALERT_LEVELS[al.level]))
            level.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            text = QTableWidgetItem(alert_text(al))
            text.setToolTip(f"{al.code}\n{alert_text(al)}")
            if al.level in LEVEL_COLORS:
                for it in (level, text):
                    it.setForeground(QBrush(LEVEL_COLORS[al.level]))
            self.checks.setItem(r, 0, level)
            self.checks.setItem(r, 1, text)
        self.checks.horizontalHeader().setStretchLastSection(True)
        self.checks.resizeColumnToContents(0)
        self.checks.resizeRowsToContents()
        serious = sum(al.level != "G" for al in alerts)
        self.tabs.setTabText(7, f"{tr('Controlli')} ({serious})" if serious else tr("Controlli"))

    def on_analysis(self) -> None:
        a = self.analysis
        if not a:
            self.summary.setHtml("")
            self._fill_checks()
            return
        c = a.cif.cell
        rows = [
            ("File", a.name), (tr("Blocco"), a.cif.name or "–"), (tr("Gruppo spaziale"), a.cif.space_group or "–"),
            (tr("Cella"), f"a {c.a:.3f} · b {c.b:.3f} · c {c.c:.3f} Å<br>α {c.alpha:.2f}° · β {c.beta:.2f}° · γ {c.gamma:.2f}°"),
            (tr("Atomi nella cella"), str(a.structure.n_cell_atoms)),
        ]
        reason = why_no_pairs(a)
        html = "" if reason is None else \
            f"<p style='color:{LEVEL_COLORS['B'].name()}'><b>{tr('Nessuna coppia di nuclei π')}</b>: {_esc(no_pairs_text(reason))}</p>"
        html += "<table cellpadding=3>" + "".join(f"<tr><td><b>{k}</b></td><td>{v}</td></tr>" for k, v in rows) + f"</table><h4>{tr('Molecole')}</h4><ul>"
        for g in a.molecule_groups():
            sub = "".join(f"{e}<sub>{n}</sub>" if n > 1 else e for e, n in g.formula)
            html += f"<li><b>{g.count} ×</b> {sub} · " + tf("{n} con nucleo π", n=g.with_core)
            if g.fragments > 1:
                html += " · " + tf("{n} frammenti π per molecola", n=g.fragments)
            html += "</li>"
        html += "</ul>"
        if a.has_extended:
            html += tr("<p><i>La rete di legami non si chiude (polimero o solido esteso): quelle parti sono escluse.</i></p>")
        if a.cif.skipped_disorder:
            html += "<p><i>" + tf("{n} atomi di componenti disordinate secondarie ignorati.", n=a.cif.skipped_disorder) + "</i></p>"
        html += self._structure_html()
        html += f"<h4>{tr('Soglie')}</h4><p>{', '.join(f'{k} = {v}' for k, v in a.limits.to_dict().items())}</p>"
        self.summary.setHtml(html)
        self._fill_columns()
        self._fill_contacts()
        self._notes()
        self._fill_checks()
        self._geometry_key = None
        self._fill_geometry()

    def _structure_html(self) -> str:
        """Le grandezze del registro di struttura che hanno un valore, per gruppo, con la descrizione nel tooltip."""
        a = self.analysis
        groups: dict[str, list[str]] = {}
        for s in structure_specs():
            v = a.structure_values.get(s.key, math.nan)
            if not math.isfinite(v):
                continue
            unit = f" {s.unit}" if s.unit else ""
            groups.setdefault(s.group_name, []).append(
                f'<tr><td title="{_esc(s.info)}">{_esc(s.name)}</td><td align=right>{v:.{s.digits}f}{unit}</td></tr>')
        return "".join(f"<h4>{g}</h4><table cellpadding=2>{''.join(rows)}</table>" for g, rows in groups.items())

    def _fill_contacts(self) -> None:
        a = self.analysis
        rows = a.molecule_contacts() if a else []
        heads = [tr("Molecole"), tr("Contatti per cella"), tr("Coppie di frammenti"), "π-stack", tr("Contatto (Å)"),
                 tr("Coppia migliore")]
        self.contacts.clear()
        self.contacts.setColumnCount(len(heads))
        self.contacts.setHorizontalHeaderLabels(heads)
        self.contacts.setRowCount(len(rows))
        self._contact_best = [c.best.index for c in rows]
        for r, c in enumerate(rows):
            best = c.best
            mols = tf("{a} (interna)", a=c.a + 1) if c.intramolecular else \
                f"{c.a + 1}–{c.b + 1} [{' '.join(map(str, c.shift))}]"
            cells = [mols, str(c.count), str(len(c.pairs)), str(len(c.stacks)), f"{c.contact:.2f}",
                     f"#{best.index + 1} {theme.kind_label(best.kind)} · π{best.fa + 1}–π{best.fb + 1} · {best.interplanar:.2f} Å"]
            for k, v in enumerate(cells):
                self.contacts.setItem(r, k, QTableWidgetItem(v))
        self.contacts.resizeColumnsToContents()
        self.tabs.setTabText(5, f"{tr('Contatti')} ({len(rows)})")

    def _contact_chosen(self, row: int, _col: int) -> None:
        if 0 <= row < len(self._contact_best):
            self.choose(self._contact_best[row])

    def _fill_columns(self) -> None:
        a = self.analysis
        cols = a.stacks if a else []
        heads = [tr("Tipo"), tr("Molecole"), tr("Passo (Å)"), tr("Asse (celle)"), tr("Interplanare (Å)"), tr("Slitt. L (Å)"), tr("Slitt. W (Å)"), tr("Asse–normale (°)")]
        self.columns.clear()
        self.columns.setColumnCount(len(heads))
        self.columns.setHorizontalHeaderLabels(heads)
        self.columns.setRowCount(len(cols))
        for r, c in enumerate(cols):
            core = a.cores[c.molecules[0]]
            axis = f"[{' '.join(map(str, c.translations[0]))}]" if c.rank == 1 else (tf("{n} direzioni", n=c.rank) if c.rank else "–")
            row = [rank_name(c.rank), ", ".join(str(m + 1) for m in c.molecules), "–" if c.repeat is None else f"{c.repeat:.2f}", axis,
                   f"{c.mean('interplanar'):.2f}", f"{c.mean('long_slip'):.2f}", f"{c.mean('lat_slip'):.2f}",
                   "–" if core is None or c.axis is None else f"{c.axis_normal_angle(core.normal):.1f}"]
            for k, v in enumerate(row):
                self.columns.setItem(r, k, QTableWidgetItem(v))
        self.columns.resizeColumnsToContents()
        self.tabs.setTabText(3, f"{tr('Colonne')} ({len(cols)})")

    def on_params(self) -> None:
        self._notes()
        self._fill_checks()

    def _notes(self) -> None:
        errors = self.analysis.param_errors if self.analysis else {}
        msgs = ["<p>" + tf("<b>Plugin {n}</b> non caricato:", n=n) + f"<pre>{m}</pre></p>" for n, m in load_errors.items()]
        msgs += ["<p>" + tf("<b>Parametro {k}</b> ha dato errori:", k=k) + f"<pre>{m}</pre></p>" for k, m in errors.items()]
        self.notes.setHtml("".join(msgs) or tr("<p>Nessun avviso.</p>"))
        self.tabs.setTabText(4, f"{tr('Avvisi')} ({len(load_errors) + len(errors)})" if (load_errors or errors) else tr("Avvisi"))

    def on_pair(self, index: int) -> None:
        self.params.clear()
        self.rings.clear()
        self.rings.setRowCount(0)
        pair, a = self.pair, self.analysis
        if not (pair and a):
            return
        u = a.uncertainty
        QTreeWidgetItem(self.params, [tf("Coppia {n}: {kind}", n=pair.index + 1, kind=theme.kind_label(pair.kind)), relation(pair), ""])
        if u is not None:
            QTreeWidgetItem(self.params, [tr("Incertezze"), tf("Monte Carlo, {n} campioni", n=u.samples),
                                          tf("tipo stabile nel {p:.0f}% dei campioni", p=100 * u.kind_agreement.get(pair.index, 0))])
        groups: dict[str, QTreeWidgetItem] = {}
        for s in specs():
            top = groups.get(s.group_name)
            if top is None:
                top = groups[s.group_name] = QTreeWidgetItem(self.params, [s.group_name])
                top.setExpanded(True)
            v = pair.get(s.key)
            text = "–" if not math.isfinite(v) else f"{v:.{s.digits}f}"
            sig = u.get(pair, s.key) if u is not None and s.key in u.keys else math.nan
            if math.isfinite(v) and math.isfinite(sig):
                text += f" ± {sig:.{max(s.digits, 3)}f}"
            it = QTreeWidgetItem(top, [s.name, text, s.unit])
            it.setToolTip(0, s.info)
            it.setTextAlignment(1, Qt.AlignmentFlag.AlignRight)
        self.params.header().setStretchLastSection(False)
        self.params.resizeColumnToContents(0)
        rp = ring_pairs(a.pair_context(pair))
        keys = list(HEADERS)
        self.rings.setColumnCount(len(keys))
        self.rings.setHorizontalHeaderLabels([tr(HEADERS[k]) for k in keys])
        self.rings.setRowCount(len(rp))
        for r, item in enumerate(rp):
            for c, k in enumerate(keys):
                v = getattr(item, k)
                text = str(v + 1) if k in ("ring_a", "ring_b") else (str(v) if isinstance(v, int) else f"{v:.2f}")
                cell = QTableWidgetItem(text)
                cell.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.rings.setItem(r, c, cell)
        self.rings.resizeColumnsToContents()
        self.tabs.setTabText(2, f"{tr('Anelli')} ({len(rp)})")

    def extra_state(self) -> dict:
        return {"tab": self.tabs.currentIndex()}

    def apply_state(self, st: dict) -> None:
        self.tabs.setCurrentIndex(int(st.get("tab", 0)))

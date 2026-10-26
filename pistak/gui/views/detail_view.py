"""Dettaglio: riepilogo del cristallo, tutti i parametri della coppia scelta, coppie di anelli (Janiak)
e avvisi (file non letti, plugin con errori)."""
from __future__ import annotations

import math

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QTableWidget, QTableWidgetItem, QTabWidget, QTextBrowser, QTreeWidget, QTreeWidgetItem

from pistak.core.ringpairs import HEADERS, ring_pairs
from pistak.gui import theme
from pistak.gui.session import Session
from pistak.gui.views.base import PairView
from pistak.i18n import tf, tr
from pistak.params import errors, load_errors, specs


class DetailView(PairView):
    TITLE = "Dettaglio"

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
        self.notes = QTextBrowser()
        self.tabs.addTab(self.summary, tr("Struttura"))
        self.tabs.addTab(self.params, tr("Coppia"))
        self.tabs.addTab(self.rings, tr("Anelli"))
        self.tabs.addTab(self.columns, tr("Colonne"))
        self.tabs.addTab(self.notes, tr("Avvisi"))
        self.body.addWidget(self.tabs)

    def on_analysis(self) -> None:
        a = self.analysis
        if not a:
            self.summary.setHtml("")
            return
        c = a.cif.cell
        rows = [
            ("File", a.name), (tr("Blocco"), a.cif.name or "–"), (tr("Gruppo spaziale"), a.cif.space_group or "–"),
            (tr("Cella"), f"a {c.a:.3f} · b {c.b:.3f} · c {c.c:.3f} Å<br>α {c.alpha:.2f}° · β {c.beta:.2f}° · γ {c.gamma:.2f}°"),
            (tr("Atomi nella cella"), str(a.structure.n_cell_atoms)),
        ]
        html = "<table cellpadding=3>" + "".join(f"<tr><td><b>{k}</b></td><td>{v}</td></tr>" for k, v in rows) + f"</table><h4>{tr('Molecole')}</h4><ul>"
        for g in a.molecule_groups():
            sub = "".join(f"{e}<sub>{n}</sub>" if n > 1 else e for e, n in g.formula)
            html += f"<li><b>{g.count} ×</b> {sub} · " + tf("{n} con nucleo π", n=g.with_core) + "</li>"
        html += "</ul>"
        if a.has_extended:
            html += tr("<p><i>La rete di legami non si chiude (polimero o solido esteso): quelle parti sono escluse.</i></p>")
        if a.cif.skipped_disorder:
            html += "<p><i>" + tf("{n} atomi di componenti disordinate secondarie ignorati.", n=a.cif.skipped_disorder) + "</i></p>"
        html += f"<h4>{tr('Soglie')}</h4><p>{', '.join(f'{k} = {v}' for k, v in a.limits.to_dict().items())}</p>"
        self.summary.setHtml(html)
        self._fill_columns()
        self._notes()

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
            row = [c.kind, ", ".join(str(m + 1) for m in c.molecules), "–" if c.repeat is None else f"{c.repeat:.2f}", axis,
                   f"{c.mean('interplanar'):.2f}", f"{c.mean('long_slip'):.2f}", f"{c.mean('lat_slip'):.2f}",
                   "–" if core is None or c.axis is None else f"{c.axis_normal_angle(core.normal):.1f}"]
            for k, v in enumerate(row):
                self.columns.setItem(r, k, QTableWidgetItem(v))
        self.columns.resizeColumnsToContents()
        self.tabs.setTabText(3, f"{tr('Colonne')} ({len(cols)})")

    def on_params(self) -> None:
        self._notes()

    def _notes(self) -> None:
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
        QTreeWidgetItem(self.params, [tf("Coppia {n}: {kind}", n=pair.index + 1, kind=theme.kind_label(pair.kind)), pair.relation(), ""])
        groups: dict[str, QTreeWidgetItem] = {}
        for s in specs():
            top = groups.get(s.group_name)
            if top is None:
                top = groups[s.group_name] = QTreeWidgetItem(self.params, [s.group_name])
                top.setExpanded(True)
            v = pair.get(s.key)
            it = QTreeWidgetItem(top, [s.name, "–" if not math.isfinite(v) else f"{v:.{s.digits}f}", s.unit])
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

"""Tabella delle coppie: ordinabile, con colonne a scelta e filtro per tipo."""
from __future__ import annotations

import math

from PySide6.QtCore import QAbstractTableModel, QModelIndex, QSortFilterProxyModel, Qt
from PySide6.QtGui import QBrush
from PySide6.QtWidgets import QCheckBox, QHeaderView, QMenu, QTableView, QToolButton

from pistak.gui.session import Session
from pistak.gui.theme import KIND_COLORS, KIND_LABELS, kind_label, qcolor
from pistak.gui.views.base import PairView
from pistak.i18n import tr
from pistak.labels import relation
from pistak.params import get_spec, specs

DEFAULT_COLUMNS = ["interplanar", "long_slip", "lat_slip", "dihedral", "centroid_distance", "contact", "overlap", "count"]
FIXED = ["#", "Tipo", "Relazione"]


class PairModel(QAbstractTableModel):
    def __init__(self) -> None:
        super().__init__()
        self.pairs: list = []
        self.keys: list[str] = list(DEFAULT_COLUMNS)

    def set_data(self, pairs: list, keys: list[str]) -> None:
        self.beginResetModel()
        self.pairs, self.keys = pairs, keys
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return 0 if parent.isValid() else len(self.pairs)

    def columnCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return 0 if parent.isValid() else len(FIXED) + len(self.keys)

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):  # noqa: N802
        if orientation != Qt.Orientation.Horizontal:
            return None
        if role == Qt.ItemDataRole.DisplayRole:
            if section < len(FIXED):
                return tr(FIXED[section])
            spec = get_spec(self.keys[section - len(FIXED)])
            return spec.title if spec else self.keys[section - len(FIXED)]
        if role == Qt.ItemDataRole.ToolTipRole and section >= len(FIXED):
            spec = get_spec(self.keys[section - len(FIXED)])
            return spec.info or None
        return None

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        p = self.pairs[index.row()]
        c = index.column()
        if c < len(FIXED):
            if role == Qt.ItemDataRole.DisplayRole:
                return (p.index + 1, kind_label(p.kind), relation(p))[c]
            if role == Qt.ItemDataRole.UserRole:
                return (p.index, ["stack", "parallel", "tilted"].index(p.kind), relation(p))[c]
            if role == Qt.ItemDataRole.ForegroundRole and c == 1:
                return QBrush(qcolor(KIND_COLORS[p.kind]))
            return None
        key = self.keys[c - len(FIXED)]
        v = p.get(key)
        if role == Qt.ItemDataRole.DisplayRole:
            spec = get_spec(key)
            return "–" if not math.isfinite(v) else f"{v:.{spec.digits if spec else 2}f}"
        if role == Qt.ItemDataRole.UserRole:
            return v if math.isfinite(v) else -1e30
        if role == Qt.ItemDataRole.TextAlignmentRole:
            return int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        return None


class SubsetProxy(QSortFilterProxyModel):
    """Ordina per il valore numerico e mostra solo le coppie del sottoinsieme (e dei tipi) scelti."""

    def __init__(self) -> None:
        super().__init__()
        self.setSortRole(Qt.ItemDataRole.UserRole)
        self.subset: set[int] | None = None
        self.hidden_kinds: set[str] = set()

    def filterAcceptsRow(self, row: int, parent) -> bool:  # noqa: N802
        p = self.sourceModel().pairs[row]
        return (self.subset is None or p.index in self.subset) and p.kind not in self.hidden_kinds

    def refresh(self) -> None:
        if hasattr(self, "beginFilterChange"):  # Qt ≥ 6.9: invalidateFilter è deprecato
            self.beginFilterChange()
            self.endFilterChange(QSortFilterProxyModel.Direction.Rows)
        else:
            self.invalidateFilter()


class TableView(PairView):
    TITLE = "Tabella"

    def __init__(self, session: Session, group: str | None = "A", parent=None):
        super().__init__(session, group, parent)
        self.model = PairModel()
        self.proxy = SubsetProxy()
        self.proxy.setSourceModel(self.model)
        self.table = QTableView()
        self.table.setModel(self.proxy)
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QTableView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QTableView.EditTrigger.NoEditTriggers)
        self.table.verticalHeader().hide()
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setAlternatingRowColors(True)
        self.table.selectionModel().currentRowChanged.connect(self._row_changed)
        self.body.addWidget(self.table)

        self.kind_boxes: dict[str, QCheckBox] = {}
        for kind, label in KIND_LABELS.items():
            cb = QCheckBox(tr(label))
            cb.setChecked(True)
            cb.toggled.connect(self._kinds_changed)
            self.kind_boxes[kind] = cb
            self.toolbar.addWidget(cb)
        self.cols_btn = QToolButton()
        self.cols_btn.setText(tr("Colonne ▾"))
        self.cols_btn.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.cols_btn.setMenu(QMenu(self.cols_btn))
        self.toolbar.addWidget(self.cols_btn)
        self.reset_btn = QToolButton()
        self.reset_btn.setText(tr("Tutte le coppie"))
        self.reset_btn.setToolTip(tr("Toglie il filtro impostato da un'altra vista (scatter, istogramma)"))
        self.reset_btn.clicked.connect(lambda: self.filter_to(None))
        self.reset_btn.setVisible(False)
        self.toolbar.addWidget(self.reset_btn)
        self._build_menu()
        self._updating = False

    # ---- colonne -------------------------------------------------------------------
    def _build_menu(self) -> None:
        menu = self.cols_btn.menu()
        menu.clear()
        groups: dict[str, list] = {}
        for s in specs():
            groups.setdefault(s.group_name, []).append(s)
        for g, items in groups.items():
            sub = menu.addMenu(g)
            for s in items:
                act = sub.addAction(s.title)
                act.setCheckable(True)
                act.setChecked(s.key in self.model.keys)
                act.setToolTip(s.info)
                act.toggled.connect(lambda on, k=s.key: self._toggle_column(k, on))

    def _toggle_column(self, key: str, on: bool) -> None:
        keys = [k for k in self.model.keys if k != key]
        if on:
            keys = [s.key for s in specs() if s.key in keys or s.key == key]
        self.model.set_data(self.model.pairs, keys)
        self._restore_selection()

    # ---- hook ----------------------------------------------------------------------
    def on_analysis(self) -> None:
        a = self.analysis
        keys = [k for k in self.model.keys if get_spec(k)]
        self.model.set_data(a.pairs if a else [], keys)
        self.proxy.subset = None
        self.proxy.refresh()
        self.table.sortByColumn(-1, Qt.SortOrder.AscendingOrder)

    def on_params(self) -> None:
        self._build_menu()
        self.model.keys = [k for k in self.model.keys if get_spec(k)]

    def on_pair(self, index: int) -> None:
        self._restore_selection()

    def on_subset(self, indices) -> None:
        self.proxy.subset = indices
        self.proxy.refresh()
        self.reset_btn.setVisible(indices is not None)
        self._restore_selection()

    def _kinds_changed(self) -> None:
        self.proxy.hidden_kinds = {k for k, cb in self.kind_boxes.items() if not cb.isChecked()}
        self.proxy.refresh()

    def _restore_selection(self) -> None:
        idx = self.index
        self._updating = True
        try:
            for r in range(self.proxy.rowCount()):
                if self.proxy.mapToSource(self.proxy.index(r, 0)).row() == idx:
                    self.table.selectRow(r)
                    self.table.scrollTo(self.proxy.index(r, 0))
                    return
            self.table.clearSelection()
        finally:
            self._updating = False

    def _row_changed(self, current: QModelIndex, _prev) -> None:
        if self._updating or not current.isValid():
            return
        self.choose(self.proxy.mapToSource(current).row())

    def extra_state(self) -> dict:
        return {"columns": list(self.model.keys), "hidden_kinds": sorted(self.proxy.hidden_kinds)}

    def apply_state(self, st: dict) -> None:
        keys = [k for k in st.get("columns", self.model.keys) if get_spec(k)]
        self.model.set_data(self.model.pairs, keys or list(DEFAULT_COLUMNS))
        self._build_menu()
        hidden = set(st.get("hidden_kinds", []))
        for kind, cb in self.kind_boxes.items():
            cb.blockSignals(True)
            cb.setChecked(kind not in hidden)
            cb.blockSignals(False)
        self._kinds_changed()

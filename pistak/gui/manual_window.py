"""Il manuale nella finestra (Aiuto ▸ Manuale, F1): indice dei capitoli a sinistra, testo a destra, ricerca.

Il manuale è Markdown (`pistak/manual/<lingua>.md`, con le figure in `img/<lingua>/`), lo stesso che si legge sul
repository. Qt importa il Markdown ma non crea àncore per i titoli: i link interni («#3-primo-esempio-guidato») si
risolvono qui, confrontando lo slug del link con quello di ogni titolo (le stesse regole di GitHub).
"""
from __future__ import annotations

import re
from pathlib import Path

from PySide6.QtCore import QEvent, QObject, Qt, QUrl
from PySide6.QtGui import QDesktopServices, QImageReader, QKeySequence, QShortcut, QTextImageFormat
from PySide6.QtWidgets import (
    QComboBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QSplitter,
    QTextBrowser,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from pistak import i18n
from pistak.gui.screen import fit_to_screen
from pistak.i18n import tr

MANUAL_DIR = Path(__file__).resolve().parents[1] / "manual"
#: Margine tra le figure e il bordo del testo (px): le figure larghe si rimpiccioliscono per starci.
IMAGE_MARGIN = 40


def manual_path(lang: str) -> Path:
    path = MANUAL_DIR / f"{lang}.md"
    return path if path.exists() else MANUAL_DIR / "it.md"


def slug(title: str) -> str:
    """L'àncora che GitHub dà a un titolo: minuscole, senza punteggiatura, spazi → trattini."""
    return re.sub(r"[^\w\- ]", "", title.strip().lower()).replace(" ", "-")


def headings(text: str) -> list[tuple[int, str]]:
    """(livello, titolo) dei titoli ## e ### del Markdown, fuori dai blocchi di codice."""
    out, code = [], False
    for line in text.splitlines():
        if line.startswith("```"):
            code = not code
        elif not code and (m := re.match(r"(#{2,3}) (.+)", line)):
            out.append((len(m.group(1)), m.group(2).strip()))
    return out


class ManualWindow(QMainWindow):
    def __init__(self, parent: QWidget | None = None, lang: str | None = None):
        super().__init__(parent)
        self.setWindowTitle(tr("Manuale di PISTAK"))
        self.lang = lang or i18n.get_language()
        self.toc = QTreeWidget()
        self.toc.setHeaderHidden(True)
        self.toc.itemClicked.connect(lambda item, _col: self.go_to(item.data(0, Qt.ItemDataRole.UserRole)))
        self.text = QTextBrowser()
        self.text.setOpenLinks(False)
        self.text.anchorClicked.connect(self._link)
        self.text.viewport().installEventFilter(self)  # le figure seguono la larghezza del testo
        self._fitted_width = 0
        self.search = QLineEdit()
        self.search.setPlaceholderText(tr("Cerca nel manuale (Invio: successivo)"))
        self.search.returnPressed.connect(self.find_next)
        self.lang_combo = QComboBox()
        for code, name in i18n.LANGS.items():
            self.lang_combo.addItem(name, code)
        self.lang_combo.setCurrentIndex(max(0, self.lang_combo.findData(self.lang)))
        self.lang_combo.currentIndexChanged.connect(lambda: self.load(self.lang_combo.currentData()))

        bar = QHBoxLayout()
        bar.addWidget(QLabel("🔍"))
        bar.addWidget(self.search, 1)
        bar.addWidget(self.lang_combo)
        split = QSplitter()
        split.addWidget(self.toc)
        split.addWidget(self.text)
        split.setStretchFactor(1, 1)
        split.setSizes([260, 760])
        root = QWidget()
        lay = QVBoxLayout(root)
        lay.addLayout(bar)
        lay.addWidget(split, 1)
        self.setCentralWidget(root)
        QShortcut(QKeySequence.StandardKey.Find, self, activated=self.search.setFocus)
        fit_to_screen(self, 1100, 820, min_width=420, min_height=320)
        self.load(self.lang)

    # ---- caricamento -------------------------------------------------------------------
    def load(self, lang: str) -> None:
        self.lang = lang
        path = self.path = manual_path(lang)
        self.source = path.read_text(encoding="utf-8")
        self.titles = headings(self.source)
        self.text.document().setBaseUrl(QUrl.fromLocalFile(str(path)))
        self.text.setSearchPaths([str(path.parent)])
        self.text.setMarkdown(self.source)
        self._fit_images()
        self.toc.clear()
        parent: QTreeWidgetItem | None = None
        for level, title in self.titles:
            item = QTreeWidgetItem([title])
            item.setData(0, Qt.ItemDataRole.UserRole, slug(title))
            if level == 2 or parent is None:
                self.toc.addTopLevelItem(item)
                parent = item
            else:
                parent.addChild(item)
        self.text.verticalScrollBar().setValue(0)

    def _fit_images(self) -> None:
        """Le figure più larghe del testo si rimpiccioliscono (in proporzione) invece di far scorrere di lato."""
        doc = self.text.document()
        width = max(200, self.text.viewport().width() - IMAGE_MARGIN)
        self._fitted_width = width
        block = doc.begin()
        while block.isValid():
            it = block.begin()
            while not it.atEnd():
                frag = it.fragment()
                fmt = frag.charFormat()
                if frag.isValid() and fmt.isImageFormat():
                    img_fmt = QTextImageFormat(fmt.toImageFormat())
                    size = QImageReader(str(self.path.parent / img_fmt.name())).size()  # solo l'intestazione del PNG
                    w0, h0 = size.width(), size.height()
                    if w0 > 0 and h0 > 0:
                        w = min(w0, width)
                        img_fmt.setWidth(w)
                        img_fmt.setHeight(h0 * w / w0)  # in proporzione
                        cursor = self.text.textCursor()
                        cursor.setPosition(frag.position())
                        cursor.setPosition(frag.position() + frag.length(), cursor.MoveMode.KeepAnchor)
                        cursor.setCharFormat(img_fmt)
                it += 1
            block = block.next()

    def eventFilter(self, obj: QObject, event: QEvent) -> bool:  # noqa: N802
        if (event.type() == QEvent.Type.Resize and obj is self.text.viewport() and hasattr(self, "source")
                and abs(self.text.viewport().width() - IMAGE_MARGIN - self._fitted_width) > 4):
            self._fit_images()
        return super().eventFilter(obj, event)

    # ---- navigazione -------------------------------------------------------------------
    def go_to(self, anchor: str) -> bool:
        """Scorre fino al titolo con quell'àncora. False se non c'è."""
        title = next((t for _lvl, t in self.titles if slug(t) == anchor), None)
        if title is None:
            return False
        block = self.text.document().begin()
        while block.isValid():
            if block.blockFormat().headingLevel() and block.text().strip() == title:
                top = self.text.document().documentLayout().blockBoundingRect(block).top()
                self.text.verticalScrollBar().setValue(int(top))
                return True
            block = block.next()
        return False

    def _link(self, url: QUrl) -> None:
        target = url.toString()
        if url.scheme() in ("http", "https", "mailto"):
            QDesktopServices.openUrl(url)
        elif target.startswith("#"):
            self.go_to(target[1:])
        elif target.endswith(".md"):  # l'altra lingua
            lang = Path(target).stem
            self.lang_combo.setCurrentIndex(max(0, self.lang_combo.findData(lang)))
        elif url.fragment():
            self.go_to(url.fragment())

    def find_next(self) -> None:
        text = self.search.text().strip()
        if text and not self.text.find(text):
            cursor = self.text.textCursor()  # da capo
            cursor.movePosition(cursor.MoveOperation.Start)
            self.text.setTextCursor(cursor)
            self.text.find(text)

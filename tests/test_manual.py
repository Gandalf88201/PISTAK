"""The manual (pistak/manual/it.md, en.md): figures present, internal links valid, the two languages aligned, the
library code examples actually run."""
import os
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MANUAL = ROOT / "pistak" / "manual"
LANGS = ("it", "en")


def text(lang: str) -> str:
    return (MANUAL / f"{lang}.md").read_text(encoding="utf-8")


def slug(title: str) -> str:
    return re.sub(r"[^\w\- ]", "", title.strip().lower()).replace(" ", "-")


def headings(md: str) -> list[tuple[int, str]]:
    out, code = [], False
    for line in md.splitlines():
        if line.startswith("```"):
            code = not code
        elif not code and (m := re.match(r"(#{2,3}) (.+)", line)):
            out.append((len(m.group(1)), m.group(2).strip()))
    return out


@pytest.mark.parametrize("lang", LANGS)
def test_figures_present(lang):
    images = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", text(lang))
    assert len(images) >= 12
    for rel in images:
        assert rel.startswith(f"img/{lang}/"), rel  # each language has its own figures
        assert (MANUAL / rel).is_file(), rel


@pytest.mark.parametrize("lang", LANGS)
def test_internal_links(lang):
    md = text(lang)
    anchors = {slug(t) for _lvl, t in headings(md)}
    links = re.findall(r"\]\(#([^)]+)\)", md)
    assert len(links) >= 20
    assert not [a for a in links if a not in anchors]


def test_both_languages_have_same_structure():
    it, en = text("it"), text("en")
    assert [lvl for lvl, _ in headings(it)] == [lvl for lvl, _ in headings(en)]
    figs = [re.findall(r"img/\w+/([\w.]+)", t) for t in (it, en)]
    assert figs[0] == figs[1]
    assert it.count("```") == en.count("```")


def test_window_finds_same_slug():
    pytest.importorskip("PySide6")  # the manual window belongs to the GUI: without PySide6 (computation-only CI) skip
    from pistak.gui.manual_window import headings as h2
    from pistak.gui.manual_window import slug as s2

    for lang in LANGS:
        assert [(lvl, s2(t)) for lvl, t in h2(text(lang))] == [(lvl, slug(t)) for lvl, t in headings(text(lang))]


@pytest.mark.parametrize("lang", LANGS)
def test_library_example(lang, tmp_path, monkeypatch, capsys):
    """The first Python block of the library chapter, run as it is (from the project folder)."""
    md = text(lang)
    section = md[md.index("## 18."):md.index("## 19.")]
    code = re.findall(r"```python\n(.*?)```", section, re.S)[0]
    code = code.replace('"examples/', f'"{ROOT}/examples/').replace("samples=200", "samples=20")
    monkeypatch.chdir(tmp_path)  # export.write_json writes a.json here
    exec(compile(code, f"manual-{lang}", "exec"), {})
    out = capsys.readouterr().out
    assert "stack" in out and "molecule dimer" in out and (tmp_path / "a.json").exists()


def test_figures_regenerable():
    """The figure script exists and uses only the CIFs in examples/ and the COD ones in tests/data/cod/ (public
    domain): CSD structures cannot be redistributed."""
    src = (ROOT / "tools" / "manual_figures.py").read_text(encoding="utf-8")
    others = re.sub(r"COD / \"\d+\.cif\"", "", src.replace("_anthracene.cif", "").replace("_naphthalene.cif", ""))
    assert "examples" in src and "DATI_REALI" not in src and not re.search(r"\.cif[\"']\)", others)
    for name in re.findall(r"COD / \"(\d+\.cif)\"", src):
        assert (ROOT / "tests" / "data" / "cod" / name).is_file(), name
    assert os.path.getsize(MANUAL / "img" / "it" / "main_window.png") < 200_000  # 256-colour PNG

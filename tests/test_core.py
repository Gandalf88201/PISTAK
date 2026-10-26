"""Porting dei test del sito (tests/pistak.test.ts) più verifiche sulle parti nuove."""
import numpy as np
import pytest

import cifgen
from pistak.core import analysis as core
from pistak.core.cif import Cell, CifError, cif_number, element_of, parse_cif, parse_symop
from pistak.core.config import Limits
from pistak.core.demo import DEMO_STACK, demo_cif
from pistak.core.elements import is_known
from pistak.core.geometry import clip_convex, convex_hull, fit_plane, polygon_area
from pistak.core.pairs import measure_pair
from pistak.core.structure import cell_vectors, to_cart, to_frac


class TestCif:
    def test_numeri_con_incertezza_e_mancanti(self):
        assert cif_number("10.3528(16)") == 10.3528
        assert cif_number("-0.15128(16)") == -0.15128
        assert np.isnan(cif_number("?"))

    def test_operazioni_di_simmetria(self):
        op = parse_symop("-x+1/2, y+1/2, -z+1/2")
        assert op.r.tolist() == [[-1, 0, 0], [0, 1, 0], [0, 0, -1]]
        assert op.t.tolist() == [0.5, 0.5, 0.5]
        assert parse_symop("1/2+X,1/2-Y,Z").t.tolist() == [0.5, 0.5, 0]
        assert parse_symop("x-y, x, z+0.25").r[0].tolist() == [1, -1, 0]
        with pytest.raises(CifError):
            parse_symop("x, y")

    def test_elemento_da_simbolo_o_etichetta(self):
        assert element_of("C", "C12A", is_known) == "C"
        assert element_of("Cl", "CL1", is_known) == "Cl"
        assert element_of(None, "Cl1", is_known) == "Cl"
        assert element_of(None, "C1A", is_known) == "C"
        assert element_of("O2-", "O1", is_known) == "O"

    def test_apici_campi_di_testo_e_disordine(self):
        cif = parse_cif("\n".join([
            "data_x", "_cell_length_a 5(1)", "_cell_length_b 6", "_cell_length_c 7",
            "_cell_angle_alpha 90", "_cell_angle_beta 90", "_cell_angle_gamma 90",
            "_publ_section_title", ";", "A title with 'quotes' and loop_ inside", ";",
            "_space_group_name_H-M_alt 'P 1'",
            "loop_", "_atom_site_label", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_disorder_group",
            "C1 0.1 0.2 0.3 .", "C2 0.2 0.2 0.3 1", "C3 0.3 0.2 0.3 2",
        ]))
        assert cif.cell.a == 5
        assert [a.label for a in cif.atoms] == ["C1", "C2"]
        assert cif.skipped_disorder == 1
        assert len(cif.symops) == 1

    def test_senza_operazioni_si_usano_quelle_del_gruppo(self):
        text = ("data_x\n_cell_length_a 5\n_cell_length_b 5\n_cell_length_c 5\n_cell_angle_alpha 90\n_cell_angle_beta 90\n"
                "_cell_angle_gamma 90\n_symmetry_space_group_name_H-M 'P 21/c'\nloop_\n_atom_site_label\n_atom_site_fract_x\n"
                "_atom_site_fract_y\n_atom_site_fract_z\nC1 0 0 0\n")
        c = parse_cif(text)
        assert c.symops_source == "table" and {o.text for o in c.symops} == {"x,y,z", "-x,-y,-z", "-x,y+1/2,-z+1/2", "x,-y+1/2,z+1/2"}

    def test_senza_operazioni_e_con_gruppo_sconosciuto_si_ferma(self):
        text = ("data_x\n_cell_length_a 5\n_cell_length_b 5\n_cell_length_c 5\n_cell_angle_alpha 90\n_cell_angle_beta 90\n"
                "_cell_angle_gamma 90\n_symmetry_space_group_name_H-M 'P 42/n c m'\nloop_\n_atom_site_label\n_atom_site_fract_x\n"
                "_atom_site_fract_y\n_atom_site_fract_z\nC1 0 0 0\n")
        with pytest.raises(CifError) as e:
            parse_cif(text)
        assert e.value.code == "noSymops"

    def test_stile_cif2_e_apice_interno(self):
        cif = parse_cif("\n".join([
            "data_y", "_cell_length_a 5", "_cell_length_b 6", "_cell_length_c 7",
            "_cell_angle_alpha 90", "_cell_angle_beta 90", "_cell_angle_gamma 90",
            "loop_", "_space_group_symop.operation_xyz", "'x, y, z'", "'-x, -y, -z'",
            "loop_", "_atom_site_label", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z", "O'Brien 0.1 0.2 0.3",
        ]))
        assert len(cif.symops) == 2
        assert cif.atoms[0].label == "O'Brien"

    def test_file_vuoto(self):
        with pytest.raises(CifError) as e:
            parse_cif("  \n")
        assert e.value.code == "empty"


class TestGeometria:
    def test_cella_frazionarie_e_cartesiane_si_corrispondono(self):
        m = cell_vectors(Cell(10.35, 13.03, 21.21, 96.3, 97.4, 98.4))
        assert np.linalg.norm(m[2]) == pytest.approx(21.21, abs=1e-6)
        assert to_frac(m, to_cart(m, np.array([0.1, -0.4, 0.7]))) == pytest.approx([0.1, -0.4, 0.7], abs=1e-10)

    def test_piano_medio_di_punti_complanari(self):
        p = fit_plane(np.array([[0, 0, 1], [4, 0, 1], [0, 1, 1], [4, 1, 1]], dtype=float))
        assert abs(p.normal[2]) == pytest.approx(1, abs=1e-9)
        assert abs(p.long[0]) == pytest.approx(1, abs=1e-9)
        assert p.rms == pytest.approx(0, abs=1e-9)

    def test_involucro_e_intersezione(self):
        sq = np.array([[0, 0], [2, 0], [2, 2], [0, 2], [1, 1]], dtype=float)
        h = convex_hull(sq)
        assert len(h) == 4 and polygon_area(h) == pytest.approx(4)
        other = h + np.array([1.0, 1.0])
        assert polygon_area(clip_convex(h, other)) == pytest.approx(1)
        assert len(clip_convex(h, h + 10)) == 0


class TestDemo:
    """Il cristallo dimostrativo del sito: la costruzione si ritrova nei valori misurati."""

    @staticmethod
    @pytest.fixture(scope="module")
    def result():
        cif = parse_cif(demo_cif())
        return cif, core.analyse_cif(cif, "demo", Limits())

    def test_il_centro_di_inversione_completa_la_molecola(self, result):
        cif, a = result
        assert len(cif.atoms) == 9  # 5 C + 4 H nell'unità asimmetrica
        assert a.structure.n_cell_atoms == 18
        assert len(a.structure.molecules) == 1
        assert len(a.structure.molecules[0].bonds) == 19  # 11 C–C + 8 C–H
        assert not a.structure.molecules[0].extended

    def test_nucleo_pi_sono_i_due_anelli_condensati(self, result):
        _, a = result
        assert a.cores[0].kind == "rings"
        assert len(a.cores[0].atoms) == 10

    def test_distanza_interplanare_e_slittamenti(self, result):
        _, a = result
        stack = next(p for p in a.pairs if p.kind == "stack")
        assert stack.interplanar == pytest.approx(DEMO_STACK["interplanar"], abs=1e-3)
        assert stack.long_slip == pytest.approx(DEMO_STACK["long_slip"], abs=1e-3)
        assert stack.lat_slip == pytest.approx(DEMO_STACK["lat_slip"], abs=1e-3)
        assert stack.dihedral == pytest.approx(0, abs=1e-3)
        assert stack.overlap > 0

    def test_il_vicino_nello_stesso_piano_non_e_uno_stack(self, result):
        _, a = result
        side = measure_pair(a.structure, a.cores, 0, 0, (0, 1, 0), a.limits)
        assert side.kind == "parallel"
        assert side.interplanar == pytest.approx(0, abs=1e-6)
        assert side.overlap == 0
        # i carboni sono a 7,5 − 2·1,4 = 4,7 Å: oltre la soglia dei contatti, quindi non è tra le coppie elencate
        assert all(p.kind == "stack" for p in a.pairs)

    def test_soglie_modificabili(self):
        a = core.analyse_text(demo_cif(), "demo", Limits(max_interplanar=3.0))
        assert all(p.kind != "stack" for p in a.pairs)  # 3,40 Å supera il nuovo limite


class TestCasiCostruiti:
    def test_catena_infinita_e_rete_estesa(self):
        a = core.analyse_text(cifgen.infinite_chain(), "chain", Limits())
        assert a.structure.molecules[0].extended
        assert a.cores == [None] and a.pairs == []

    def test_catena_senza_anelli_usa_la_molecola(self):
        a = core.analyse_text(cifgen.chain_molecule(), "zigzag", Limits())
        assert len(a.structure.molecules) == 2
        assert all(c is not None and c.kind == "molecule" for c in a.cores)

    def test_lisca_di_pesce_ha_coppie_inclinate(self):
        a = core.analyse_text(cifgen.anthracene_p21a_like(), "herringbone", Limits())
        assert len(a.structure.molecules) == 2
        assert {"stack", "tilted"} <= {p.kind for p in a.pairs}

    def test_atomi_in_cella_e_molecole_intere(self):
        a = core.analyse_text(cifgen.tilted_sandwich(), "tilted", Limits())
        n_per_mol = len(a.structure.molecules[0])
        assert n_per_mol == 24 and all(len(m) == n_per_mol for m in a.structure.molecules)
        assert a.structure.n_cell_atoms == len(a.structure.molecules) * n_per_mol

    def test_troppi_atomi(self):
        from pistak.core.structure import StructureError

        with pytest.raises(StructureError):
            core.analyse_text(cifgen.anthracene_p21a_like(), "x", Limits(max_cell_atoms=20))

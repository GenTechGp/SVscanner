"""Tests for src/rm_align_chain.py. Run: python -m unittest discover -s test -p 'test_rm_align_chain.py'

Fixtures in test/rm_align_chain/ are RepeatMasker 4.2.6 output (-e rmblast -species human,
Dfam 4.0, -a): sva_reproducer.* is test/sva_strand_reproducer.fa (one SVA, forward and
reverse complement); elements.* holds three elements cut from the PB_human run: an L1 split
over L1M4_5end/L1M4_orf2 (DEL.57 ID 34), an ERV joined from THE1D-int and its MSTA LTR
(DEL.96 ID 56) and an AluY with a (T)n tail of another class and strand (INS.119 ID 142).
"""
import io
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import unittest.mock
from collections import defaultdict

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(TEST_DIR, "..", "src")
FIXTURES = os.path.join(TEST_DIR, "rm_align_chain")
sys.path.insert(0, SRC_DIR)

import rm_align_chain as rc  # noqa: E402

REPRO_LEN = 2724  # both reproducer sequences
FORWARD = "hg38_chr2_41746002-41748725_forward"
REVCOMP = "hg38_chr2_41746002-41748725_revcomp"
WHOLE = (0, 10**9)  # SV region covering the whole query


def load(prefix):
    """(blocks by (query, ID), .out lines by (query, ID)) for one fixture."""
    blocks, outs = defaultdict(list), defaultdict(list)
    for b in rc.parse_align(os.path.join(FIXTURES, prefix + ".align")):
        blocks[(b["query"], b["id"])].append(b)
    for r in rc.parse_out(os.path.join(FIXTURES, prefix + ".out")):
        outs[(r["query"], r["id"])].append(r)
    return blocks, outs


def block(qbeg, qend, cbeg, strand="+", name="X", cls="LINE/L1", score=100, rm_id=1, query="q"):
    """A gap-free, identical synthetic block: query qbeg..qend on consensus cbeg..cbeg+len-1."""
    n = qend - qbeg + 1
    return {"score": score, "query": query, "qbeg": qbeg, "qend": qend, "strand": strand,
            "cbeg": cbeg, "cend": cbeg + n - 1, "name": name, "cls": cls, "id": rm_id,
            "q_row": "A" * n, "c_row": "A" * n}


def chain(blocks, sv=WHOLE, jump_cost=rc.DEFAULT_JUMP_COST):
    return rc.chain_element(blocks, [], sv[0], sv[1], jump_cost)


class ParseRepeatMaskerOutput(unittest.TestCase):
    def test_align_blocks_carry_out_ids(self):
        blocks, outs = load("sva_reproducer.fa")
        self.assertEqual(sum(len(v) for v in blocks.values()), 36)
        # Every .align ID is an .out ID; the VNTR tiles stay under one ID per strand.
        self.assertEqual(set(blocks), set(outs))
        self.assertEqual(len(blocks[(FORWARD, 2)]), 17)
        self.assertEqual(len(blocks[(REVCOMP, 3)]), 18)

    def test_minus_strand_rows(self):
        # Consensus rows of minus-strand blocks start with "C "; coordinates run downwards.
        blocks, _ = load("sva_reproducer.fa")
        b = max(blocks[(REVCOMP, 3)], key=lambda b: b["score"])
        self.assertEqual(b["strand"], "-")
        self.assertEqual(len(b["q_row"]), len(b["c_row"]))
        mapped = rc.base_map(b, 0)
        self.assertEqual(mapped[b["qbeg"]][0], b["cend"])

    def test_out_without_hits(self):
        with tempfile.NamedTemporaryFile("w", suffix=".out") as f:
            f.write("\nThere were no repetitive sequences detected in x.fa\n")
            f.flush()
            self.assertEqual(rc.parse_out(f.name), [])


class ChainRealElements(unittest.TestCase):
    def test_sva_reproducer_both_strands(self):
        # One expanded SVA_F: every consensus position matched once, the extra VNTR copies
        # inserted as one block at the same consensus position on both strands.
        blocks, outs = load("sva_reproducer.fa")
        fwd = rc.chain_element(blocks[(FORWARD, 2)], outs[(FORWARD, 2)], *WHOLE, rc.DEFAULT_JUMP_COST)
        rev = rc.chain_element(blocks[(REVCOMP, 3)], outs[(REVCOMP, 3)], *WHOLE, rc.DEFAULT_JUMP_COST)
        for res in (fwd, rev):
            self.assertEqual(res["repeat"], "SVA_F")
            self.assertEqual(res["consensus_len"], 1375)
            self.assertEqual(res["matched"], 1375)
            self.assertEqual(res["consensus_hit"], 1375)
            self.assertEqual(len(res["ins_blocks"]), 1)
        self.assertEqual(fwd["ins_blocks"], [(540, 1778, 480)])
        # The reverse complement's insertion is the same stretch, mirrored.
        s, e, c = rev["ins_blocks"][0]
        self.assertEqual((REPRO_LEN - e + 1, REPRO_LEN - s + 1, c), fwd["ins_blocks"][0])

    def test_l1_split_models_on_one_frame(self):
        # L1M4_5end and L1M4_orf2 are numbered separately; .out puts them on L1M4's frame.
        blocks, outs = load("elements.fa")
        bl = blocks[("DEL.57", 34)]
        self.assertEqual({b["name"] for b in bl}, {"L1M4_5end", "L1M4_orf2"})
        maps = [rc.base_map(b, off) for b, off in rc.frame_offsets(bl, outs[("DEL.57", 34)], False)]
        path = rc.chain_path(maps, False, rc.DEFAULT_JUMP_COST)
        self.assertEqual((path[0][1], path[-1][1]), (1389, 5483))  # the .out span
        res = rc.chain_element(bl, outs[("DEL.57", 34)], *WHOLE, rc.DEFAULT_JUMP_COST)
        self.assertEqual(res["matched"], aligned_bases(bl))  # nothing dropped
        self.assertEqual(res["repeat"], "L1M4")

    def test_erv_ltr_and_int(self):
        # .out names this element MSTA throughout (THE1D-int as "MSTA 1-907"); the -int and
        # LTR blocks are placed one after the other instead of on top of each other.
        blocks, outs = load("elements.fa")
        bl = blocks[("DEL.96", 56)]
        self.assertEqual({b["name"] for b in bl}, {"THE1D-int", "MSTA"})
        res = rc.chain_element(bl, outs[("DEL.96", 56)], *WHOLE, rc.DEFAULT_JUMP_COST)
        self.assertEqual(res["matched"], aligned_bases(bl))

    def test_tail_of_another_class_left_out(self):
        blocks, outs = load("elements.fa")
        bl = blocks[("INS.119", 142)]
        self.assertEqual({b["cls"] for b in bl}, {"SINE/Alu", "Simple_repeat"})
        res = rc.chain_element(bl, outs[("INS.119", 142)], *WHOLE, rc.DEFAULT_JUMP_COST)
        self.assertEqual(res["n_left_out"], 1)
        self.assertEqual(res["models"], ["AluY"])
        self.assertEqual(res["matched"], 281)


def aligned_bases(blocks):
    """Query bases aligned to the consensus by blocks of the best block's class and strand."""
    top = max(blocks, key=lambda b: b["score"])
    q = set()
    for b in blocks:
        if b["cls"] == top["cls"] and b["strand"] == top["strand"]:
            q.update(p for p, (c, _) in rc.base_map(b, 0).items() if c is not None)
    return len(q)


class ChainRules(unittest.TestCase):
    def test_tandem_copy_is_inserted(self):
        # Two copies of consensus 1-100, then 101-200: one copy is matched, the other inserted.
        res = chain([block(1, 100, 1), block(101, 200, 1), block(201, 300, 101)])
        self.assertEqual((res["matched"], res["inserted"]), (200, 100))
        self.assertEqual(len(res["ins_blocks"]), 1)
        self.assertEqual(res["ins_blocks"][0][1] - res["ins_blocks"][0][0] + 1, 100)

    def test_minus_strand_mirrors_plus(self):
        plus = chain([block(1, 100, 1), block(101, 200, 1), block(201, 300, 101)])
        # The same element on the minus strand: the query runs 300 -> 1 along the consensus.
        minus = chain([block(201, 300, 1, "-"), block(101, 200, 1, "-"), block(1, 100, 101, "-")])
        self.assertEqual((minus["matched"], minus["inserted"]), (plus["matched"], plus["inserted"]))

    def test_truncated_first_block_does_not_block_the_chain(self):
        # The first block starts mid-way (consensus 101-150); a later full copy (1-150)
        # followed by 151-250 matches more. A first-come rule would keep 50 + 100 bases.
        res = chain([block(1, 50, 101), block(51, 200, 1), block(201, 300, 151)])
        self.assertEqual(res["matched"], 250)
        self.assertEqual(res["inserted"], 50)  # the truncated start, inside the element

    def test_jump_cost(self):
        blocks = [block(1, 100, 1), block(151, 170, 101)]
        self.assertEqual(chain(blocks, jump_cost=30)["matched"], 100)  # 20 bases < cost
        self.assertEqual(chain(blocks, jump_cost=0)["matched"], 120)

    def test_sv_inside_element(self):
        # RepeatMasker joined an element across the SV (query 101-200): all of it inserted.
        res = chain([block(1, 100, 1), block(201, 300, 101)], sv=(100, 200))
        self.assertEqual((res["matched"], res["inserted"]), (0, 100))
        self.assertEqual(res["ins_blocks"], [(101, 200, 100)])

    def test_bases_outside_element_not_counted(self):
        res = chain([block(51, 150, 1)], sv=(0, 200))
        self.assertEqual((res["matched"], res["inserted"]), (100, 0))
        self.assertEqual((res["element_start"], res["element_end"]), (51, 150))


class TableAndDiagram(unittest.TestCase):
    def test_table_round_trip(self):
        res = chain([block(1, 100, 1), block(101, 200, 1), block(201, 300, 101)])
        res.update(sv_id="INS.1", rm_id=7)
        with tempfile.NamedTemporaryFile("w", suffix=".tab", delete=False) as f:
            f.write("#" + "\t".join(rc.TABLE_COLS) + "\n" + rc.format_row(res) + "\n")
        self.addCleanup(os.remove, f.name)
        row = rc.load_chain_table(f.name)["INS.1"][7]
        for col in ("matched", "identical", "inserted", "element_start", "element_end", "ins_blocks"):
            self.assertEqual(row[col], res[col], col)

    def row(self, start, end, ins_blocks=(), matched=10, hit=10, cons_len=100):
        return {"element_start": start, "element_end": end, "ins_blocks": list(ins_blocks),
                "matched": matched, "identical": matched, "inserted": 0, "consensus_hit": hit,
                "consensus_len": cons_len, "repeat": "X", "n_blocks": 1}

    def test_diagram_row(self):
        # SV = query 101-200 on a 10-character diagram (10 bp per character). Same scaling as
        # repeat_annotation.py's create_rm(): (position - sv_start) * scale, truncated.
        self.assertEqual(rc.chain_diagram_row(self.row(101, 199, [(131, 159, 5)]), 100, 0.1, 10),
                         "[..***...]")
        # An element ending on the SV's last base scales to the column past the end: no bracket.
        self.assertEqual(rc.chain_diagram_row(self.row(101, 200), 100, 0.1, 10), "[.........")
        # An element running past both ends has no brackets.
        self.assertEqual(rc.chain_diagram_row(self.row(1, 300, [(101, 200, 5)]), 100, 0.1, 10),
                         "**********")

    def test_block_rows_one_per_element(self):
        # Two selected .out lines of one joined element share an ID: one chain row.
        elements = [{"rm_id": 3, "classification": "LINE"}, {"rm_id": 3, "classification": "LINE"},
                    {"rm_id": 4, "classification": "SINE"}]
        sv_chain = {3: self.row(101, 200), 4: self.row(150, 180)}
        rows = rc.chain_block_rows(sv_chain, elements, 100, 100, 10)
        self.assertEqual(len(rows), 2)
        self.assertTrue(rows[0].startswith(" LINE"))


class WrongBlockCoordinates(unittest.TestCase):
    """rm_offbyone.*: an AluSx block (DEL.9856 ID 290, PB_human on Gadi) that RepeatMasker
    4.2.6 labels 1001-1280 (header, rows and .out) but whose 281 bases sit at 1000-1280."""
    PREFIX = os.path.join(FIXTURES, "rm_offbyone.fa")
    WHOLE_SV = {"DEL.9856": (0, 2268)}

    def chain(self, out_path):
        with unittest.mock.patch("sys.stderr", new_callable=io.StringIO) as err:
            rows = rc.chain_rm_output(out_path, self.PREFIX + ".align", self.WHOLE_SV, rc.DEFAULT_JUMP_COST)
        return rows, err.getvalue()

    def test_placed_from_the_query_sequence(self):
        rows, err = self.chain(self.PREFIX + ".out")
        self.assertEqual(len(rows), 1)
        self.assertEqual((rows[0]["element_start"], rows[0]["element_end"]), (1000, 1280))
        self.assertEqual(rows[0]["matched"], 281)
        self.assertIn("they sit at 1000-1280", err)

    def test_left_out_without_the_query_sequence(self):
        # No N.fa next to the .out: the block cannot be placed, the run goes on without it.
        with tempfile.TemporaryDirectory() as tmp:
            out_path = os.path.join(tmp, "0.fa.out")
            shutil.copy(self.PREFIX + ".out", out_path)
            rows, err = self.chain(out_path)
        self.assertEqual(rows, [])
        self.assertIn("block left out", err)

    def test_masked_bases_match_any_base(self):
        self.assertTrue(rc.rows_match("ACNNT", "GACGTT", 2))
        self.assertFalse(rc.rows_match("ACNNA", "GACGTT", 2))


class Script(unittest.TestCase):
    SCRIPT = os.path.join(SRC_DIR, "rm_align_chain.py")
    INFO = ("chrom\tquerystart\tqueryend\tpos\tend\tsvlen\trelativeID\tcallerID\n"
            "chr1\t1\t3000\t1001\t2000\t-1000\tDEL.57\tcaller.DEL.57\n"
            "chr1\t1\t3000\t1001\t2000\t-1000\tDEL.96\tcaller.DEL.96\n"
            "chr1\t1\t6000\t4001\t4001\t500\tINS.119\tcaller.INS.119\n")

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp)
        self.rm_dir = os.path.join(self.tmp, "rm")
        os.mkdir(self.rm_dir)
        self.info = os.path.join(self.tmp, "info.tab")
        with open(self.info, "w") as f:
            f.write(self.INFO)

    def run_script(self, *extra):
        out = os.path.join(self.tmp, "rm_chain.tab")
        proc = subprocess.run([sys.executable, self.SCRIPT, "--rm_dir", self.rm_dir, "--info", self.info,
                               "--out", out, *extra], capture_output=True, text=True)
        return proc, out

    def test_writes_table(self):
        for ext in ("out", "align"):
            shutil.copy(os.path.join(FIXTURES, f"elements.fa.{ext}"), os.path.join(self.rm_dir, f"0.fa.{ext}"))
        proc, out = self.run_script("--jump_cost", "30")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        table = rc.load_chain_table(out)
        self.assertEqual({(sv, i) for sv in table for i in table[sv]},
                         {("DEL.57", 34), ("DEL.96", 56), ("INS.119", 142)})

    def test_missing_align_is_an_error(self):
        shutil.copy(os.path.join(FIXTURES, "elements.fa.out"), os.path.join(self.rm_dir, "0.fa.out"))
        proc, _ = self.run_script()
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("-a", proc.stderr)


if __name__ == "__main__":
    unittest.main()

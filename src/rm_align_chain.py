"""
Chain RepeatMasker alignments (.align, from RepeatMasker -a) into one alignment per element.

RepeatMasker reports one element as several .align blocks under one .out ID, e.g. one block
per copy of an expanded SVA VNTR, or an L1 split over its _5end/_orf2/_3end models, without
saying how they fit together. For each ID this picks one colinear path through the blocks:
query and consensus both increase in the element's direction, matched bases are maximised and
every jump between blocks costs jump_cost. Query bases inside the element that are not on the
path are inserted, placed at the consensus position reached; e.g. the extra VNTR copies of an
SVA, or an SV that RepeatMasker joined a reference element across.

Consensus frames: blocks are put on the coordinates .out reports for the element. A model
anchored to an .out line of the same ID (sharing a query end with one of its blocks) takes
that line's frame; .out already joins an L1's _5end/_orf2/_3end models into one full-length
frame (e.g. L1M5_orf2 -> L1M5 +2109). Otherwise sibling subfamilies (SVA_A..F tiles under an
SVA_F line) share coordinates, and separately numbered segment models (ERV LTR and -int,
which .out keeps in their own frames) are placed one after another in the element's
direction. Blocks of another class or strand than the element's best block (an Alu's (T)n
tail) are left out of the chain.

Run as a script after RepeatMasker to write the chain table (one row per SV and .out ID whose
blocks touch the SV); repeat_annotation.py reads it to draw the RM chain diagram block.
"""
import argparse
import glob
import os
import re
import sys
import time
from collections import defaultdict
from multiprocessing import Pool

HEADER_RE = re.compile(r'^\s*\d+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+\S+\s+\d+\s+\d+\s+\(\d+\)')
SEQ_RE = re.compile(r'^(C )?\s*(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s*$')
SEGMENT_RE = re.compile(r'_(5end|orf2|3end)$')
INS_BLOCK_MIN = 50  # inserted runs this long (bp) or longer are listed and drawn as '*'
DEFAULT_JUMP_COST = 30

TABLE_COLS = ['sv_id', 'rm_id', 'repeat', 'class_family', 'strand', 'models', 'n_blocks', 'n_left_out',
              'element_start', 'element_end', 'consensus_len', 'matched', 'identical', 'inserted',
              'consensus_hit', 'ins_blocks']


# --- RepeatMasker output parsing ---

def parse_align_header(line):
    f = line.split()
    d = {"score": int(f[0]), "query": f[4], "qbeg": int(f[5]), "qend": int(f[6])}
    if f[8] == "C":
        name = f[9]
        d.update(strand="-", cbeg=int(f[12]), cend=int(f[11]))
        rest = f[13:]
    else:
        name = f[8]
        d.update(strand="+", cbeg=int(f[9]), cend=int(f[10]))
        rest = f[12:]
    d["name"], _, d["cls"] = name.partition("#")
    # The .out ID is the last numeric field (after the internal alignment name).
    d["id"] = int([x for x in rest if x.isdigit()][-1])
    return d


def parse_align(path):
    """
    Blocks of a RepeatMasker .align file: header fields plus q_row and c_row, the gapped
    query and consensus rows. On minus-strand blocks the consensus rows start with "C " and
    the consensus coordinates run downwards.
    """
    blocks = []
    cur = None
    expect_query = True
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if HEADER_RE.match(line) and not SEQ_RE.match(line):
                cur = parse_align_header(line)
                cur["q_row"], cur["c_row"] = [], []
                blocks.append(cur)
                expect_query = True
                continue
            m = SEQ_RE.match(line) if cur is not None else None
            if not m:
                continue
            (cur["q_row"] if expect_query else cur["c_row"]).append(m.group(4))
            expect_query = not expect_query
    for b in blocks:
        b["q_row"] = "".join(b["q_row"])
        b["c_row"] = "".join(b["c_row"])
    return blocks


def parse_out(path):
    """Hits of a RepeatMasker .out file (header lines and 'no repetitive sequences' skipped)."""
    rows = []
    with open(path) as f:
        for line in f:
            fields = line.split()
            if len(fields) < 15 or not fields[0].isdigit():
                continue
            minus = fields[8] == "C"
            cend = int(fields[12])
            left = int((fields[11] if minus else fields[13]).strip("()"))
            rows.append({"query": fields[4], "qbeg": int(fields[5]), "qend": int(fields[6]),
                         "strand": "-" if minus else "+", "name": fields[9], "cls": fields[10],
                         "cbeg": int(fields[13]) if minus else int(fields[11]), "cend": cend,
                         "cons_len": cend + left, "id": int(fields[14])})
    return rows


# --- consensus frames ---

def block_anchor(b, outs):
    """
    (offset, .out name) placing block b on the frame of an .out line of the same ID that
    shares a query end with it (that end's consensus position is exact on both); None if no
    line does. Plus strand: query start <-> consensus start; minus: query start <-> consensus end.
    """
    plus = b["strand"] == "+"
    for r in sorted(outs, key=lambda r: (r["qbeg"] != b["qbeg"]) + (r["qend"] != b["qend"])):
        if r["strand"] != b["strand"]:
            continue
        if r["qbeg"] == b["qbeg"]:
            return (r["cbeg"] - b["cbeg"] if plus else r["cend"] - b["cend"]), r["name"]
        if r["qend"] == b["qend"]:
            return (r["cend"] - b["cend"] if plus else r["cbeg"] - b["cbeg"]), r["name"]
    return None


def segment_group(name):
    """Models numbered separately along one element: L1 _5end/_orf2/_3end, ERV LTR vs -int."""
    m = SEGMENT_RE.search(name)
    if m:
        return m.group(1)
    return "int" if "-int" in name else "main"


def frame_offsets(used, outs, plus):
    """(block, consensus offset) for each block, in the element's direction, on one frame."""
    anchors = {}
    for b in sorted(used, key=lambda b: -b["score"]):
        a = block_anchor(b, outs)
        if a is not None:
            anchors.setdefault(b["name"], a)
    order = sorted(used, key=lambda b: b["qbeg"] if plus else -b["qend"])
    off = {b["name"]: anchors.get(b["name"], (0, None))[0] for b in used}
    # Segment group from the .out name a model is anchored to (.out puts an L1's segment
    # models on one full-length name; HERVH -> HERVH-int), except that an -int model stays
    # -int: .out sometimes names a whole ERV after its LTR (THE1D-int -> MSTA 1-907).
    group = {}
    for b in used:
        g = segment_group(b["name"])
        if g != "int" and b["name"] in anchors:
            g = segment_group(anchors[b["name"]][1])
        group[b["name"]] = g
    if len(set(group.values())) == 1:
        return [(b, off[b["name"]]) for b in order]
    # Several groups left (e.g. 5' LTR, -int, 3' LTR): each run of one group starts right
    # after the consensus reached by the runs before it.
    placed, reached, run_group, shift = [], 0, None, 0
    for b in order:
        if group[b["name"]] != run_group:
            run_group = group[b["name"]]
            shift = reached - (b["cbeg"] + off[b["name"]]) + 1
        o = off[b["name"]] + shift
        placed.append((b, o))
        reached = max(reached, b["cend"] + o)
    return placed


# --- chaining ---

def base_map(b, off):
    """query position -> (consensus position on the element's frame or None, identical?)."""
    out = {}
    q = b["qbeg"]
    step = 1 if b["strand"] == "+" else -1
    c = (b["cbeg"] if b["strand"] == "+" else b["cend"]) + off
    for qa, ca in zip(b["q_row"], b["c_row"]):
        if qa != "-" and ca != "-":
            out[q] = (c, qa.upper() == ca.upper())
            q += 1
            c += step
        elif qa != "-":
            out[q] = (None, False)
            q += 1
        else:
            c += step
    if q - 1 != b["qend"]:
        raise ValueError(f"{b['query']} ID {b['id']}: .align block ends at {q - 1}, header says {b['qend']}")
    return out


def chain_path(maps, plus, jump_cost):
    """
    Best path as (query, consensus, identical) points in walking order. Strictly increasing
    in walking-direction query and in consensus; score = matched bases - jump_cost per jump.
    Continuing to the same block's previous aligned column is free (RM's own indels included).
    """
    pts = []
    for m in maps:
        prev = None
        for q in (sorted(m) if plus else sorted(m, reverse=True)):
            c, same = m[q]
            if c is not None:
                pts.append((q if plus else -q, c, same, q, prev))
                prev = len(pts) - 1
    if not pts:
        return []
    cmin = min(p[1] for p in pts)
    size = max(p[1] for p in pts) - cmin + 2
    tree = [(-1, None)] * (size + 1)  # Fenwick prefix-max over consensus position

    def best_below(i):
        best = (-1, None)
        while i > 0:
            if tree[i][0] > best[0]:
                best = tree[i]
            i -= i & -i
        return best

    def update(i, val):
        while i <= size:
            if val[0] > tree[i][0]:
                tree[i] = val
            i += i & -i

    order = sorted(range(len(pts)), key=lambda i: pts[i][0])
    score, parent = [0] * len(pts), [None] * len(pts)
    i = 0
    while i < len(order):
        # Points at the same query position cannot follow each other: score them all first.
        j = i
        while j < len(order) and pts[order[j]][0] == pts[order[i]][0]:
            j += 1
        for k in order[i:j]:
            c, prev = pts[k][1], pts[k][4]
            best, arg = 1, None
            jump = best_below(c - cmin)  # consensus positions strictly below c
            if jump[1] is not None and jump[0] + 1 - jump_cost > best:
                best, arg = jump[0] + 1 - jump_cost, jump[1]
            if prev is not None and score[prev] + 1 >= best:
                best, arg = score[prev] + 1, prev
            score[k], parent[k] = best, arg
        for k in order[i:j]:
            update(pts[k][1] - cmin + 1, (score[k], k))
        i = j
    k = max(range(len(pts)), key=lambda x: score[x])
    path = []
    while k is not None:
        path.append((pts[k][3], pts[k][1], pts[k][2]))
        k = parent[k]
    return path[::-1]


def chain_element(blocks, outs, sv_start, sv_end, jump_cost):
    """
    Chain one .out ID's blocks; counts are over the SV (query positions sv_start < q <= sv_end,
    the region repeat_annotation.py uses). outs = the ID's .out lines.
    """
    top = max(blocks, key=lambda b: b["score"])
    used = [b for b in blocks if b["cls"] == top["cls"] and b["strand"] == top["strand"]]
    plus = top["strand"] == "+"
    maps = [base_map(b, off) for b, off in frame_offsets(used, outs, plus)]
    path = chain_path(maps, plus, jump_cost)
    # The element is the union of its chained blocks (RM joined them). Inside it, an SV base
    # off the path is inserted at the consensus position reached (bases before the first
    # path point: at the element's start); outside it, an SV base is not part of the element.
    elem_start, elem_end = min(b["qbeg"] for b in used), max(b["qend"] for b in used)
    on_path = {q: (c, same) for q, c, same in path}
    first_c = path[0][1] if path else None
    matched = identical = 0
    cons_hit = set()
    reached, runs, run = first_c, [], None
    for q in (range(elem_start, elem_end + 1) if plus else range(elem_end, elem_start - 1, -1)):
        in_sv = sv_start < q <= sv_end
        if q in on_path:
            reached = on_path[q][0]
            if in_sv:
                matched += 1
                identical += on_path[q][1]
                cons_hit.add(reached)
            if run:
                runs.append(run)
                run = None
        elif in_sv:
            run = [q, q, reached] if run is None else [min(run[0], q), max(run[1], q), run[2]]
    if run:
        runs.append(run)
    outs_main = max(outs, key=lambda r: r["qend"] - r["qbeg"]) if outs else None
    return {
        "repeat": outs_main["name"] if outs_main else top["name"],
        "class_family": outs_main["cls"] if outs_main else top["cls"],
        "strand": top["strand"],
        "models": sorted({b["name"] for b in used}),
        "n_blocks": len(used),
        "n_left_out": len(blocks) - len(used),
        "element_start": elem_start,
        "element_end": elem_end,
        "consensus_len": outs_main["cons_len"] if outs_main else None,
        "matched": matched,
        "identical": identical,
        "inserted": sum(e - s + 1 for s, e, _ in runs),
        "consensus_hit": len(cons_hit),
        "ins_blocks": [(s, e, c) for s, e, c in runs if e - s + 1 >= INS_BLOCK_MIN],
    }


# --- chain table ---

def chain_rm_output(out_path, align_path, sv_regions, jump_cost):
    """Table rows for one RepeatMasker run. sv_regions: sv_id -> (sv_start, sv_end)."""
    outs_by_id = defaultdict(list)
    for r in parse_out(out_path):
        outs_by_id[(r["query"], r["id"])].append(r)
    blocks_by_id = defaultdict(list)
    for b in parse_align(align_path):
        if b["query"] in sv_regions:
            blocks_by_id[(b["query"], b["id"])].append(b)
    rows = []
    for (sv_id, rm_id), blocks in sorted(blocks_by_id.items()):
        sv_start, sv_end = sv_regions[sv_id]
        if not any(b["qbeg"] <= sv_end and b["qend"] > sv_start for b in blocks):
            continue
        res = chain_element(blocks, outs_by_id[(sv_id, rm_id)], sv_start, sv_end, jump_cost)
        res.update(sv_id=sv_id, rm_id=rm_id)
        rows.append(res)
    return rows


def format_row(res):
    vals = dict(res)
    vals["models"] = ",".join(res["models"])
    vals["ins_blocks"] = ";".join(f"{s}-{e}@{c}" for s, e, c in res["ins_blocks"]) or "."
    if vals["consensus_len"] is None:
        vals["consensus_len"] = "."
    return "\t".join(str(vals[c]) for c in TABLE_COLS)


def load_chain_table(path):
    """sv_id -> rm_id (the .out ID, as in column 15 of rm.tab) -> row dict."""
    table = defaultdict(dict)
    with open(path) as f:
        header = f.readline().rstrip("\n").lstrip("#").split("\t")
        for line in f:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            for c in ("rm_id", "n_blocks", "n_left_out", "element_start", "element_end", "matched",
                      "identical", "inserted", "consensus_hit"):
                row[c] = int(row[c])
            row["consensus_len"] = None if row["consensus_len"] == "." else int(row["consensus_len"])
            row["ins_blocks"] = [] if row["ins_blocks"] == "." else [
                tuple(int(x) for x in re.match(r"(\d+)-(\d+)@(-?\d+)", b).groups())
                for b in row["ins_blocks"].split(";")]
            table[row["sv_id"]][row["rm_id"]] = row
    return table


# --- diagram block (called from repeat_annotation.py's draw_diagrams) ---

def chain_metrics(row, sv_len):
    """cov, ident, ins and ref_cov of one chain table row, as percentages / bp."""
    cov = row["matched"] / sv_len * 100
    ident = row["identical"] / row["matched"] * 100 if row["matched"] else 0.0
    ref_cov = min(100.0, row["consensus_hit"] / row["consensus_len"] * 100) if row["consensus_len"] else 0.0
    return cov, ident, row["inserted"], ref_cov


def chain_diagram_row(row, sv_start, sv_scale, diagram_length):
    """
    '[...***...]' over the SV: the element's extent, '*' over inserted blocks >= INS_BLOCK_MIN bp.
    Same scaling and bracket rules as repeat_annotation.py's create_rm(): an end outside the
    SV has no bracket, the dots run to the edge.
    """
    diagram = [' '] * diagram_length
    a = int((row["element_start"] - sv_start) * sv_scale)
    b = int((row["element_end"] - sv_start) * sv_scale)
    for i in range(max(0, a), min(diagram_length - 1, b) + 1):
        diagram[i] = '.'
    if 0 <= a < diagram_length:
        diagram[a] = '['
    if 0 <= b < diagram_length:
        diagram[b] = ']' if b != a else '|'
    for s, e, _ in row["ins_blocks"]:
        lo = max(0, int((s - sv_start) * sv_scale))
        hi = min(diagram_length - 1, int((e - sv_start) * sv_scale))
        for i in range(lo, hi + 1):
            if diagram[i] == '.':
                diagram[i] = '*'
    return ''.join(diagram)


def chain_block_rows(sv_chain, elements, sv_start, sv_len, diagram_length, field_width=20):
    """
    Rows of the RM chain diagram block for one SV: one per .out ID among the selected RM
    elements (several .out lines of one joined element share an ID and get one row), laid
    out like the RM rows of repeat_annotation.py's SV Diagram block. sv_chain: rm_id -> row.
    """
    rows = []
    seen = set()
    sv_scale = diagram_length / sv_len
    for element in elements:
        rm_id = element['rm_id']
        if rm_id in seen or rm_id not in sv_chain:
            continue
        seen.add(rm_id)
        row = sv_chain[rm_id]
        cov, ident, ins, ref_cov = chain_metrics(row, sv_len)
        label = f" {element['classification']}".ljust(field_width)
        diagram = chain_diagram_row(row, sv_start, sv_scale, diagram_length).ljust(diagram_length + 1)
        rows.append(f"{label}{diagram}\t{row['repeat']}\tcov={cov:.1f}%\tident={ident:.1f}%\t"
                    f"ins={ins}\tref_cov={ref_cov:.1f}%\tblocks={row['n_blocks']}\n")
    return rows


# --- script ---

def _chain_one(job):
    out_path, align_path, sv_regions, jump_cost = job
    return chain_rm_output(out_path, align_path, sv_regions, jump_cost)


def read_sv_regions(info_path):
    """sv_id -> (sv_start, sv_end) in query coordinates, as repeat_annotation.py defines the SV."""
    # Imported here, not at the top: repeat_annotation.py imports this module.
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from repeat_annotation import read_sv_info
    return {sv_id: (sv["rel_start"], sv["rel_end"]) for sv_id, sv in read_sv_info(info_path).items()}


def argparser():
    def positive_int(value):
        ivalue = int(value)
        if ivalue <= 0:
            raise argparse.ArgumentTypeError(f"{value} is invalid. It must be greater than zero.")
        return ivalue

    def non_negative_float(value):
        fvalue = float(value)
        if fvalue < 0:
            raise argparse.ArgumentTypeError(f"{value} is invalid. It must be zero or greater.")
        return fvalue

    parser = argparse.ArgumentParser(
        description="Chain RepeatMasker .align blocks into one alignment per element and write the chain table",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument('--rm_dir', required=True, type=str, help="Directory with RepeatMasker *.fa.out and *.fa.align files")
    required_args.add_argument('-i', '--info', required=True, type=str, help="Path to SV info file")
    required_args.add_argument('-o', '--out', required=True, type=str, help="Path to write the chain table")

    optional_args = parser.add_argument_group("optional arguments")
    optional_args.add_argument('--jump_cost', required=False, type=non_negative_float, default=DEFAULT_JUMP_COST,
                               help="Score cost of each jump between alignment blocks when chaining (in matched bases)")
    optional_args.add_argument('-t', '--threads', required=False, type=positive_int, default=1, help="Number of RepeatMasker outputs to chain in parallel")
    optional_args.add_argument('-h', '--help', action='help', help="Show this help message and exit")
    return parser


if __name__ == "__main__":
    start = time.time()
    args = argparser().parse_args()

    out_files = sorted(glob.glob(os.path.join(args.rm_dir, "*.fa.out")))
    if not out_files:
        print(f"Error: no *.fa.out files in {args.rm_dir}", file=sys.stderr)
        sys.exit(1)
    sv_regions = read_sv_regions(args.info)

    jobs = []
    for out_path in out_files:
        align_path = out_path[:-len(".out")] + ".align"
        if not os.path.exists(align_path):
            # RepeatMasker writes no .align when it finds nothing; any hit means a missing -a.
            if parse_out(out_path):
                print(f"Error: {align_path} not found (RepeatMasker must be run with -a)", file=sys.stderr)
                sys.exit(1)
            continue
        jobs.append((out_path, align_path, sv_regions, args.jump_cost))

    with Pool(args.threads) as pool:
        results = pool.map(_chain_one, jobs)

    n_rows = 0
    with open(args.out, "w") as f:
        f.write("#" + "\t".join(TABLE_COLS) + "\n")
        for rows in results:
            for res in rows:
                f.write(format_row(res) + "\n")
                n_rows += 1
    print(f"Info: {n_rows} chained elements from {len(jobs)} RepeatMasker output(s) written to {args.out}")
    print(f"Run time: {time.time() - start:.3f} seconds")

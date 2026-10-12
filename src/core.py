# vim: set expandtab shiftwidth=4 softtabstop=4:
"""Pure-numpy helpers for the "compare" command.

Nothing in this module imports ChimeraX, so it can be unit-tested with a
plain Python interpreter.
"""
import numpy

SS_LETTERS = {0: "-", 1: "H", 2: "E"}


def rmsd(dists):
    """Root-mean-square of an array of pair distances."""
    d = numpy.asarray(dists, dtype=float)
    if d.size == 0:
        return float("nan")
    return float(numpy.sqrt((d * d).mean()))


def tm_score(dists, norm_len):
    """TM-score of an existing superposition.

    dists are the distances between aligned CA pairs and norm_len is the
    length of the protein the score is normalised by.  The superposition is
    taken as given (it is not re-optimised the way TM-align does), so the
    value is a lower bound on the true TM-score.
    """
    d = numpy.asarray(dists, dtype=float)
    if norm_len <= 0 or d.size == 0:
        return float("nan")
    d0 = 1.24 * max(norm_len - 15, 1) ** (1.0 / 3.0) - 1.8
    d0 = max(d0, 0.5)
    return float((1.0 / (1.0 + (d / d0) ** 2)).sum() / norm_len)


def percent_identity(seq_a, seq_b):
    """Percent of aligned pairs with the same one-letter code."""
    if not seq_a:
        return float("nan")
    same = sum(1 for a, b in zip(seq_a, seq_b) if a == b)
    return 100.0 * same / len(seq_a)


def ss_agreement(ss_a, ss_b):
    """Percent of aligned pairs with the same secondary-structure class."""
    if not ss_a:
        return float("nan")
    same = sum(1 for a, b in zip(ss_a, ss_b) if a == b)
    return 100.0 * same / len(ss_a)


def ss_composition(ss):
    """(% helix, % strand, % coil) of a sequence of H/E/- letters."""
    n = len(ss)
    if n == 0:
        return (float("nan"),) * 3
    helix = 100.0 * sum(1 for s in ss if s == "H") / n
    strand = 100.0 * sum(1 for s in ss if s == "E") / n
    return helix, strand, 100.0 - helix - strand


def gapped_columns(n_ref, n_match, pairs):
    """Rebuild a gapped alignment from paired indices.

    pairs is a list of (ref_index, match_index), increasing in both.
    Returns a list of (ref_index or None, match_index or None) columns that
    also covers the unpaired residues of both chains.
    """
    cols = []
    i = j = 0
    for pi, pj in pairs:
        while i < pi:
            cols.append((i, None))
            i += 1
        while j < pj:
            cols.append((None, j))
            j += 1
        cols.append((pi, pj))
        i, j = pi + 1, pj + 1
    while i < n_ref:
        cols.append((i, None))
        i += 1
    while j < n_match:
        cols.append((None, j))
        j += 1
    return cols


def format_alignment(ref_name, match_name, ref_seq, ref_ss, match_seq, match_ss,
        pairs, width=60):
    """Text alignment with a secondary-structure line above/below each
    sequence and a midline marking identities (|) and SS disagreement (^)."""
    cols = gapped_columns(len(ref_seq), len(match_seq), pairs)
    rows = {"rss": [], "rseq": [], "mid": [], "mseq": [], "mss": [], "flag": []}
    for i, j in cols:
        rows["rseq"].append(ref_seq[i] if i is not None else "-")
        rows["rss"].append(ref_ss[i] if i is not None else " ")
        rows["mseq"].append(match_seq[j] if j is not None else "-")
        rows["mss"].append(match_ss[j] if j is not None else " ")
        if i is not None and j is not None:
            rows["mid"].append("|" if ref_seq[i] == match_seq[j] else " ")
            rows["flag"].append(" " if ref_ss[i] == match_ss[j] else "^")
        else:
            rows["mid"].append(" ")
            rows["flag"].append(" ")
    label_w = max(len(ref_name), len(match_name), 2) + 4
    lines = []
    for start in range(0, len(cols), width):
        end = start + width
        for key, label in (("rss", ref_name + " SS"), ("rseq", ref_name),
                ("mid", ""), ("mseq", match_name), ("mss", match_name + " SS"),
                ("flag", "SS differs")):
            lines.append(label.ljust(max(label_w, 11)) + " " + "".join(rows[key][start:end]))
        lines.append("")
    return "\n".join(lines)


def substitutions(ref_seq, ref_numbers, match_seq, match_numbers, pairs, dists):
    """List the aligned positions whose residue type differs.

    Returns tuples (ref_code, ref_number, match_code, match_number, distance).
    """
    subs = []
    for (i, j), d in zip(pairs, dists):
        if ref_seq[i] != match_seq[j]:
            subs.append((ref_seq[i], ref_numbers[i], match_seq[j], match_numbers[j], float(d)))
    return subs


def star_columns(ref_len, others):
    """Merge pairwise alignments to one reference into a multiple alignment.

    others is a list of (match_len, pairs), one per non-reference structure,
    where pairs is a list of (ref_index, match_index) increasing in both.
    Returns a list of columns (ref_index or None, [match_index or None, ...]).
    Residues a structure has between two reference positions go into insert
    columns; inserts from different structures share those columns but are
    not aligned to each other, only to the reference.
    """
    inserts, paired = [], []
    for match_len, pairs in others:
        to_ref = {j: i for i, j in pairs}
        ins = [[] for _ in range(ref_len + 1)]
        pending = []
        for j in range(match_len):
            if j in to_ref:
                ins[to_ref[j]] = pending
                pending = []
            else:
                pending.append(j)
        ins[ref_len] = pending
        inserts.append(ins)
        paired.append({i: j for i, j in pairs})
    columns = []
    for slot in range(ref_len + 1):
        n_ins = max([len(ins[slot]) for ins in inserts] or [0])
        for k in range(n_ins):
            columns.append((None, [ins[slot][k] if k < len(ins[slot]) else None
                for ins in inserts]))
        if slot < ref_len:
            columns.append((slot, [p.get(slot) for p in paired]))
    return columns


def multiple_rows(seqs, sss, columns):
    """Gapped sequence and SS strings for every structure, plus two marker
    lines: '*' where every structure has the same residue, '=' where every
    structure has the same secondary-structure class.  seqs[0] / sss[0] are
    the reference."""
    n = len(seqs)
    seq_rows = [[] for _ in range(n)]
    ss_rows = [[] for _ in range(n)]
    same_res, same_ss = [], []
    for ref_i, match_js in columns:
        idx = [ref_i] + list(match_js)
        for k, i in enumerate(idx):
            seq_rows[k].append(seqs[k][i] if i is not None else "-")
            ss_rows[k].append(sss[k][i] if i is not None else " ")
        full = all(i is not None for i in idx)
        same_res.append("*" if full and len({seqs[k][i] for k, i in enumerate(idx)}) == 1
            else " ")
        same_ss.append("=" if full and len({sss[k][i] for k, i in enumerate(idx)}) == 1
            else " ")
    return (["".join(r) for r in seq_rows], ["".join(r) for r in ss_rows],
        "".join(same_res), "".join(same_ss))


def multiple_stats(columns, same_res, same_ss):
    """(columns where all structures are aligned, of those: identical residue,
    same SS class)."""
    full = [c for c, (ref_i, js) in enumerate(columns)
        if ref_i is not None and all(j is not None for j in js)]
    return (len(full), sum(1 for c in full if same_res[c] == "*"),
        sum(1 for c in full if same_ss[c] == "="))


def format_multiple(names, seq_rows, ss_rows, same_res, same_ss, width=60):
    """Blocked text of a multiple alignment with an SS line under each
    sequence."""
    label_w = max(max(len(n) for n in names) + 4, 14)
    lines = []
    for start in range(0, len(same_res), width):
        end = start + width
        for name, seq, ss in zip(names, seq_rows, ss_rows):
            lines.append(name.ljust(label_w) + seq[start:end])
            lines.append((name + " SS").ljust(label_w) + ss[start:end])
        lines.append("identical".ljust(label_w) + same_res[start:end])
        lines.append("same SS".ljust(label_w) + same_ss[start:end])
        lines.append("")
    return "\n".join(lines)

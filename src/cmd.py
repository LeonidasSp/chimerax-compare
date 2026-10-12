# vim: set expandtab shiftwidth=4 softtabstop=4:
"""The "compare" command: compare the secondary (2D) and tertiary (3D)
structure of two or more proteins in one step."""
import os

import numpy

from chimerax.core.errors import UserError

from . import core

DEV_ATTR = "compare_dev"


def _label(m):
    return "%s (#%s)" % (m.name, m.id_string)


def _file_stem(m):
    keep = [c if c.isalnum() or c in "-_" else "_" for c in m.name]
    return "%s_%s" % (m.id_string.replace(".", "_"), "".join(keep))


def _one_letter(r):
    code = getattr(r, "one_letter_code", None)
    if code:
        return code
    from chimerax.atomic import Sequence
    return Sequence.rname3to1(r.name)


def _chain_info(residue):
    """Sequence, SS string and numbering of the chain a residue belongs to."""
    chain = residue.chain
    residues = [r for r in chain.existing_residues]
    seq = [_one_letter(r) for r in residues]
    ss = [core.SS_LETTERS.get(r.ss_type, "-") for r in residues]
    numbers = [r.number for r in residues]
    index = {r: i for i, r in enumerate(residues)}
    return residues, seq, ss, numbers, index


def _matchmaker(session, match, ref):
    """Superpose match onto ref with matchmaker and return its result dict."""
    from chimerax.core.commands import run
    results = run(session, "matchmaker %s to %s" % (match.atomspec, ref.atomspec), log=False)
    while isinstance(results, (list, tuple)) and results and isinstance(results[0], (list, tuple)):
        results = results[0]
    if not results:
        raise UserError("matchmaker could not pair %s with %s" % (_label(match), _label(ref)))
    return results[0]


def _pair_stats(session, match, ref):
    """Superpose match onto ref and measure the pair.

    The match model is left superposed on ref.
    """
    res = _matchmaker(session, match, ref)
    ref_atoms, match_atoms = res["full ref atoms"], res["full match atoms"]
    if len(ref_atoms) == 0:
        raise UserError("No residues could be aligned between %s and %s"
            % (_label(match), _label(ref)))
    dists = numpy.linalg.norm(ref_atoms.scene_coords - match_atoms.scene_coords, axis=1)
    ref_res, ref_seq, ref_ss, ref_nums, ref_index = _chain_info(ref_atoms[0].residue)
    m_res, m_seq, m_ss, m_nums, m_index = _chain_info(match_atoms[0].residue)
    pairs = [(ref_index[ra.residue], m_index[ma.residue])
        for ra, ma in zip(ref_atoms, match_atoms)]
    pair_ref_seq = [ref_seq[i] for i, j in pairs]
    pair_m_seq = [m_seq[j] for i, j in pairs]
    pair_ref_ss = [ref_ss[i] for i, j in pairs]
    pair_m_ss = [m_ss[j] for i, j in pairs]
    return {
        "ref": ref, "match": match,
        "pruned_rmsd": res["final RMSD"], "pruned_pairs": len(res["final ref atoms"]),
        "full_rmsd": core.rmsd(dists), "full_pairs": len(pairs),
        "tm_score": core.tm_score(dists, len(ref_seq)),
        "identity": core.percent_identity(pair_ref_seq, pair_m_seq),
        "ss_agreement": core.ss_agreement(pair_ref_ss, pair_m_ss),
        "dists": dists, "pairs": pairs,
        "ref_residues": ref_res, "ref_seq": ref_seq, "ref_ss": ref_ss, "ref_nums": ref_nums,
        "match_residues": m_res, "match_seq": m_seq, "match_ss": m_ss, "match_nums": m_nums,
    }


def _html_table(header, rows):
    cells = ["<tr>" + "".join("<th>%s</th>" % h for h in header) + "</tr>"]
    for row in rows:
        cells.append("<tr>" + "".join("<td>%s</td>" % c for c in row) + "</tr>")
    return '<table border="1" cellpadding="4" cellspacing="0">%s</table>' % "".join(cells)


def _write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(str(c) for c in row) + "\n")


def _plot(path, ref, to_ref, max_deviation):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot
    fig, ax = pyplot.subplots(figsize=(10, 3.5))
    # Disordered tails can be tens of angstroms away; clip the axis so the
    # ordered regions stay readable.
    top = 4 * max_deviation
    clipped = False
    for st in to_ref:
        x = [st["ref_nums"][i] for i, j in st["pairs"]]
        clipped = clipped or bool((st["dists"] > top).any())
        ax.plot(x, numpy.minimum(st["dists"], top), linewidth=1, label=st["match"].name)
    ax.set_ylim(0, top * 1.02)
    ax.set_xlabel("Residue number in %s" % ref.name)
    ax.set_ylabel(u"Cα deviation (Å)" + (", clipped at %g" % top if clipped else ""))
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    pyplot.close(fig)


def compare(session, structures, to=None, dssp=True, color=True, max_deviation=5.0,
        substitutions=None, report=None):
    """Compare the secondary and tertiary structure of several proteins.

    Every structure is superposed on every other one with matchmaker to fill
    pairwise tables, then all are left superposed on the reference.
    """
    from chimerax.core.commands import run
    from chimerax.atomic import Residue
    models = list(structures)
    if to is None:
        to = models[0]
    if to not in models:
        models.insert(0, to)
    others = [m for m in models if m is not to]
    if not others:
        raise UserError("Need at least two structures to compare")

    if dssp:
        for m in models:
            run(session, "dssp %s" % m.atomspec, log=False)

    # All-against-all first, reference pairs last so that every model ends
    # up superposed on the reference.
    stats = {}
    for i, a in enumerate(models):
        for b in models[i + 1:]:
            if a is to or b is to:
                continue
            stats[(a, b)] = _pair_stats(session, b, a)
    to_ref = [_pair_stats(session, m, to) for m in others]
    for st in to_ref:
        stats[(to, st["match"])] = st

    def pair(a, b):
        return stats.get((a, b)) or stats.get((b, a))

    # --- per-residue deviation as a residue attribute, for colouring -------
    Residue.register_attr(session, DEV_ATTR, "compare", attr_type=float)
    ref_dev = {}
    for st in to_ref:
        for (i, j), d in zip(st["pairs"], st["dists"]):
            setattr(st["match_residues"][j], DEV_ATTR, float(d))
            ref_dev.setdefault(st["ref_residues"][i], []).append(float(d))
    for r, values in ref_dev.items():
        setattr(r, DEV_ATTR, sum(values) / len(values))
    if color:
        for m in models:
            run(session, "color byattribute r:%s %s target cs palette blue:white:red "
                "range 0,%g noValueColor gray" % (DEV_ATTR, m.atomspec, max_deviation), log=False)

    # --- tables ------------------------------------------------------------
    summary_header = ["Structure", "Residues", "% helix", "% strand", "% coil",
        "Aligned to reference", "% identity", "Pruned RMSD (pairs)", "All-pairs RMSD",
        "TM-score", "% SS agreement"]
    h, e, c = core.ss_composition(to_ref[0]["ref_ss"])
    summary = [[_label(to) + " [reference]", len(to_ref[0]["ref_seq"]),
        "%.1f" % h, "%.1f" % e, "%.1f" % c, "-", "-", "-", "-", "-", "-"]]
    for st in to_ref:
        h, e, c = core.ss_composition(st["match_ss"])
        summary.append([_label(st["match"]), len(st["match_seq"]),
            "%.1f" % h, "%.1f" % e, "%.1f" % c, st["full_pairs"], "%.1f" % st["identity"],
            "%.2f (%d)" % (st["pruned_rmsd"], st["pruned_pairs"]), "%.2f" % st["full_rmsd"],
            "%.3f" % st["tm_score"], "%.1f" % st["ss_agreement"]])

    def matrix(key, fmt):
        rows = []
        for a in models:
            row = [_label(a)]
            for b in models:
                row.append("-" if a is b else fmt % pair(a, b)[key])
            rows.append(row)
        return rows

    matrix_header = [""] + [_label(m) for m in models]
    matrices = [
        ("Pruned RMSD (Å)", "pruned_rmsd", "%.2f"),
        ("All-pairs RMSD (Å)", "full_rmsd", "%.2f"),
        ("Sequence identity over aligned pairs (%)", "identity", "%.1f"),
        ("Secondary-structure agreement (%)", "ss_agreement", "%.1f"),
    ]

    # --- substitutions (mutant versus wild type) ---------------------------
    sub_tables = []
    for st in to_ref:
        subs = core.substitutions(st["ref_seq"], st["ref_nums"], st["match_seq"],
            st["match_nums"], st["pairs"], st["dists"])
        show = substitutions if substitutions is not None else (0 < len(subs) <= 25)
        sub_tables.append((st, subs, show))

    # --- multiple alignment, anchored on the reference ---------------------
    names = [to.name] + [st["match"].name for st in to_ref]
    seqs = [to_ref[0]["ref_seq"]] + [st["match_seq"] for st in to_ref]
    sss = [to_ref[0]["ref_ss"]] + [st["match_ss"] for st in to_ref]
    columns = core.star_columns(len(seqs[0]),
        [(len(st["match_seq"]), st["pairs"]) for st in to_ref])
    seq_rows, ss_rows, same_res, same_ss = core.multiple_rows(seqs, sss, columns)
    n_full, n_ident, n_ss = core.multiple_stats(columns, same_res, same_ss)

    log = session.logger
    log.info("<h3>Structure comparison, reference %s</h3>" % _label(to), is_html=True)
    log.info(_html_table(summary_header, summary), is_html=True)
    for title, key, fmt in matrices:
        log.info("<p><b>%s</b></p>" % title + _html_table(matrix_header, matrix(key, fmt)),
            is_html=True)
    for st, subs, show in sub_tables:
        if not show:
            continue
        rows = [["%s%d" % (a, an), "%s%d" % (b, bn), "%.2f" % d] for a, an, b, bn, d in subs]
        log.info("<p><b>Substitutions in %s relative to %s</b></p>"
            % (_label(st["match"]), _label(to))
            + _html_table([to.name, st["match"].name, u"Cα deviation (Å)"], rows),
            is_html=True)
    if n_full:
        log.info("Multiple alignment of %d structures: %d positions aligned in all; "
            "%d (%.1f%%) have the same residue and %d (%.1f%%) the same secondary structure "
            "in all." % (len(names), n_full, n_ident, 100.0 * n_ident / n_full,
            n_ss, 100.0 * n_ss / n_full))
    log.info("TM-score is computed on the matchmaker superposition and normalised by the "
        "reference length; it is a lower bound on the TM-align value. Residue attribute "
        "'%s' holds the per-residue Cα deviation." % DEV_ATTR)

    # --- files -------------------------------------------------------------
    if report:
        os.makedirs(report, exist_ok=True)
        _write_tsv(os.path.join(report, "summary.tsv"), summary_header, summary)
        for title, key, fmt in matrices:
            _write_tsv(os.path.join(report, "matrix_%s.tsv" % key), matrix_header,
                matrix(key, fmt))
        for st, subs, show in sub_tables:
            stem = _file_stem(st["match"])
            rows = []
            for (i, j), d in zip(st["pairs"], st["dists"]):
                rows.append([st["ref_nums"][i], st["ref_seq"][i], st["ref_ss"][i],
                    st["match_nums"][j], st["match_seq"][j], st["match_ss"][j], "%.3f" % d])
            _write_tsv(os.path.join(report, "per_residue_%s.tsv" % stem),
                ["ref_number", "ref_residue", "ref_ss", "match_number", "match_residue",
                 "match_ss", "ca_deviation"], rows)
            _write_tsv(os.path.join(report, "substitutions_%s.tsv" % stem),
                ["ref_residue", "ref_number", "match_residue", "match_number", "ca_deviation"],
                [[a, an, b, bn, "%.3f" % d] for a, an, b, bn, d in subs])
            with open(os.path.join(report, "alignment_%s.txt" % stem), "w",
                    encoding="utf-8") as f:
                f.write(core.format_alignment(to.name, st["match"].name, st["ref_seq"],
                    st["ref_ss"], st["match_seq"], st["match_ss"], st["pairs"]) + "\n")
        with open(os.path.join(report, "multiple_alignment.txt"), "w", encoding="utf-8") as f:
            f.write("Structure-based multiple alignment anchored on %s.\n"
                "Rows are aligned to the reference; residues in insert columns are not "
                "aligned to each other.\n"
                "Positions aligned in all %d structures: %d; identical residue: %d; "
                "same secondary structure: %d\n\n"
                % (to.name, len(names), n_full, n_ident, n_ss))
            f.write(core.format_multiple(names, seq_rows, ss_rows, same_res, same_ss) + "\n")
        with open(os.path.join(report, "multiple_alignment.fasta"), "w",
                encoding="utf-8") as f:
            for name, row in zip(names, seq_rows):
                f.write(">%s\n%s\n" % (name, row))
        dev = [dict(zip((i for i, j in st["pairs"]), st["dists"])) for st in to_ref]
        rows = []
        for c, (ref_i, js) in enumerate(columns):
            if ref_i is None:
                continue
            ds = [d[ref_i] for d in dev if ref_i in d]
            rows.append([to_ref[0]["ref_nums"][ref_i]]
                + [r[c] for r in seq_rows] + [r[c].strip() for r in ss_rows]
                + ["yes" if same_res[c] == "*" else "", "yes" if same_ss[c] == "=" else "",
                   "%.3f" % max(ds) if ds else ""])
        _write_tsv(os.path.join(report, "conservation.tsv"),
            ["ref_number"] + names + [n + "_ss" for n in names]
            + ["identical_in_all", "same_ss_in_all", "max_ca_deviation"], rows)
        _plot(os.path.join(report, "deviation.png"), to, to_ref, max_deviation)
        log.info("Comparison report written to %s" % report)
    return stats


def register_command(logger):
    from chimerax.core.commands import CmdDesc, register, BoolArg, FloatArg, SaveFolderNameArg
    from chimerax.atomic import AtomicStructuresArg, AtomicStructureArg
    desc = CmdDesc(
        required=[("structures", AtomicStructuresArg)],
        keyword=[("to", AtomicStructureArg), ("dssp", BoolArg), ("color", BoolArg),
            ("max_deviation", FloatArg), ("substitutions", BoolArg),
            ("report", SaveFolderNameArg)],
        synopsis="Compare secondary and tertiary structure of two or more proteins")
    register("compare", desc, compare, logger=logger)

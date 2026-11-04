# Changelog

## 0.1.0 (2026-11-04)

First release.

- `compare` command: superposes two or more protein structures with
  matchmaker and reports sequence identity, pruned and all-pairs Cα RMSD,
  TM-score, secondary-structure content and agreement.
- All-against-all matrices of RMSD, identity and secondary-structure agreement.
- Per-residue Cα deviation coloured on the structures and stored as the
  residue attribute `compare_dev`.
- Structure-based multiple alignment with secondary structure, identity and
  same-SS markers (text and aligned FASTA).
- Substitution list against the reference, for mutant versus wild type.
- `report <folder>` writes everything as TSV / text / FASTA plus a deviation
  plot.

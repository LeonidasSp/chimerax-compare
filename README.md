# ChimeraX-CompareStructures

A UCSF ChimeraX bundle that adds one command, `compare`, to compare the
secondary (2D) and tertiary (3D) structure of two or more proteins in a single
step: orthologues from several species, or a mutant against its wild type.

ChimeraX can already superpose structures (`matchmaker`), but it reports one
pair at a time and gives no combined view of sequence, secondary structure and
3D deviation. `compare` runs every pair and puts the results in one place.

## What it reports

- Summary per structure: length, % helix / strand / coil, residues aligned to
  the reference, % identity, pruned and all-pairs CA RMSD, TM-score,
  % secondary-structure agreement.
- All-against-all matrices of RMSD, identity and secondary-structure agreement.
- Per-residue CA deviation from the reference, stored as the residue attribute
  `compare_dev` and coloured on the cartoons (blue = same, red = different).
- The substitutions between each structure and the reference with the local
  deviation at each one (shown automatically for close variants such as point
  mutants).
- A structure-based multiple alignment of all structures, anchored on the
  reference, with the secondary structure of each underneath and markers for
  positions that have the same residue or the same secondary structure in all.
  Rows are aligned to the reference; residues in insert columns are not
  aligned to each other.
- With `report <folder>`: TSV tables, pairwise and multiple alignments (text
  and aligned FASTA), per-position conservation, and a per-residue deviation
  plot.

## Install

In ChimeraX (1.4 or newer; developed on 1.12):

    devel install "path/to/chimerax-compare"

## Use

    compare #1,2,3                      # reference is the first model
    compare #2,3 to #1 report results   # choose the reference, write files
    compare #2 to #1 substitutions true # mutant against wild type

Options: `to`, `dssp` (recompute secondary structure, default true), `color`
(default true), `maxDeviation` (red end of the colour scale in angstroms,
default 5), `substitutions`, `report`.

## Mutant against wild type

`swapaa` only replaces the side chain, so a `swapaa` mutant has the same
backbone as the wild type and every deviation is 0. To see a structural effect,
predict the mutant sequence separately (AlphaFold / ColabFold / ESMFold), open
both models and run `compare`.

## Notes

- TM-score is computed on the matchmaker superposition and normalised by the
  reference length. It is not re-optimised as TM-align does, so it is a lower
  bound.
- All-pairs RMSD is dominated by disordered, low-confidence regions of
  predicted models. Use the pruned RMSD and the per-residue plot for the folded
  core.

## Tests

    python tests/test_core.py

## Example

Compare three AlphaFold models of orthologous proteins and write a report:

    alphafold fetch P69905
    alphafold fetch P01942
    alphafold fetch P01958
    compare #1,2,3 report hemoglobin_alpha

## Licence

MIT. See `LICENSE`.

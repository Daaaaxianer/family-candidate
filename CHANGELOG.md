# Changelog

## Unreleased

### Features

- HMMER-based candidate search with optional species-specific HMM refinement.
- Optional BLASTP search or import of existing similarity search results.
- InterProScan TSV import with required and forbidden signature rules.
- Candidate evidence tables and accepted/review/rejected exports.
- Optional protein–gene–CDS mapping and gene representative sequences.
- Input hashes, external tool information and stage logs for reproducibility.

### Fixes

- Apply the configured final E-value threshold.
- Match FASTA identifiers exactly and reject duplicate identifiers.
- Select qualifying seed domains by independent E-value, score and coverage.
- Validate arguments and stop when an external command fails.
- Reject nonempty output directories and existing legacy outputs.

### Compatibility

- Retain six positional arguments and legacy output filenames for Shell entrypoints.
- Route legacy Perl filenames through the shared Python implementation.
- Require MUSCLE v5 for refinement; the separate MUSCLE v3 workflow is retired.
- Distinguish unvalidated candidates from members passing configured domain rules.

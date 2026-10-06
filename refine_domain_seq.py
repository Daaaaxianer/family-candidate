#!/usr/bin/env python3
"""Compatibility command: select the best qualifying domain per protein/model."""
import argparse
import sys
from family_core import read_fasta, read_domtbl, best_domains, domain_records, write_fasta, positive_float


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('hmmout')
    p.add_argument('proteins')
    p.add_argument('outfile')
    p.add_argument('evalue', type=positive_float)
    p.add_argument('--domain-evalue', type=positive_float)
    p.add_argument('--min-coverage', type=float, default=0.0)
    args = p.parse_args()
    if not 0 <= args.min_coverage <= 1:
        p.error('--min-coverage must be between 0 and 1')
    try:
        hits = best_domains(read_domtbl(args.hmmout), args.evalue,
                            args.domain_evalue or args.evalue, args.min_coverage)
        records = domain_records(hits, read_fasta(args.proteins))
        write_fasta(args.outfile, records)
        print(f'Extracted {len(records)} qualifying protein/model domains')
    except (ValueError, OSError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())

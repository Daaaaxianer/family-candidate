#!/usr/bin/env python3
"""Extract FASTA records by exact ID; fail before writing if any ID is missing."""
import argparse
import sys
from family_core import read_ids, read_fasta, extract_exact, write_fasta


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('idfile')
    p.add_argument('infastafile')
    p.add_argument('outfastafile')
    args = p.parse_args()
    try:
        records = extract_exact(read_ids(args.idfile), read_fasta(args.infastafile))
        write_fasta(args.outfastafile, records)
        print(f'Extracted {len(records)} sequences by exact ID')
    except (ValueError, OSError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())

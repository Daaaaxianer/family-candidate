#!/usr/bin/env bash
# Legacy filename retained; uses the corrected Python implementation.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$script_dir/batch_pfam_hmmer_to_family_candidate_search.sh" "$@"

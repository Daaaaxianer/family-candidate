#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 6 ]; then
    echo "Usage: bash $0 <hmm_seed> <protein.fa> <cds.fa> <seed_evalue> <final_evalue> <outname>" >&2
    exit 2
fi
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python_exe="${PYTHON:-python3}"
name="$6"
for target in "${name}.id.txt" "${name}.protein.fasta" "${name}.cds.fasta"; do
    if [ -e "$target" ]; then
        echo "Refusing to overwrite existing output: $target" >&2
        exit 2
    fi
done
"$python_exe" "$script_dir/family_candidate.py" \
    --hmm "$1" --proteins "$2" --cds "$3" \
    --seed-evalue "$4" --seed-domain-evalue "$4" \
    --sequence-evalue "$5" --domain-evalue "$5" \
    --refine --outdir "${name}.run"
cp -- "${name}.run/candidates.ids.txt" "${name}.id.txt"
cp -- "${name}.run/candidates.protein.fasta" "${name}.protein.fasta"
cp -- "${name}.run/candidates.cds.fasta" "${name}.cds.fasta"
echo "Candidates: ${name}.id.txt, ${name}.protein.fasta, ${name}.cds.fasta"
echo "Evidence and logs: ${name}.run (unvalidated candidates require review)"

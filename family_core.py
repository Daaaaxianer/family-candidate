"""Validated parsers shared by the pipeline and compatibility commands (stdlib only)."""
from dataclasses import dataclass
from pathlib import Path
import hashlib
import math
import re


def positive_float(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise ValueError('threshold must be finite and greater than zero')
    return number


def read_fasta(path):
    records = {}
    header, chunks = None, []
    def save():
        if header is None:
            return
        identifier = header.split()[0]
        if identifier in records:
            raise ValueError(f'{path}: duplicate FASTA ID {identifier}')
        sequence = ''.join(chunks)
        if not sequence:
            raise ValueError(f'{path}: empty sequence {identifier}')
        records[identifier] = (header, sequence)
    with open(path, encoding='utf-8-sig') as handle:
        for number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                save()
                header, chunks = line[1:].strip(), []
                if not header:
                    raise ValueError(f'{path}:{number}: empty FASTA header')
            else:
                if header is None:
                    raise ValueError(f'{path}:{number}: sequence before FASTA header')
                chunks.append(''.join(line.split()))
    save()
    if not records:
        raise ValueError(f'{path}: no FASTA records')
    return records


def write_fasta(path, records):
    with open(path, 'w', encoding='utf-8', newline='\n') as handle:
        for header, sequence in records:
            handle.write('>' + header + '\n')
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index:index+80] + '\n')


def read_ids(path):
    # Preserve the requested order and remove duplicate/blank lines.
    return list(dict.fromkeys(line.strip() for line in
                             Path(path).read_text(encoding='utf-8-sig').splitlines()
                             if line.strip()))


def extract_exact(ids, records):
    missing = [identifier for identifier in ids if identifier not in records]
    if missing:
        raise ValueError('IDs not found: ' + ', '.join(missing[:20]))
    return [records[identifier] for identifier in ids]


@dataclass(frozen=True)
class DomainHit:
    target: str
    target_length: int
    model: str
    accession: str
    model_length: int
    sequence_evalue: float
    independent_evalue: float
    score: float
    hmm_start: int
    hmm_end: int
    start: int
    end: int
    source: str = 'pfam_hmm'

    @property
    def coverage(self):
        return (self.hmm_end - self.hmm_start + 1) / self.model_length

    @property
    def model_id(self):
        return normalize_domain(self.accession if self.accession != '-' else self.model)


def normalize_domain(identifier):
    # Remove Pfam accession version, without altering other domain identifiers.
    return re.sub(r'^(PF\d+)\.\d+$', r'\1', identifier.strip())


def read_domtbl(path, source='pfam_hmm'):
    hits = []
    with open(path, encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split(maxsplit=22)
            try:
                if len(fields) < 22:
                    raise ValueError('expected at least 22 columns')
                hit = DomainHit(fields[0], int(fields[2]), fields[3], fields[4],
                                int(fields[5]), float(fields[6]), float(fields[12]),
                                float(fields[13]), int(fields[15]), int(fields[16]),
                                int(fields[17]), int(fields[18]), source)
                if not (1 <= hit.start <= hit.end <= hit.target_length and
                        1 <= hit.hmm_start <= hit.hmm_end <= hit.model_length):
                    raise ValueError('domain coordinates outside sequence/model')
                if any(not math.isfinite(x) or x < 0 for x in
                       (hit.sequence_evalue, hit.independent_evalue)):
                    raise ValueError('invalid E-value')
                if not math.isfinite(hit.score):
                    raise ValueError('invalid domain score')
            except (ValueError, IndexError) as error:
                raise ValueError(f'{path}:{number}: invalid domtblout: {error}') from error
            hits.append(hit)
    return hits


def passes(hit, sequence_evalue, domain_evalue, coverage):
    return (hit.sequence_evalue <= sequence_evalue and
            hit.independent_evalue <= domain_evalue and hit.coverage >= coverage)


def best_domains(hits, sequence_evalue, domain_evalue, coverage):
    # One best qualifying domain per protein AND model; input order is irrelevant.
    best = {}
    for hit in hits:
        if not passes(hit, sequence_evalue, domain_evalue, coverage):
            continue
        key = (hit.target, hit.model_id)
        quality = (hit.independent_evalue, -hit.score, -hit.coverage, hit.start, hit.end)
        old = best.get(key)
        if old is None or quality < (old.independent_evalue, -old.score, -old.coverage,
                                     old.start, old.end):
            best[key] = hit
    return [best[key] for key in sorted(best)]


def domain_records(hits, proteins):
    records = []
    for hit in hits:
        if hit.target not in proteins:
            raise ValueError(f'HMM target not found in protein FASTA: {hit.target}')
        sequence = proteins[hit.target][1]
        if len(sequence) != hit.target_length:
            raise ValueError(f'HMM target length mismatch: {hit.target}')
        records.append((f'{hit.target}|{hit.model_id}|{hit.start}-{hit.end}',
                        sequence[hit.start-1:hit.end]))
    return records


def read_interpro(path, proteins):
    """Read standard InterProScan TSV; collect signature + InterPro accessions."""
    domains = {}
    with open(path, encoding='utf-8-sig') as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\r\n').split('\t')
            try:
                if len(fields) < 11:
                    raise ValueError('expected at least 11 InterProScan TSV columns')
                identifier = fields[0]
                if identifier not in proteins:
                    raise ValueError(f'unknown protein ID: {identifier}')
                if int(fields[2]) != len(proteins[identifier][1]):
                    raise ValueError(f'protein length mismatch: {identifier}')
                expected_md5 = hashlib.md5(proteins[identifier][1].upper().encode('ascii')).hexdigest()
                if fields[1].lower() != expected_md5:
                    raise ValueError(f'protein MD5 mismatch: {identifier}; use a scan of this exact proteome')
                if not 1 <= int(fields[6]) <= int(fields[7]) <= int(fields[2]):
                    raise ValueError('invalid domain coordinates')
                if fields[9] != 'T':
                    continue
                identifiers = [fields[4]]
                if len(fields) >= 13 and fields[11] != '-':
                    identifiers.append(fields[11])
                domains.setdefault(identifier, set()).update(
                    normalize_domain(x) for x in identifiers if x and x != '-')
            except (ValueError, IndexError) as error:
                raise ValueError(f'{path}:{number}: invalid InterProScan TSV: {error}') from error
    return domains

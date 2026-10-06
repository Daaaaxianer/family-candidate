#!/usr/bin/env python3
"""Evidence-led family candidate search. Python 3.9+; no Python dependencies."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone

from family_core import (read_fasta, write_fasta, read_domtbl, read_interpro,
                         positive_float, best_domains, domain_records, passes,
                         normalize_domain)


def fraction(value):
    value = float(value)
    if not 0 <= value <= 1:
        raise argparse.ArgumentTypeError('must be between 0 and 1')
    return value


def positive_int(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError('must be at least 1')
    return value


def threshold(value):
    try:
        return positive_float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hmm', required=True, type=Path)
    p.add_argument('--proteins', required=True, type=Path)
    p.add_argument('--cds', type=Path, help='optional CDS FASTA; exact IDs unless --id-map is given')
    p.add_argument('--outdir', required=True, type=Path, help='new or empty directory')
    p.add_argument('--sequence-evalue', type=threshold, default=1e-10)
    p.add_argument('--domain-evalue', type=threshold, default=1e-10)
    p.add_argument('--min-hmm-coverage', type=fraction, default=0.0)
    p.add_argument('--cpu', type=positive_int, default=1)
    p.add_argument('--domtbl', type=Path, help='import existing original-model hmmsearch domtblout')
    p.add_argument('--refine', action='store_true', help='also search an optional species HMM')
    p.add_argument('--seed-evalue', type=threshold, default=1e-20)
    p.add_argument('--seed-domain-evalue', type=threshold, default=1e-20)
    p.add_argument('--seed-min-coverage', type=fraction, default=0.6)
    p.add_argument('--min-seeds', type=positive_int, default=3)
    p.add_argument('--reference-proteins', type=Path, help='optional validated references for BLASTP')
    p.add_argument('--blast-tsv', type=Path, help='import BLASTP TSV with the documented 12 columns')
    p.add_argument('--blast-evalue', type=threshold, default=1e-5)
    p.add_argument('--blast-min-query-coverage', type=fraction, default=0.5)
    p.add_argument('--interpro-tsv', type=Path, help='independent InterProScan TSV for these proteins')
    p.add_argument('--required-domain', action='append', default=[], help='repeat for ALL required signature/IPR IDs')
    p.add_argument('--forbidden-domain', action='append', default=[])
    p.add_argument('--id-map', type=Path, help='TSV: protein_id, gene_id, cds_id (header required)')
    p.add_argument('--hmmsearch-exe', default='hmmsearch')
    p.add_argument('--hmmbuild-exe', default='hmmbuild')
    p.add_argument('--muscle-exe', default='muscle', help='MUSCLE v5')
    p.add_argument('--blastp-exe', default='blastp')
    p.add_argument('--makeblastdb-exe', default='makeblastdb')
    return p


def table(path, columns, rows):
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def run_command(command, log, manifest):
    manifest['commands'].append([str(x) for x in command])
    with open(log, 'w', encoding='utf-8') as handle:
        result = subprocess.run([str(x) for x in command], stdout=handle,
                                stderr=subprocess.STDOUT, check=False)
    if result.returncode:
        raise RuntimeError(f'{command[0]} exited {result.returncode}; see {log}')


def require_output(path):
    if not path.is_file():
        raise RuntimeError(f'command did not create expected output: {path}')


def mapping(path, proteins):
    if not path:
        return {identifier: ('', identifier) for identifier in proteins}
    result = {}
    with open(path, encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'protein_id', 'gene_id', 'cds_id'} <= set(reader.fieldnames or []):
            raise ValueError('ID map requires protein_id, gene_id, cds_id columns')
        for row in reader:
            identifier = row['protein_id']
            if identifier not in proteins or identifier in result:
                raise ValueError(f'unknown or duplicate mapping protein: {identifier}')
            if not row['gene_id']:
                raise ValueError(f'empty mapped gene_id: {identifier}')
            result[identifier] = (row['gene_id'], row['cds_id'])
    missing = set(proteins) - result.keys()
    if missing:
        raise ValueError(f'ID map must cover the input proteome; missing {len(missing)} IDs')
    return result


def read_blast(path, proteins, args):
    # qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen
    from math import isfinite
    best = {}
    with open(path, encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.split()
            try:
                if len(fields) != 12:
                    raise ValueError('expected exactly 12 BLAST columns')
                identifier = fields[1]
                if identifier not in proteins:
                    raise ValueError(f'unknown target {identifier}')
                qstart, qend, sstart, send = map(int, fields[4:8])
                evalue, score = map(float, fields[8:10])
                qlen, slen = map(int, fields[10:12])
                if not (1 <= qstart <= qend <= qlen and 1 <= sstart <= send <= slen):
                    raise ValueError('invalid BLAST coordinates')
                if slen != len(proteins[identifier][1]):
                    raise ValueError(f'target length mismatch: {identifier}')
                if not isfinite(evalue) or evalue < 0 or not isfinite(score):
                    raise ValueError('invalid BLAST score/E-value')
                coverage = (qend-qstart+1)/qlen
                if evalue > args.blast_evalue or coverage < args.blast_min_query_coverage:
                    continue
                value = (evalue, -score, fields[0], coverage)
                if identifier not in best or value[:2] < best[identifier][:2]:
                    best[identifier] = value
            except (ValueError, IndexError) as error:
                raise ValueError(f'{path}:{number}: invalid BLAST TSV: {error}') from error
    return best


def execute(args):
    if args.reference_proteins and args.blast_tsv:
        raise ValueError('choose --reference-proteins or --blast-tsv')
    for key in ('hmm', 'proteins', 'cds', 'domtbl', 'reference_proteins',
                'blast_tsv', 'interpro_tsv', 'id_map'):
        value = getattr(args, key)
        if value and not value.is_file():
            raise ValueError(f'missing input: {value}')
    proteins = read_fasta(args.proteins)
    cds = read_fasta(args.cds) if args.cds else {}
    id_map = mapping(args.id_map, proteins)
    hmm_text = args.hmm.read_text(encoding='utf-8')
    if not hmm_text.startswith('HMMER3/') or '\nNAME ' not in hmm_text:
        raise ValueError('expected an HMMER3 profile file')
    if args.refine and sum(line.startswith('NAME ') for line in hmm_text.splitlines()) != 1:
        raise ValueError('--refine supports one model at a time; split multi-model input')
    model_lengths = {}
    for block in hmm_text.split('//'):
        metadata = {}
        for line in block.splitlines():
            fields = line.split()
            if len(fields) >= 2 and fields[0] in ('NAME', 'ACC', 'LENG'):
                metadata[fields[0]] = fields[1]
        if 'NAME' in metadata:
            if 'LENG' not in metadata or int(metadata['LENG']) <= 0:
                raise ValueError('HMM model has no valid LENG field')
            model_lengths[metadata['NAME']] = int(metadata['LENG'])
    executables = []
    if not args.domtbl or args.refine:
        executables.append(args.hmmsearch_exe)
    if args.reference_proteins:
        read_fasta(args.reference_proteins)
        executables += [args.blastp_exe, args.makeblastdb_exe]
    for executable in executables:
        if not shutil.which(executable):
            raise ValueError(f'executable not found: {executable}')
    out = args.outdir
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError(f'output directory must be new or empty: {out}')
    out.mkdir(parents=True, exist_ok=True)
    manifest = {'status': 'running', 'started_utc': datetime.now(timezone.utc).isoformat(),
                'python': sys.version, 'parameters': {k: str(v) if isinstance(v, Path) else v
                                                     for k, v in vars(args).items()},
                'inputs': {}, 'commands': [], 'warnings': [], 'tools': {}}
    def save_manifest():
        (out/'run.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding='utf-8')
    try:
        for key in ('hmm', 'proteins', 'cds', 'domtbl', 'reference_proteins',
                    'blast_tsv', 'interpro_tsv', 'id_map'):
            value = getattr(args, key)
            if value:
                manifest['inputs'][key] = {'path': str(value.resolve()), 'sha256': sha256(value)}
        def tool_version(executable, option):
            if executable in manifest['tools']:
                return
            result = subprocess.run([executable, option], capture_output=True, text=True,
                                    errors='replace', timeout=20)
            manifest['tools'][executable] = {'path': shutil.which(executable),
                'version_output': (result.stdout + result.stderr)[:2000], 'exit': result.returncode}
        for executable in executables:
            tool_version(executable, '-version' if executable in
                         (args.blastp_exe, args.makeblastdb_exe) else '-h')
        original = out/'original.domtblout'
        if args.domtbl:
            shutil.copyfile(args.domtbl, original)
        else:
            run_command([args.hmmsearch_exe, '--cpu', args.cpu, '--domtblout', original,
                         '-E', max(10, args.sequence_evalue, args.seed_evalue),
                         '--domE', max(10, args.domain_evalue, args.seed_domain_evalue),
                         args.hmm.resolve(), args.proteins.resolve()], out/'original.log', manifest)
            require_output(original)
        hits = read_domtbl(original)
        for hit in hits:
            if hit.target not in proteins or hit.target_length != len(proteins[hit.target][1]):
                raise ValueError(f'domtblout does not match input proteome: {hit.target}')
            if model_lengths.get(hit.model) != hit.model_length:
                raise ValueError(f'domtblout does not match original HMM model: {hit.model}')
        if args.refine:
            seeds = best_domains(hits, args.seed_evalue, args.seed_domain_evalue, args.seed_min_coverage)
            manifest['seed_count'] = len(seeds)
            write_fasta(out/'seed_domains.fasta', domain_records(seeds, proteins))
            if len(seeds) < args.min_seeds:
                manifest['warnings'].append('refinement skipped: too few qualifying seed proteins')
            else:
                for executable in (args.muscle_exe, args.hmmbuild_exe):
                    if not shutil.which(executable):
                        raise ValueError(f'executable not found: {executable}')
                    tool_version(executable, '-version' if executable == args.muscle_exe else '-h')
                alignment, profile = out/'seed_alignment.fasta', out/'species.hmm'
                run_command([args.muscle_exe, '-align', out/'seed_domains.fasta', '-output', alignment],
                            out/'muscle.log', manifest)
                require_output(alignment)
                aligned = read_fasta(alignment)
                expected = dict(domain_records(seeds, proteins))
                if (set(aligned) != set(expected) or len({len(x[1]) for x in aligned.values()}) != 1
                    or any(seq.replace('-', '').replace('.', '').upper() != expected[key].upper()
                           for key, (_, seq) in aligned.items())):
                    raise ValueError('alignment IDs, lengths or ungapped sequences do not match seeds')
                run_command([args.hmmbuild_exe, '--amino', profile, alignment], out/'hmmbuild.log', manifest)
                require_output(profile)
                refined = out/'species.domtblout'
                run_command([args.hmmsearch_exe, '--cpu', args.cpu, '--domtblout', refined,
                             '-E', max(10, args.sequence_evalue), '--domE', max(10, args.domain_evalue),
                             profile, args.proteins.resolve()], out/'species.log', manifest)
                require_output(refined)
                hits.extend(read_domtbl(refined, 'species_hmm'))
        # Also validates lengths/coordinates of any second-round hits.
        for hit in hits:
            if hit.target not in proteins or hit.target_length != len(proteins[hit.target][1]):
                raise ValueError(f'domtblout does not match input proteome: {hit.target}')
        blast_path = args.blast_tsv
        if args.reference_proteins:
            database, blast_path = out/'protein_db', out/'blast.tsv'
            run_command([args.makeblastdb_exe, '-in', args.proteins.resolve(), '-dbtype', 'prot',
                         '-out', database], out/'makeblastdb.log', manifest)
            run_command([args.blastp_exe, '-query', args.reference_proteins.resolve(), '-db', database,
                         '-evalue', args.blast_evalue, '-num_threads', args.cpu,
                         '-max_target_seqs', len(proteins), '-outfmt',
                         '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen',
                         '-out', blast_path], out/'blast.log', manifest)
            require_output(blast_path)
        if args.blast_tsv:
            shutil.copyfile(args.blast_tsv, out/'blast.tsv')
        blast = read_blast(blast_path, proteins, args) if blast_path else {}
        validation = read_interpro(args.interpro_tsv, proteins) if args.interpro_tsv else {}
        if args.interpro_tsv:
            shutil.copyfile(args.interpro_tsv, out/'interpro.tsv')
        required = {normalize_domain(x) for x in args.required_domain}
        forbidden = {normalize_domain(x) for x in args.forbidden_domain}
        if required & forbidden:
            raise ValueError('a domain cannot be both required and forbidden')
        if not args.interpro_tsv or not required:
            manifest['warnings'].append('no complete validation rule/input: candidates remain review')
        grouped = {}
        domain_rows = []
        for hit in hits:
            grouped.setdefault(hit.target, []).append(hit)
            domain_rows.append({'protein_id': hit.target, 'source': hit.source, 'model_id': hit.model_id,
                'sequence_evalue': hit.sequence_evalue, 'domain_ievalue': hit.independent_evalue,
                'domain_score': hit.score, 'hmm_coverage': hit.coverage,
                'domain_start': hit.start, 'domain_end': hit.end})
        domain_columns = ['protein_id','source','model_id','sequence_evalue','domain_ievalue',
                          'domain_score','hmm_coverage','domain_start','domain_end']
        table(out/'domain_evidence.tsv', domain_columns, domain_rows)
        evidence = []
        for identifier in sorted(set(grouped) | set(blast)):
            protein_hits = grouped.get(identifier, [])
            good = [h for h in protein_hits if passes(h, args.sequence_evalue, args.domain_evalue,
                                                      args.min_hmm_coverage)]
            significant = [h for h in protein_hits if passes(h, args.sequence_evalue, args.domain_evalue, 0)]
            selected = min(good or significant or protein_hits,
                           key=lambda h: (h.independent_evalue, -h.score, -h.coverage)) if protein_hits else None
            reasons = []
            if not significant and identifier not in blast:
                status, reasons = 'rejected', ['search_threshold_not_met']
            elif significant and not good and identifier not in blast:
                status, reasons = 'review', ['hmm_coverage_below_threshold']
            else:
                status = 'review'
            observed = validation.get(identifier, set())
            if status != 'rejected':
                if not args.interpro_tsv:
                    reasons.append('domain_validation_not_run')
                elif not required:
                    reasons.append('required_domain_rule_not_set')
                elif not observed:
                    reasons.append('no_validation_record')
                elif observed & forbidden:
                    status = 'rejected'
                    reasons.append('forbidden_domain:' + ','.join(sorted(observed & forbidden)))
                elif not required <= observed:
                    status = 'rejected'
                    reasons.append('required_domain_missing:' + ','.join(sorted(required-observed)))
                elif reasons:
                    pass  # partial HMM evidence still requires manual review
                else:
                    status = 'accepted'
                    reasons.append('search_and_domain_rules_passed')
            gene_id, cds_id = id_map[identifier]
            row = {'protein_id': identifier, 'gene_id': gene_id, 'cds_id': cds_id if args.cds else '',
                   'source': '+'.join(sorted({h.source for h in protein_hits} |
                                            ({'blast'} if identifier in blast else set()))),
                   'best_domain_source': selected.source if selected else '',
                   'domains': ';'.join(sorted(observed)), 'status': status, 'reason': ';'.join(reasons),
                   'blast_query': blast[identifier][2] if identifier in blast else '',
                   'blast_evalue': blast[identifier][0] if identifier in blast else '',
                   'blast_query_coverage': blast[identifier][3] if identifier in blast else '',
                   'cds_found': str(cds_id in cds).lower() if args.cds else '', 'gene_representative': ''}
            row.update({key: '' for key in domain_columns if key not in ('protein_id','source')})
            if selected:
                row.update(model_id=selected.model_id, sequence_evalue=selected.sequence_evalue,
                           domain_ievalue=selected.independent_evalue, domain_score=selected.score,
                           hmm_coverage=selected.coverage, domain_start=selected.start, domain_end=selected.end)
            evidence.append(row)
        # Gene representatives are separate; never merge distinct loci by sequence identity.
        representatives = {}
        if args.id_map:
            for row in evidence:
                if row['status'] == 'rejected':
                    continue
                quality = (row['status'] == 'accepted', len(proteins[row['protein_id']][1]))
                old = representatives.get(row['gene_id'])
                if old is None or quality > old[0]:
                    representatives[row['gene_id']] = (quality, row['protein_id'])
            representative_ids = {x[1] for x in representatives.values()}
            for row in evidence:
                row['gene_representative'] = str(row['protein_id'] in representative_ids).lower()
            write_fasta(out/'gene_representatives.protein.fasta', [proteins[x] for x in sorted(representative_ids)])
        columns = ['protein_id','gene_id','cds_id','source','best_domain_source'] + domain_columns[2:] + [
                   'blast_query','blast_evalue','blast_query_coverage','domains','status','reason',
                   'cds_found','gene_representative']
        table(out/'candidate_evidence.tsv', columns, evidence)
        missing_cds = []
        for label, predicate in [('candidates', lambda r: r['status'] != 'rejected'),
                                 ('accepted', lambda r: r['status'] == 'accepted'),
                                 ('review', lambda r: r['status'] == 'review'),
                                 ('rejected', lambda r: r['status'] == 'rejected')]:
            rows = [row for row in evidence if predicate(row)]
            (out/(label+'.ids.txt')).write_text(''.join(r['protein_id']+'\n' for r in rows), encoding='utf-8')
            write_fasta(out/(label+'.protein.fasta'), [proteins[r['protein_id']] for r in rows])
            if args.cds:
                # Preserve actual CDS IDs; their relation is explicit in the evidence table.
                keys = list(dict.fromkeys(r['cds_id'] for r in rows if r['cds_id'] in cds))
                write_fasta(out/(label+'.cds.fasta'), [cds[x] for x in keys])
        if args.cds:
            missing_cds = [r for r in evidence if r['status'] != 'rejected' and r['cds_id'] not in cds]
            table(out/'missing_cds.tsv', ['protein_id','cds_id'],
                  [{k:r[k] for k in ('protein_id','cds_id')} for r in missing_cds])
            if missing_cds:
                manifest['warnings'].append(f'{len(missing_cds)} candidate proteins have no mapped CDS')
        manifest['counts'] = {label: sum(r['status'] == label for r in evidence)
                              for label in ('accepted','review','rejected')}
        manifest['counts']['candidate_proteins'] = sum(r['status'] != 'rejected' for r in evidence)
        if args.id_map:
            manifest['counts']['candidate_gene_loci'] = len(representatives)
        manifest['status'] = 'success'
        manifest['finished_utc'] = datetime.now(timezone.utc).isoformat()
        save_manifest()
        print(json.dumps({'outdir':str(out.resolve()), 'counts':manifest['counts'],
                          'warnings':manifest['warnings']}, ensure_ascii=False))
    except Exception as error:
        manifest['status'], manifest['error'] = 'failed', str(error)
        manifest['finished_utc'] = datetime.now(timezone.utc).isoformat()
        save_manifest()
        raise


def main():
    args = parser().parse_args()
    try:
        execute(args)
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())

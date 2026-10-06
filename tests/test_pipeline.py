"""Regression + offline integration tests. External tool outputs are explicitly simulated."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import family_candidate as pipeline
from family_core import (read_fasta, read_domtbl, best_domains, domain_records,
                         extract_exact, read_interpro, read_ids)


def domain(target='p1', model='model', accession='PF00001.1', length=20,
           seq_e='1e-30', dom_e='1e-25', score=80, start=1, end=10,
           hmm_start=1, hmm_end=10, model_length=10, number=1):
    return ' '.join(map(str, [target,'-',length,model,accession,model_length,
        seq_e,100,0,number,2,dom_e,dom_e,score,0,hmm_start,hmm_end,start,end,start,end,0.9,'example']))


def interpro(identifier, sequence, signature='PF00001', ipr='IPR000001'):
    checksum = hashlib.md5(sequence.upper().encode()).hexdigest()
    return '\t'.join(map(str, [identifier,checksum,len(sequence),'Pfam',signature,'domain',
                              1,len(sequence),'1e-20','T','06-10-2026',ipr,'description']))


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='family test ')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        quiet = patch.object(pipeline, 'print', create=True)
        quiet.start()
        self.addCleanup(quiet.stop)
        self.pro = self.file('proteins.fa', '>p1\nAAAACCCCCCGGGGGGGGGG\n>p2\nMMMMMMMMMMMMMMMMMMMM\n')
        self.hmm = self.file('model.hmm', 'HMMER3/f\nNAME  model\nACC   PF00001.1\nLENG  10\n//\n')
        self.dom = self.file('original.txt', domain()+'\n')

    def file(self, name, content):
        path = self.root/name
        path.write_text(content, encoding='utf-8')
        return path

    def args(self, *extra):
        return pipeline.parser().parse_args(['--hmm',str(self.hmm), '--proteins',str(self.pro),
            '--domtbl',str(self.dom), '--outdir',str(self.root/'result'), *map(str, extra)])

    def execute(self, *extra):
        pipeline.execute(self.args(*extra))
        with open(self.root/'result/candidate_evidence.tsv', newline='', encoding='utf-8') as handle:
            return list(csv.DictReader(handle, delimiter='\t'))

    def test_exact_id_not_substring(self):
        proteins = read_fasta(self.file('match.fa','>gene1\nAAAA\n>gene10\nCCCC\n>other_gene1\nGGGG\n'))
        self.assertEqual(extract_exact(['gene1'], proteins), [('gene1','AAAA')])

    def test_regex_characters_literal(self):
        proteins = read_fasta(self.file('match.fa','>gene.1\nAAAA\n>geneX1\nCCCC\n'))
        self.assertEqual(extract_exact(['gene.1'], proteins), [('gene.1','AAAA')])

    def test_missing_id_fails(self):
        with self.assertRaisesRegex(ValueError, 'IDs not found'):
            extract_exact(['missing'], read_fasta(self.pro))

    def test_duplicate_fasta_id_fails(self):
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            read_fasta(self.file('duplicate.fa','>p1\nAA\n>p1\nCC\n'))

    def test_blank_duplicate_requested_ids(self):
        self.assertEqual(read_ids(self.file('ids','p2\n\np1\np2\n')), ['p2','p1'])

    def test_empty_sequence_fails(self):
        with self.assertRaisesRegex(ValueError, 'empty sequence'):
            read_fasta(self.file('empty.fa','>p1\n'))

    def test_best_domain_not_first(self):
        path = self.file('domains', domain(dom_e='0.1',score=1,start=1,end=4)+'\n'+
                         domain(start=11,end=20,number=2)+'\n')
        hits = best_domains(read_domtbl(path),1e-20,1e-20,0)
        self.assertEqual(len(hits),1)
        self.assertEqual(domain_records(hits,read_fasta(self.pro))[0][1], 'GGGGGGGGGG')

    def test_domain_order_invariance(self):
        first, second = domain(dom_e='1e-22'), domain(dom_e='1e-25',start=11,end=20)
        a = best_domains(read_domtbl(self.file('a',first+'\n'+second)),1e-20,1e-20,0)
        b = best_domains(read_domtbl(self.file('b',second+'\n'+first)),1e-20,1e-20,0)
        self.assertEqual(a,b)

    def test_multiple_models_kept_separate(self):
        path = self.file('multi',domain()+'\n'+domain(model='other',accession='PF00002.1'))
        self.assertEqual(len(best_domains(read_domtbl(path),1e-20,1e-20,0)),2)

    def test_threshold_boundary_and_zero(self):
        path = self.file('edge',domain(seq_e='1e-20',dom_e='1e-20')+'\n'+domain(target='p2',seq_e=0,dom_e=0))
        self.assertEqual(len(best_domains(read_domtbl(path),1e-20,1e-20,1)),2)

    def test_bad_domtbl_columns_fail(self):
        with self.assertRaisesRegex(ValueError, 'invalid domtblout'):
            read_domtbl(self.file('bad','p1 not enough columns\n'))

    def test_bad_coordinates_fail(self):
        with self.assertRaisesRegex(ValueError, 'coordinates'):
            read_domtbl(self.file('bad',domain(end=21)))

    def test_unvalidated_remains_review(self):
        rows = self.execute()
        self.assertEqual(rows[0]['status'],'review')
        self.assertEqual((self.root/'result/accepted.ids.txt').read_text(),'')

    def test_final_threshold_changes_members(self):
        self.dom.write_text(domain(seq_e='1e-8',dom_e='1e-8'))
        rows = self.execute('--sequence-evalue','1e-10','--domain-evalue','1e-10')
        self.assertEqual(rows[0]['status'],'rejected')
        self.assertEqual(rows[0]['source'],'pfam_hmm')
        other = self.root/'looser'
        args = self.args('--sequence-evalue','1e-5','--domain-evalue','1e-5')
        args.outdir = other
        pipeline.execute(args)
        self.assertEqual((other/'candidates.ids.txt').read_text(),'p1\n')

    def test_validation_accepts(self):
        scan = self.file('interpro.tsv',interpro('p1',read_fasta(self.pro)['p1'][1]))
        rows = self.execute('--interpro-tsv',scan,'--required-domain','PF00001.99')
        self.assertEqual(rows[0]['status'],'accepted')

    def test_validation_rule_missing(self):
        scan = self.file('interpro.tsv',interpro('p1',read_fasta(self.pro)['p1'][1]))
        self.assertEqual(self.execute('--interpro-tsv',scan)[0]['status'],'review')

    def test_required_and_forbidden_domains(self):
        scan = self.file('interpro.tsv',interpro('p1',read_fasta(self.pro)['p1'][1]))
        rows = self.execute('--interpro-tsv',scan,'--required-domain','PF00001',
                            '--forbidden-domain','IPR000001')
        self.assertEqual(rows[0]['status'],'rejected')
        self.assertIn('forbidden_domain',rows[0]['reason'])

    def test_required_domain_missing(self):
        scan = self.file('interpro.tsv',interpro('p1',read_fasta(self.pro)['p1'][1],signature='PF99999'))
        self.assertEqual(self.execute('--interpro-tsv',scan,'--required-domain','PF00001')[0]['status'],'rejected')

    def test_validation_md5_mismatch_fails(self):
        scan = self.file('interpro.tsv',interpro('p1','Z'*20))
        with self.assertRaisesRegex(ValueError,'MD5 mismatch'):
            read_interpro(scan,read_fasta(self.pro))

    def test_partial_hit_remains_review(self):
        self.dom.write_text(domain(hmm_end=4,end=4))
        rows = self.execute('--min-hmm-coverage','0.6')
        self.assertEqual(rows[0]['status'],'review')
        self.assertIn('coverage',rows[0]['reason'])

    def test_blast_only_candidate_added(self):
        blast = self.file('blast.tsv','ref\tp2\t80\t20\t1\t20\t1\t20\t1e-20\t100\t20\t20\n')
        rows = self.execute('--blast-tsv',blast)
        self.assertEqual([row['protein_id'] for row in rows],['p1','p2'])
        self.assertEqual(rows[1]['source'],'blast')
        self.assertEqual(rows[1]['status'],'review')

    def test_missing_cds_reported(self):
        cds = self.file('cds.fa','>p2\nATG\n')
        rows = self.execute('--cds',cds)
        self.assertEqual(rows[0]['cds_found'],'false')
        self.assertIn('p1', (self.root/'result/missing_cds.tsv').read_text())

    def test_no_hits_is_successful_empty_result(self):
        self.dom.write_text('# no hits\n\n')
        self.assertEqual(self.execute(),[])
        self.assertEqual(json.loads((self.root/'result/run.json').read_text())['status'],'success')

    def test_nonempty_output_refused(self):
        (self.root/'result').mkdir()
        sentinel = self.file('result/sentinel','do not overwrite')
        with self.assertRaisesRegex(ValueError,'new or empty'):
            pipeline.execute(self.args())
        self.assertEqual(sentinel.read_text(),'do not overwrite')

    def test_mismatched_model_rejected(self):
        self.dom.write_text(domain(model='wrong'))
        with self.assertRaisesRegex(ValueError,'original HMM model'):
            pipeline.execute(self.args())

    def test_gene_mapping_representative_and_cds(self):
        self.dom.write_text(domain()+'\n'+domain(target='p2'))
        idmap = self.file('map.tsv','protein_id\tgene_id\tcds_id\np1\tg1\tc1\np2\tg1\tc2\n')
        cds = self.file('cds.fa','>c1\nATG\n>c2\nATGATG\n')
        rows = self.execute('--id-map',idmap,'--cds',cds)
        self.assertEqual(sum(row['gene_representative']=='true' for row in rows),1)
        self.assertEqual(set(read_fasta(self.root/'result/candidates.cds.fasta')),{'c1','c2'})
        self.assertEqual(json.loads((self.root/'result/run.json').read_text())['counts']['candidate_gene_loci'],1)

    def test_cli_works_from_different_directory_with_spaces(self):
        args = self.args()
        result = subprocess.run([sys.executable,str(ROOT/'family_candidate.py'), '--hmm',str(args.hmm),
            '--proteins',str(args.proteins),'--domtbl',str(args.domtbl), '--outdir',str(args.outdir)],
            cwd=self.root,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)

    def test_external_failure_manifest(self):
        args = self.args('--hmmsearch-exe',sys.executable)
        args.domtbl = None
        with self.assertRaisesRegex(RuntimeError,'exited'):
            pipeline.execute(args)  # Python cannot act as hmmsearch: deliberate command failure.
        manifest = json.loads((self.root/'result/run.json').read_text())
        self.assertEqual(manifest['status'],'failed')
        self.assertFalse((self.root/'result/candidates.ids.txt').exists())

    def test_refinement_union_with_simulated_external_tools(self):
        self.dom.write_text(domain()+'\n'+domain(target='p2'))
        args = self.args('--refine','--min-seeds','2','--hmmsearch-exe',sys.executable,
                         '--muscle-exe',sys.executable,'--hmmbuild-exe',sys.executable)
        def simulate(command, log, manifest):
            manifest['commands'].append([str(x) for x in command])
            Path(log).write_text('SIMULATED EXTERNAL TOOL OUTPUT\n')
            if '-align' in command:
                Path(command[command.index('-output')+1]).write_text(
                    Path(command[command.index('-align')+1]).read_text())
            elif '--amino' in command:
                Path(command[2]).write_text('SIMULATED species HMM\n')
            elif '--domtblout' in command:
                # Only p2 remains in the species-model output; p1 must survive via original HMM.
                Path(command[command.index('--domtblout')+1]).write_text(
                    domain(target='p2',model='species',accession='-'))
            else:
                self.fail('unexpected simulated command')
        with patch.object(pipeline,'run_command',side_effect=simulate):
            pipeline.execute(args)
        self.assertEqual((self.root/'result/candidates.ids.txt').read_text(),'p1\np2\n')
        self.assertEqual(json.loads((self.root/'result/run.json').read_text())['seed_count'],2)

    def test_seed_shortage_skips_refinement(self):
        args = self.args('--refine','--hmmsearch-exe',sys.executable)
        pipeline.execute(args)
        manifest = json.loads((self.root/'result/run.json').read_text())
        self.assertTrue(any('too few' in x for x in manifest['warnings']))
        self.assertFalse((self.root/'result/species.hmm').exists())

    def test_helper_missing_id_does_not_write(self):
        ids = self.file('ids.txt','missing\n')
        output = self.root/'extract.fa'
        result = subprocess.run([sys.executable,str(ROOT/'retrieve.seq.from.all.fasta.py'),
                                 str(ids),str(self.pro),str(output)],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertFalse(output.exists())

    def test_blast_only_validation_can_accept(self):
        blast = self.file('blast.tsv','ref\tp2\t80\t20\t1\t20\t1\t20\t1e-20\t100\t20\t20\n')
        scan = self.file('interpro.tsv',interpro('p2',read_fasta(self.pro)['p2'][1]))
        rows = self.execute('--blast-tsv',blast,'--interpro-tsv',scan,'--required-domain','PF00001')
        self.assertEqual(rows[1]['status'],'accepted')

    def test_missing_validation_record_is_review(self):
        scan = self.file('interpro.tsv',interpro('p2',read_fasta(self.pro)['p2'][1]))
        rows = self.execute('--interpro-tsv',scan,'--required-domain','PF00001')
        self.assertEqual(rows[0]['status'],'review')
        self.assertEqual(rows[0]['reason'],'no_validation_record')

    def test_failed_pipeline_records_error(self):
        self.dom.write_text('malformed domtblout')
        with self.assertRaises(ValueError):
            pipeline.execute(self.args())
        manifest = json.loads((self.root/'result/run.json').read_text())
        self.assertEqual(manifest['status'],'failed')
        self.assertIn('invalid domtblout',manifest['error'])

    def test_blank_domtbl_lines_safe(self):
        self.assertEqual(read_domtbl(self.file('blank','\n  \n# comment\n')),[])

    def test_cli_invalid_thresholds_rejected(self):
        for value in ['0','-1','nan','inf']:
            result = subprocess.run([sys.executable,str(ROOT/'family_candidate.py'),
                '--hmm',str(self.hmm),'--proteins',str(self.pro),'--domtbl',str(self.dom),
                '--outdir',str(self.root/'unused'),'--sequence-evalue',value],
                capture_output=True,text=True)
            self.assertEqual(result.returncode,2,value)

    def test_refinement_alignment_failure_stops_search(self):
        self.dom.write_text(domain()+'\n'+domain(target='p2'))
        args = self.args('--refine','--min-seeds','2','--hmmsearch-exe',sys.executable,
                         '--muscle-exe',sys.executable,'--hmmbuild-exe',sys.executable)
        def invalid_alignment(command, log, manifest):
            Path(command[command.index('-output')+1]).write_text('>wrong_id\nAAAA\n')
        with patch.object(pipeline,'run_command',side_effect=invalid_alignment):
            with self.assertRaisesRegex(ValueError,'alignment IDs'):
                pipeline.execute(args)
        self.assertFalse((self.root/'result/species.hmm').exists())

    def test_missing_executable_does_not_create_results(self):
        args = self.args('--hmmsearch-exe','definitely_no_such_executable_20261006')
        args.domtbl = None
        with self.assertRaisesRegex(ValueError,'executable not found'):
            pipeline.execute(args)
        self.assertFalse(args.outdir.exists())

    def test_blast_command_with_simulated_tools(self):
        args = self.args('--reference-proteins',self.pro, '--blastp-exe',sys.executable,
                         '--makeblastdb-exe',sys.executable)
        commands = []
        def simulate(command, log, manifest):
            commands.append(command)
            if '-query' in command:
                Path(command[command.index('-out')+1]).write_text(
                    'ref\tp2\t80\t20\t1\t20\t1\t20\t1e-20\t100\t20\t20\n')
        with patch.object(pipeline,'run_command',side_effect=simulate):
            pipeline.execute(args)
        self.assertEqual(len(commands),2)
        self.assertIn('-max_target_seqs',commands[1])
        self.assertEqual((self.root/'result/candidates.ids.txt').read_text(),'p1\np2\n')


if __name__ == '__main__':
    unittest.main()

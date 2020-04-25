# Simply introduction

Use Pfam Hmm Seed to Search Gene Family Candidates
1. Identification of family seed candidates with high confidence using Pfam HMM seeds --- Strict threshold
2. Generating HMM seeds of specific species
3. Identification of family candidates by specific HMM seeds --- threshold

# Dependency package

1.Perl https://www.perl.org/

2.Bioperl https://bioperl.org/

3.Muscle http://drive5.com/muscle/

4.Hmmer http://hmmer.org/download.html

# How to use

$ sh batch_pfam_hmmer_to_family_candidate_search.sh <hmm.seed> <protein.fasta> <Strict threshold> <threshold> <out.name>

example:

$ sh batch_pfam_hmmer_to_family_candidate_search.sh NB-ARC.hmm Athaliana_167_TAIR10.protein.fa 1e-20 1e-10 NB-ARC.in.At

# Result file

1.NB-ARC.in.At.id.txt

The id list of searched family candidate

2.NB-ARC.in.At.protein.fasta

The file contains the protein sequence, and the sequence name is consistent with the ID file

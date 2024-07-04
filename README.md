# Simply introduction

## First, Use Pfam Hmm Seed to Search Gene Family Candidates 
1. Identification of family seed candidates with high confidence using Pfam HMM seeds --- Strict threshold
2. Generating HMM seeds of specific species
3. Identification of family candidates by specific HMM seeds --- threshold

**Provides running methods for perl and python versions respectively**

## 1. perl version
### Dependency package

1.Perl https://www.perl.org/

2.Bioperl https://bioperl.org/

3.Muscle3 http://drive5.com/muscle/

4.Hmmer3 http://hmmer.org/download.html

###  How to use
```Bash
sh batch_pfam_hmmer_to_family_candidate_search.perl.sh hmm.seed protein.fasta cds.fasta Strict_threshold threshold outname    
```
Example:
```Bash
sh batch_pfam_hmmer_to_family_candidate_search.perl.sh NB-ARC.hmm Athaliana_167_TAIR10.protein.fa Athaliana_167_TAIR10.cds.fa 1e-20 1e-10 NB-ARC.in.At
```

## 2. python version
### Dependency package

1.Python https://www.python.org/

2.Biopython https://biopython.org/

3.Muscle5 http://drive5.com/muscle/

4.Hmmer3 http://hmmer.org/download.html

###  How to use
```Bash
sh batch_pfam_hmmer_to_family_candidate_search.sh hmm.seed protein.fasta cds.fasta Strict_threshold threshold outname    
```
Example:
```Bash
sh batch_pfam_hmmer_to_family_candidate_search.sh NB-ARC.hmm Athaliana_167_TAIR10.protein.fa Athaliana_167_TAIR10.cds.fa 1e-20 1e-10 NB-ARC.in.At
```

Example in terminal:

![code_example.gif](https://github.com/Daaaaxianer/family-candidate/blob/master/code_example.gif)


## Result file

1.NB-ARC.in.At.id.txt

The id list of searched family candidate

2.NB-ARC.in.At.protein.fasta

The file contains the protein sequence, and the sequence name is consistent with the id file

3.NB-ARC.in.At.cds.fasta

The file contains the cds sequence, and the sequence name is consistent with the id file

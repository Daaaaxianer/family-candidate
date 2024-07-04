#! /bin/bash

### use pfam hmm seeds serach gene family candidates
### hmm_seed : hmm seed file of pfam
### protein.fa : fasta file of protein
### cds.fa : fasta file of cds (nucleotide)
### e-value1 : threshold of e-value to build specie hmm seed
### e-value1 : threshold of e-value to filter family candidate
### outname : outfile name for id file and fasta file like: NBSS-LRR.in.at 
### xian 2020-04-01 ###

if [ $# -lt 5 ];then
    echo "Five parameters needed, but only $# given!"
    echo "Usage: sh $0 <hmm_seed> <protein.fa> <cds.fa> <e-value1> <e-value2> <outname>"
exit 1;
fi

seed=$1
pro=$2
nuc=$3
evalue1=$4
evalue2=$5
name=$6

### search use hmmseed  
hmmsearch --domtblout ${seed}_domain_search.txt ${seed} ${pro}
### refine domian sequence
perl refine_domain_seq.pl ${seed}_domain_search.txt ${pro} ${name}.domain.fasta ${evalue1}
### mutiple alignment for hmmer
muscle -in ${name}.domain.fasta -phyiout ${name}.phyi -clwout ${name}.aln
### build Specificity hmmseed
hmmbuild ${name}.hmm ${name}.phyi
### search use Specificity hmmseed
hmmsearch --domtblout ${name}.serach.txt ${name}.hmm ${pro}
### filter id use awk sort and uniq
awk -F '[ ]+' '{if($1 !~ /#/ && $7 < 1e-10 ) {print $1}}' ${name}.serach.txt |sort |uniq > ${name}.id.txt
### retrieve id sequence from all sequence
perl retrieve.seq.from.all.fasta.pl ${name}.id.txt ${pro} ${name}.protein.fasta
perl retrieve.seq.from.all.fasta.pl ${name}.id.txt ${nuc} ${name}.cds.fasta

###
echo "${name} id file : ${name}.id.txt"
echo "${name} protein file : ${name}.protein.fa"
echo "${name} cds file : ${name}.cds.fa"
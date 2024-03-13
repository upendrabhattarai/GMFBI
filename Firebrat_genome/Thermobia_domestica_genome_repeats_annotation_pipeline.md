This is the pipeline used for repeats annotation in Thermobia domestica genome assembly. We used Haplotype-1 as our reference genome because it was the most complete genome compared to the other haplotype or primary assembly.
## 1. Repeatmodeler
Used `RepeatModeler` `V.2.0.5`
Search Engine = rmblast 2.14.1+
Dependencies: TRF 4.09, RECON 1.08, RepeatScout 1.0.6, RepeatMasker 4.1.5
LTR Structural Analysis: Disabled
Random Number Seed: 1701396861
```
export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software:$PATH

## Step-1 build database
singularity exec tetools_latest.sif BuildDatabase -name T.dom.Hap1.final.fa -engine ncbi ../T.dom.Hap1.final.fa

## Step-2 run repeatmoduler
singularity exec tetools_latest.sif RepeatModeler -engine ncbi -threads 46 -database T.dom.Hap1.final.fa
```

repmodeler.out
```
Number of families returned by RECON: 30831
Processing families with greater than 15 elements
Instance Gathering: 00:00:16 (hh:mm:ss) Elapsed Time
Refining 2357 families
Family Refinement: 00:30:29 (hh:mm:ss) Elapsed Time
Round Time: 11:38:11 (hh:mm:ss) Elapsed Time : 2233 families discovered.

RepeatScout/RECON discovery complete: 4310 families found
```

## 2. TransposonPSI

Used `TransposonPSI` `V.2.2.26`G
```
module load ncf/1.0.0-fasrc01
module load perl_modules/5.10.1-ncf
export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/mambaforge/bin:$PATH
export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/TransposonPSI_08222010:$PATH

# Step 1
transposonPSI.pl ../T.dom.Hap1.final.fa nuc

# Step 2
transposonPSI_2_fasta.pl ../T.dom.Hap1.final.fa T.dom.Hap1.final.fa.TPSI.allHits.chains.bestPerLocus
```

## 3. LTRHarvest
`GenomeTools` `V.1.5.9` was used to run the ltr program

find the rule file `filter_protein_match.lua` [here](filter_protein_match.lua)
```
export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/mambaforge/bin:$PATH

# Run following commands sequentially
# Step-1: create index
gt suffixerator -dna -db ../T.dom.Hap1.final.fa -lcp -ssp -suf -tis -des -lossless

# Step-2: find candidate elements
gt -j 60 ltrharvest -index T.dom.Hap1.final.fa -tabout no -seqids -md5 > T.dom.Hap1.final.fa.ltrh.gff3

# Step-3: quick check of results
gt stat T.dom.final.Pri.fasta.ltrh.gff3 > stat.T.dom.final.Pri.fasta.ltrh.gff3.out

# Step-4: Sort the gff3 file produced
gt gff3 -sort T.dom.Hap1.final.fa.ltrh.gff3 > T.dom.Hap1.final.fa.ltrh.sorted.gff3

# Step-5
#Download the GyDB collection (downloaded December-2, 2023, 09:49) and unzip it and run ltrdigest

gt -j 60 ltrdigest -hmms GyDB_collection/profiles/*.hmm \
                  -encseq T.dom.Hap1.final.fa -matchdescstart \
                  < T.dom.Hap1.final.fa.ltrh.sorted.gff3 > T.dom.Hap1.final.fa.ltrh.sorted.ltrd.gff3

# Step-6
gt select -rule_files filter_protein_match.lua < T.dom.Hap1.final.fa.ltrh.sorted.ltrd.gff3 > T.dom.Hap1.final.fa.ltrh.sorted.ltrd.filtered.gff3

# Step-7
gt stat T.dom.Hap1.final.fa.ltrh.sorted.ltrd.filtered.gff3 >> stat.T.dom.final.Pri.fasta.ltrh.gff3.out 

# Step-8
gt extractfeat -type LTR_retrotransposon -encseq T.dom.Hap1.final.fa \
                -matchdescstart < T.dom.Hap1.final.fa.ltrh.sorted.ltrd.filtered.gff3 \
                > T.dom.Hap1.final.fa.ltrh.sorted.ltrd.filtered.sequences.fasta
```

Stats:
```
parsed genome node DAGs: 105493
sequence regions: 2120 (total length: 5637872114)
LTR_retrotransposons: 101253
long_terminal_repeats: 202506
repeat_regions: 101253
target_site_duplications: 202506
parsed genome node DAGs: 18140
sequence regions: 2120 (total length: 5637872114)
LTR_retrotransposons: 13900
RR_tracts: 3488
long_terminal_repeats: 27800
protein_matchs: 530639
repeat_regions: 13900
target_site_duplications: 27800
```

## 4. Repbase
I extracted Insect library from the Whole Repbasedatabase.
RepBaseRepeatMaskerEdition-20181026 was used as a database
Older version of `RepeatMasker-4.1.0` was used for the `queryRepeatDatabase.pl` script. You can find the script [here](queryRepeatDatabase.pl)

## 5. SINE 
Database was downloaded March 7, 2023 to use with the repeat library

## 6. SRF
We ran Satellite repeat finder (SRF) in the QC'ed pacbio reads.
`KMC` `V.3.2.1` (2022-01-04) and `minimap2` `V.2.26-r1175`was used
```
export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/mambaforge/bin:$PATH

# Step -1
kmc -fq -k151 -t36 -ci200 -cs100000 -m100 @fofn.txt count.kmc tmp_dir

# Step -2
kmc_dump count.kmc count.txt

# Step -3
/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/srf/srf -p srf.T.dom.hifi count.txt > srf.fa

# Step -4
minimap2 -c -N1000000 -f1000 -r100,100 <(./srfutils.js enlong srf.fa) *.fastq.gz > srf-aln.paf

# Step -5
srfutils.js paf2bed srf-aln.paf > srf-aln.bed   # filter and extract non-overlapping regions

# Step -6
#srfutils.js bed2abun srf-aln.bed > srf-aln.len  # calculate abundance of each srf contig
```
## 7. RepeatClassifier
We concatenated all the libraries together and ran RepeatClassifier

```
# Step 1: Filtering to remove sequences less than 50bp
seqtk seq -L 50 Combined.rep.library.fasta > Combined.rep.library.minlen50.fasta

# Step 2: remove redundant sequences using Usearch
usearch -cluster_fast Combined.rep.library.minlen50.fasta -id 0.8 -consout Combined.rep.library.minlen50.usearch.fasta

# Step 3: Ran repeatclassifier
export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software:$PATH

singularity exec tetools_latest.sif RepeatClassifier -consensi Combined.rep.library.minlen50.usearch.fasta

# Step 4: Filtering out Unknown hits
grep "Unknown" Combined.rep.library.minlen50.usearch.fasta.classified > UnknownIDs.txt
sed -i 's/>//g' UnknownIDs.txt
seqtk subseq Combined.rep.library.minlen50.usearch.fasta.classified UnknownIDs.txt > Unknown_repeats.fasta

# Step 5
# Download Uniprot reviewed insect database from Uniprot website. our data downloaded 2023 Dec 03 and blast the unknown repeats against this database
makeblastdb -in uniprotkb_insecta_AND_reviewed_true_2023_12_03.fasta -dbtype prot

blastx -query ../Unknown_repeats.fasta -db  uniprotkb_insecta_AND_reviewed_true_2023_12_03.fasta -evalue 1e-10 \
        -num_threads 36 -max_target_seqs 1 -outfmt '6 qseqid sseqid evalue bitscore sgi sacc stitle' -out Blast_out.txt

# Step 6
# Filter the repeat library from the blastx report

awk -F "\t" '{print $1,$7}' Blast_out.txt  | sort | uniq | grep -i -v "transposon" | grep -i -v "Copia protein" | grep -i -v "mobile element" | \
grep -i -v "transposable"  | grep -i -v "transposase" | awk '{print $1}' > Unknowns_with_Port_hit.txt


faSomeRecords -exclude Combined.rep.library.minlen50.usearch.fasta.classified Unknowns_with_Port_hit.txt Combined.rep.library.minlen50.usearch.fasta.classified.filtered
```
The final Repeat library is `Combined.rep.library.minlen50.usearch.fasta.classified.filtered`
Find `faSomeRecords` script [here](https://github.com/santiagosnchez/faSomeRecords)

## 8. RepeatMasker

We used the finalized repeat library in above step to mask the repeats in our genome assembly with `RepeatMasker` 
```
#Repeatmasker
#export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software:$PATH

## Step-1 build database
#tetools_latest.sif RepeatMasker -pa 64 --xsmall -lib Combined.rep.library.minlen50.usearch.fasta.classified.filtered T.dom.Hap1.final.fa

#To obtain gff file
singularity exec /n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/tetools_latest.sif rmOutToGFF3.pl T.dom.Hap1.final.fa.out T.dom.Hap1.final.fa.gff
```
RepeatMasker.out.tbl
```
file name: T.dom.Hap1.final.fa      
sequences:          2531
total length: 5651003049 bp  (5650743249 bp excl N/X-runs)
GC level:         35.77 %
bases masked: 4181935130 bp ( 74.00 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
Retroelements       5043437   1488194435 bp   26.34 %
   SINEs:           581341    108439969 bp    1.92 %
   Penelope:           713       131670 bp    0.00 %
   LINEs:           2118921    497840224 bp    8.81 %
    CRE/SLACS        23955     11164290 bp    0.20 %
     L2/CR1/Rex     567973    103624063 bp    1.83 %
     R1/LOA/Jockey  288766     82194795 bp    1.45 %
     R2/R4/NeSL       7458      1965453 bp    0.03 %
     RTE/Bov-B      859872    205238867 bp    3.63 %
     L1/CIN4         24663      4061995 bp    0.07 %
   LTR elements:    2343175    881914242 bp   15.61 %
     BEL/Pao        277318     79470783 bp    1.41 %
     Ty1/Copia       75876     22541418 bp    0.40 %
     Gypsy/DIRS1    1778783    754984774 bp   13.36 %
       Retroviral    15924      2105623 bp    0.04 %

DNA transposons     3455801   1029429975 bp   18.22 %
   hobo-Activator   643710    179824938 bp    3.18 %
   Tc1-IS630-Pogo   1835244    492765939 bp    8.72 %
   En-Spm                0            0 bp    0.00 %
   MULE-MuDR         81266     21742295 bp    0.38 %
   PiggyBac         177207     29624824 bp    0.52 %
   Tourist/Harbinger 34046      6909409 bp    0.12 %
   Other (Mirage,    28151      4653715 bp    0.08 %
    P-element, Transib)

Rolling-circles     105015     26706290 bp    0.47 %

Unclassified:       6084928   1535698957 bp   27.18 %

Total interspersed repeats:  4053455037 bp   71.73 %


Small RNA:          161318     57146225 bp    1.01 %

Satellites:          16711      1788949 bp    0.03 %
Simple repeats:     1367386     70771015 bp    1.25 %
Low complexity:      95193      4464186 bp    0.08 %
==================================================

* most repeats fragmented by insertions or deletions
  have been counted as one element
                                                      

RepeatMasker version 4.1.5 , default mode
                                        
run with rmblastn version 2.14.1+
The query was compared to classified sequences in "Combined.rep.library.minlen50.usearch.fasta.classified.filtered"
FamDB:
```


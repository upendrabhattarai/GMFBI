We used a combination of ab initio and evidence-based annotation of the Thermobia genome assembly.
Haplotype 1 was used as the assembly because of its better quality compared to other haplotype or primary assembly.
All the protein data from NCBI and uniprot for Zygentoma were downloaded on 2023 Dec 02.

## 1. RNA-seq data preprocessing
Trimmomatic used for QC
```
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 HEADCROP:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35
```
## 2. De novo transcriptome assembly
Ran rcorrector and filtered-out uncorrected reads. Then `Trinity` `V.2.15.1` was used to assemble.
```
TRINITY=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/trinityrnaseq_latest.sif

singularity exec --cleanenv ${TRINITY} Trinity --CPU ${SLURM_CPUS_PER_TASK} --max_memory 100G \
  --seqType fq \
  --SS_lib_type RF \
  --left comma_separated_list_of_R1 \
  --right comma_separated_list_of_R2 \
  --output Trinity_out
```

Busco result on Trinity output
```
BUSCO version is: 5.1.3 
# The lineage dataset is: arthropoda_odb10 (Creation date: 2020-09-10, number of genomes: 90, number of BUSCOs: 1013)
# BUSCO was run in mode: transcriptome

        ***** Results: *****

        C:99.1%[S:30.1%,D:69.0%],F:0.5%,M:0.4%,n:1013      
        1004    Complete BUSCOs (C)                        
        305     Complete and single-copy BUSCOs (S)        
        699     Complete and duplicated BUSCOs (D)         
        5       Fragmented BUSCOs (F)                      
        4       Missing BUSCOs (M)                         
        1013    Total BUSCO groups searched                

Dependencies and versions:
        hmmsearch: 3.2
        metaeuk: 9dee7a78db0f2a8d6aafe7dbf18ac06bb6e23bf0-MPI
```


## 3. Genome-guided RNA-seq assembly
We used `HISAT2` `V.2.2.1` to map the in-house generated RNA-seq data to the assembly and then used `Stringtie` `V.2.2.1` to generate the RNA-seq assembly.

```
export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/mambaforge/bin:$PATH
# Step 1
hisat2-build -p 20 ../genome/T.dom.Hap1.final.fa Index.H1

# Step 2
SAMPLES="PREP0367_UBhat16381A_A01v1_T_d_Emb_0_S1 PREP0367_UBhat16381A_A02v1_T_d_M_Go_1_S9 PREP0367_UBhat16381A_A03v1_T_d_F_Nc_2_S17 PREP0367_UBhat16381A_B01v1_T_d_Nym_0_S2 PREP0367_UBhat16381A_B02v1_T_d_M_Go_2_S10 PREP0367_UBhat16381A_B03v1_T_d_F_Nc_3_S18 PREP0367_UBhat16381A_C01v1_T_d_M_Ant_1_S3 PREP0367_UBhat16381A_C02v1_T_d_M_Go_3_S11 PREP0367_UBhat16381A_C03v1_T_d_F_Go_1_S19 PREP0367_UBhat16381A_D01v1_T_d_M_Ant_2_S4 PREP0367_UBhat16381A_D02v1_T_d_M_Wb_0_S12 PREP0367_UBhat16381A_D03v1_T_d_F_Go_2_S20 PREP0367_UBhat16381A_E01v1_T_d_M_Ant_3_S5 PREP0367_UBhat16381A_E02v1_T_d_F_Ant_1_S13 PREP0367_UBhat16381A_E03v1_T_d_F_Go_3_S21 PREP0367_UBhat16381A_F01v1_T_d_M_Nc_1_S6 PREP0367_UBhat16381A_F02v1_T_d_F_Ant_2_S14 PREP0367_UBhat16381A_F03v1_T_d_F_Wb_0_S22 PREP0367_UBhat16381A_G01v1_T_d_M_Nc_2_S7 PREP0367_UBhat16381A_G02v1_T_d_F_Ant_3_S15 PREP0367_UBhat16381A_H01v1_T_d_M_Nc_3_S8 PREP0367_UBhat16381A_H02v1_T_d_F_Nc_1_S16"

for SAMPLES in $SAMPLES
do
hisat2 -p 20 --dta -x Index.H1 \
       -1 ../../2.T.dom.RNA_seq/C.trimmomatic/trim/${SAMPLES}_L004_R1_001_trim.fastq.gz \
       -2  ../../2.T.dom.RNA_seq/C.trimmomatic/trim/${SAMPLES}_L004_R2_001_trim.fastq.gz \
       -S ${SAMPLES}.Hap1.sam
done

# Step 3 Convert to bam

for SAMPLES in $SAMPLES
do
samtools sort -@ 20 -o ${SAMPLES}_Hap1.bam ${SAMPLES}.Hap1.sam
done


# Step 4 index bam files using samtools

for SAMPLES in $SAMPLES
do
samtools index ${SAMPLES}_Hap1.bam ${SAMPLES}_Hap1.bam.bai
done

# Step 5 Running stringtie
export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/mambaforge/bin:$PATH

for SAMPLES in $SAMPLES
do
stringtie ../${SAMPLES}_Hap1.bam -l ${SAMPLES} -p 20 --rf -o assembly.rf/${SAMPLES}_Hap1.gtf
done

# Step 6 Merge all gtf files to one.
stringtie --merge --rf -p 20 -o stringtie_merged.gtf assembly.rf/*.gtf

# Step 7 Get the assembly fasta from stringtie.gtf file using Agat tool

singularity exec /n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/agat_1.0.0--pl5321hdfd78af_0.sif agat_sp_extract_sequences.pl --cdna --gff stringtie_merged.gtf --fasta ../../genome/T.dom.Hap1.final.fa -o Stringtie_assembly.fasta
```

## 4. Pasa
`PASA` `V.2.5.3` was used with Trinity de novo assembly
Sqlite `V.3.45.0` was used for database management
Just a tip: To create a new .sqlite database
`.open --new T.dom.pasa.sqlite`

```
PASAHOME="/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/PASApipeline"

singularity exec /n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/pasapipeline_latest.sif $PASAHOME/Launch_PASA_pipeline.pl \
-c alignAssembly.config -C -R --ALIGNER gmap -g ../genome/final.soft.masked.TdH1.fasta.masked -t ../../2.T.dom.RNA_seq/E.trinity/Trinity_out.Trinity.fasta \
--transcribed_is_aligned_orient --CPU 64
```

Busco analysis on Pasa returned assembly:
```
# BUSCO version is: 5.1.3 
# The lineage dataset is: arthropoda_odb10 (Creation date: 2024-01-08, number of genomes: 90, number of BUSCOs: 1013)
# Summarized benchmarking in BUSCO notation for file /n/holyscratch01/extavour_lab/Everyone/upendrabhattarai/3.Tdom.H1.annota/pasa/T.dom.pasa.sqlite.assemblies.fasta
# BUSCO was run in mode: genome
# Gene predictor used: metaeuk

        ***** Results: *****

        C:86.7%[S:51.9%,D:34.8%],F:2.9%,M:10.4%,n:1013     
        879     Complete BUSCOs (C)                        
        526     Complete and single-copy BUSCOs (S)        
        353     Complete and duplicated BUSCOs (D)         
        29      Fragmented BUSCOs (F)                      
        105     Missing BUSCOs (M)                         
        1013    Total BUSCO groups searched                

Dependencies and versions:
        hmmsearch: 3.2
        metaeuk: 9dee7a78db0f2a8d6aafe7dbf18ac06bb6e23bf0-MPI
```

## 5. Transdecoder

`Transdecoder` `V.5.7.1` was used. We followed the script here to run transdecoder after Pasa
https://github.com/TransDecoder/TransDecoder/blob/master/sample_data/pasa_example/runMe.sh

```
singularity exec --cleanenv /n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/transdecoder_latest.sif \
        TransDecoder.LongOrfs -t ../T.dom.pasa.sqlite.assemblies.fasta \
        --gene_trans_map pasa.gene_trnas_map.txt -O pasa.transdecoder_workdir

```
## 6. Braker3
`Braker` `V.3.0.6` was used with the protein and Hisat2 aligned RNA-seq reads
```
BRAKER_PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/braker3.sif
PATH_TO_BAM=/n/holyscratch01/extavour_lab/Everyone/upendrabhattarai/3.Tdom.H1.annota/rna.seq.assembly.genome.guided

singularity exec --cleanenv $BRAKER_PATH braker.pl \
        --genome=../genome.2.masked/T.dom.Hap1.final.fa.masked \
        --prot_seq=../protein/merged.ncbi.uniprot.zygentoma.proteins.2023_12_02.fasta \
        --bam=$PATH_TO_BAM/PREP0367_UBhat16381A_A01v1_T_d_Emb_0_S1_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_A02v1_T_d_M_Go_1_S9_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_A03v1_T_d_F_Nc_2_S17_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_B01v1_T_d_Nym_0_S2_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_B02v1_T_d_M_Go_2_S10_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_B03v1_T_d_F_Nc_3_S18_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_C01v1_T_d_M_Ant_1_S3_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_C02v1_T_d_M_Go_3_S11_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_C03v1_T_d_F_Go_1_S19_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_D01v1_T_d_M_Ant_2_S4_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_D02v1_T_d_M_Wb_0_S12_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_D03v1_T_d_F_Go_2_S20_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_E01v1_T_d_M_Ant_3_S5_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_E02v1_T_d_F_Ant_1_S13_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_E03v1_T_d_F_Go_3_S21_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_F01v1_T_d_M_Nc_1_S6_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_F02v1_T_d_F_Ant_2_S14_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_F03v1_T_d_F_Wb_0_S22_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_G01v1_T_d_M_Nc_2_S7_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_G02v1_T_d_F_Ant_3_S15_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_H01v1_T_d_M_Nc_3_S8_Hap1.bam,$PATH_TO_BAM/PREP0367_UBhat16381A_H02v1_T_d_F_Nc_1_S16_Hap1.bam \
        --species=Braker.firebrat --softmasking --threads=36
```

Statistics of braker output
```
Type (3rd column)       Number  Size total (kb) Size mean (bp)  % of the genome /!\Results are rounding to two decimal places 
cds     152601  33392.40        218.82  0.59
exon    152601  33392.40        218.82  0.59
gene    19128   1179917.01      61685.33        20.88
intron  129685  1503949.52      11596.94        26.61
mrna    1372    32584.61        23749.71        0.58
start_codon     22913   68.70   3.00    0.00
stop_codon      22914   68.72   3.00    0.00
transcript      22916   1537341.93      67085.96        27.20
Total   524130  4320715.29      8243.59 76.46
```
BUSCO on the protein annotation from brake
```
# BUSCO version is: 5.1.3 
# The lineage dataset is: arthropoda_odb10 (Creation date: 2024-01-08, number of genomes: 90, number of BUSCOs: 1013)
# Summarized benchmarking in BUSCO notation for file /n/holyscratch01/extavour_lab/Everyone/upendrabhattarai/3.Tdom.H1.annota/Braker.2/braker/braker.aa
# BUSCO was run in mode: proteins

        ***** Results: *****

        C:94.9%[S:76.1%,D:18.8%],F:0.9%,M:4.2%,n:1013      
        961     Complete BUSCOs (C)                        
        771     Complete and single-copy BUSCOs (S)        
        190     Complete and duplicated BUSCOs (D)         
        9       Fragmented BUSCOs (F)                      
        43      Missing BUSCOs (M)                         
        1013    Total BUSCO groups searched                

Dependencies and versions:
        hmmsearch: 3.2
```
## 7. SNAP
Snap was trained with Braker output following the tutorial [here](https://github.com/KorfLab/SNAP)
`ZOE library version 2017-03-01` Software downloaded Jan 25 2024 from Github.
```
gff3_to_zff.pl genome.fasta braker.out.gff3 > snap.zff

fathom -validate snap.zff genome.fasta > snap.validate

fathom -categorize 1000 snap.zff genome.fasta

fathom -export 1000 -plus uni.*

fathom -validate export.ann export.dna

# make a subfolder: param, cd inside param and run the script below because it will produce many files
mkdir param
cd param
forge ../export.ann ../export.dna

# Finally run hmm-assembler.pl to glue the various models together to form an hmm parameter file.

hmm-assembler.pl snap.braker . > snap.braker.hmm
```

## 8. Miniprot
Aligned the downloaded protein set against the genome with `Miniprot` `V.0.12-r237`
```
export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/miniprot:$PATH

miniprot -Iut16 \
        --gff ../genome/T.dom.Hap1.final.fa \
        ../maker/proteins/merged.ncbi.uniprot.zygentoma.proteins.2023_12_02.fasta > prot.align.gff
```

## 9. EVidenceModeler
`EVidenceModeler` `V.2.1.0` was used.
Output from braker, snap, and transdecoder were converted to EVM format and merged together.
output from miniprot alignment was converted to EVM format.
output from pasa and stringtie were converted to EVM format and merged together.

Weight file for EVM was created.
`weight.txt`
```
ABINITIO_PREDICTION     SNAP    1
ABINITIO_PREDICTION     AUGUSTUS        1
ABINITIO_PREDICTION     GeneMark.hmm3   1
ABINITIO_PREDICTION     gmst    1
OTHER_PREDICTION        transdecoder    5
PROTEIN miniprot        2
TRANSCRIPT      exon    10
```

EVM was run with the following commands:
```
EVM_PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/evidencemodeler_latest.sif
PATH_TO_BAM=/n/holyscratch01/extavour_lab/Everyone/upendrabhattarai/3.Tdom.H1.annota/rna.seq.assembly.genome.guided

singularity exec --cleanenv $EVM_PATH EVidenceModeler \
        --sample_id T.dom_H1 \
        --genome ../genome.2.masked/T.dom.Hap1.final.fa.masked \
        --weights ./weight.txt \
        --gene_predictions input/braker.snap.transdecoder.merged.2EVMformat.gff3 \
        --protein_alignments input/miniprot_prot_align_2EVMformat.gff3 \
        --transcript_alignments input/pasa_stringtie_merged_2EVM.gff3 \
        --segmentSize 7000000 \
        --overlapSize 700000 \
        --CPU 36
```

EVM output statistics.
`agat_sq_stat_basic`
```
Type (3rd column)       Number  Size total (kb) Size mean (bp)  % of the genome /!\Results are rounding to two decimal places 
cds     135149  31144.91        230.45  0.55
exon    135149  31144.91        230.45  0.55
gene    25914   1497045.14      57769.74        26.49
mrna    25914   1497045.14      57769.74        26.49
Total   322126  3056380.10      9488.15 54.09
```
`agat_sp_statistics`
```
Compute mrna with isoforms if any

Number of gene                               25914
Number of mrna                               25914
Number of cds                                25914
Number of exon                               135149
Number of exon in cds                        135149
Number of intron in cds                      109235
Number of intron in exon                     109235
Number gene overlapping                      0
Number of single exon gene                   8712
Number of single exon mrna                   8712
mean mrnas per gene                          1.0
mean cdss per mrna                           1.0
mean exons per mrna                          5.2
mean exons per cds                           5.2
mean introns in cdss per mrna                4.2
mean introns in exons per mrna               4.2
Total gene length                            1497045143
Total mrna length                            1497045143
Total cds length                             31144909
Total exon length                            31144909
Total intron length per cds                  1465900234
Total intron length per exon                 1465900234
mean gene length                             57769
mean mrna length                             57769
mean cds length                              1201
mean exon length                             230
mean cds piece length                        230
mean intron in cds length                    13419
mean intron in exon length                   13419
% of genome covered by gene                  26.5
% of genome covered by mrna                  26.5
% of genome covered by cds                   0.6
% of genome covered by exon                  0.6
% of genome covered by intron from cds       25.9
% of genome covered by intron from exon      25.9
Longest gene                                 2312799
Longest mrna                                 2312799
Longest cds                                  49521
Longest exon                                 36630
Longest cds piece                            36630
Longest intron into cds part                 498132
Longest intron into exon part                498132
Shortest gene                                150
Shortest mrna                                150
Shortest cds                                 150
Shortest exon                                2
Shortest cds piece                           2
Shortest intron into cds part                21
Shortest intron into exon part               21

Re-compute mrna without isoforms asked. We remove shortest isoforms if any
Number of gene                               25914
Number of mrna                               25914
Number of cds                                25914
Number of exon                               135149
Number of exon in cds                        135149
Number of intron in cds                      109235
Number of intron in exon                     109235
Number gene overlapping                      0
Number of single exon gene                   8712
Number of single exon mrna                   8712
mean mrnas per gene                          1.0
mean cdss per mrna                           1.0
mean exons per mrna                          5.2
mean exons per cds                           5.2
mean introns in cdss per mrna                4.2
mean introns in exons per mrna               4.2
Total gene length                            1497045143
Total mrna length                            1497045143
Total cds length                             31144909
Total exon length                            31144909
Total intron length per cds                  1465900234
Total intron length per exon                 1465900234
mean gene length                             57769
mean mrna length                             57769
mean cds length                              1201
mean exon length                             230
mean cds piece length                        230
mean intron in cds length                    13419
mean intron in exon length                   13419
% of genome covered by gene                  26.5
% of genome covered by mrna                  26.5
% of genome covered by cds                   0.6
% of genome covered by exon                  0.6
% of genome covered by intron from cds       25.9
% of genome covered by intron from exon      25.9
Longest gene                                 2312799
Longest mrna                                 2312799
Longest cds                                  49521
Longest exon                                 36630
Longest cds piece                            36630
Longest intron into cds part                 498132
Longest intron into exon part                498132
Shortest gene                                150
Shortest mrna                                150
Shortest cds                                 150
Shortest exon                                2
Shortest cds piece                           2
Shortest intron into cds part                21
Shortest intron into exon part               21
```

BUSCO assessment of coded proteins from EVM (arthropoda database)
```
# BUSCO version is: 5.1.3 
# The lineage dataset is: arthropoda_odb10 (Creation date: 2024-01-08, number of genomes: 90, number of BUSCOs: 1013)
# Summarized benchmarking in BUSCO notation for file /n/holyscratch01/extavour_lab/Everyone/upendrabhattarai/3.Tdom.H1.annota/EVM.4/T.dom_H1.EVM.pep
# BUSCO was run in mode: proteins

        ***** Results: *****

        C:97.0%[S:93.3%,D:3.7%],F:0.8%,M:2.2%,n:1013       
        982     Complete BUSCOs (C)                        
        945     Complete and single-copy BUSCOs (S)        
        37      Complete and duplicated BUSCOs (D)         
        8       Fragmented BUSCOs (F)                      
        23      Missing BUSCOs (M)                         
        1013    Total BUSCO groups searched                

Dependencies and versions:
        hmmsearch: 3.2
```
BUSCO assessment of coded proteins from EVM (insecta database)
```
# BUSCO version is: 5.1.3 
# The lineage dataset is: insecta_odb10 (Creation date: 2024-01-08, number of genomes: 75, number of BUSCOs: 1367)
# Summarized benchmarking in BUSCO notation for file /n/holyscratch01/extavour_lab/Everyone/upendrabhattarai/3.Tdom.H1.annota/EVM.4/T.dom_H1.EVM.pep
# BUSCO was run in mode: proteins

        ***** Results: *****

        C:96.2%[S:92.0%,D:4.2%],F:0.4%,M:3.4%,n:1367       
        1316    Complete BUSCOs (C)                        
        1258    Complete and single-copy BUSCOs (S)        
        58      Complete and duplicated BUSCOs (D)         
        5       Fragmented BUSCOs (F)                      
        46      Missing BUSCOs (M)                         
        1367    Total BUSCO groups searched                

Dependencies and versions:
        hmmsearch: 3.2
```

## 10. Functional annotation with Blast
Ran `Blastp` `V.2.25.0` with uniprot_sprot database
uniprot_sprot.fasta downloaded from https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
updated dataset on 2024-01-24 10:00

```
singularity exec --cleanenv /n/holylabs/LABS/extavour_lab/Everyone/bin/blast_2.15.0.sif \
blastp -query ../EVM.4/T.dom_H1.EVM.pep -db uniprot.database/uniprot_sprot.fasta -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -num_threads 36 -out output.blastp
```

## 11. Functional annotation with interproscan
Ran `Interproscan` `V.5.66.98.0`

```
module load jdk/20.0.1-fasrc01
export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/interproscan-5.66-98.0:$PATH

interproscan.sh -appl pfam -dp -f TSV -goterms -iprlookup -pa -t p --cpu 36 -i T.dom_H1.EVM.pep -o EVM.4.iprscan
```

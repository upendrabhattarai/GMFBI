# Genome assembly pipeline for Thermobia genome
Pacbio HiFi sequencing data were obtained from Bauer core in bam format.
## Removing adapters
Used `HiFiAdapterFilt` `V.2.0.1`

```
module load bamtools/2.3.0-fasrc01
module load GCC/8.2.0-2.31.1 OpenMPI/3.1.3 BLAST+/2.9.0
module load GCCcore/8.2.0 pigz/2.4

export PATH=$PATH:/n/home12/upendrabhattarai/bin/HiFiAdapterFilt
export PATH=$PATH:/n/home12/upendrabhattarai/bin/HiFiAdapterFilt/DB

bash pbadapterfilt.sh -p ../1_A01/m64408e_230104_191438 -t 6 -o 1_A_filt
```
Output
`m64408e_230104_191438.hifi_reads.stats`
```
Started on Mon May 22 19:45:26 EDT 2023
For the m64408e_230104_191438.hifi_reads dataset:
Removing reads containing adapters a minimum of 44 bp in length and 97% match.

Number of ccs reads: 2715798
Number of adapter contaminated ccs reads: 7 (0.000257751% of total)
Number of ccs reads retained: 2715791 (99.9997% of total)

Finished on Mon May 22 20:46:10 EDT 2023
```
`m64408e_230106_060914.hifi_reads.stats`
```
Started on Mon May 22 20:46:14 EDT 2023
For the m64408e_230106_060914.hifi_reads dataset:
Removing reads containing adapters a minimum of 44 bp in length and 97% match.

Number of ccs reads: 1897718
Number of adapter contaminated ccs reads: 2 (0.00010539% of total)
Number of ccs reads retained: 1897716 (99.9999% of total)

Finished on Mon May 22 21:28:52 EDT 2023
```
`m64408e_230107_170457.hifi_reads.stats`
```
Started on Mon May 22 21:28:58 EDT 2023
For the m64408e_230107_170457.hifi_reads dataset:
Removing reads containing adapters a minimum of 44 bp in length and 97% match.

Number of ccs reads: 2800335
Number of adapter contaminated ccs reads: 8 (0.00028568% of total)
Number of ccs reads retained: 2800327 (99.9997% of total)

Finished on Mon May 22 22:32:53 EDT 2023
```
`m64408e_230109_021420.hifi_reads.stats`
```
Started on Mon May 22 22:32:58 EDT 2023
For the m64408e_230109_021420.hifi_reads dataset:
Removing reads containing adapters a minimum of 44 bp in length and 97% match.

Number of ccs reads: 2492527
Number of adapter contaminated ccs reads: 2 (8.02399e-05% of total)
Number of ccs reads retained: 2492525 (99.9999% of total)

Finished on Mon May 22 23:30:02 EDT 2023
```
The Multiqc report after adapter filtering is [here](multiqc_report_hifiadapterfilt.html) (Need to download the html file to visualize it properly)

# Removing contaminant reads
`Kraken2` `V2.1.3` was used with the standard Kraken database to filter out contaminants. The extract_kraken_reads.py script is a part of Kraken tools. you can also find it (here)[extract_kraken_reads.py]
```
export PATH="/n/holyscratch01/extavour_lab/Everyone/upendrabhattarai/bin/kraken2/kraken2:$PATH"
export KRAKEN2_DB_PATH=/n/holyscratch01/extavour_lab/Everyone/upendrabhattarai/bin/kraken2/kraken-std-database

## Step 1 Reads classification
kraken2 --db kraken2-std --threads 20 --classified-out classified.out --unclassified-out unclassified.out --report reports.txt --use-names \
m64408e_230106_060914.hifi_reads.filt.fasta \
m64408e_230109_021420.hifi_reads.filt.fasta \
m64408e_230104_191438.hifi_reads.filt.fasta \
m64408e_230107_170457.hifi_reads.filt.fasta
```
output:
```
Loading database information... done.
9906367 sequences (106095.59 Mbp) processed in 459.887s (1292.5 Kseq/m, 13841.94 Mbp/m).
  9395789 sequences classified (94.85%)
  510578 sequences unclassified (5.15%)
```

Step 2 Extract the clean reads with no bacterial contaminants.
```
export PATH="/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/mambaforge/bin/:$PATH"
python ./extract_kraken_reads.py \
       -k kraken2_61553288.out \
       -s m64408e_230107_170457.hifi_reads.filt.fasta \
       --output m64408e_230107_170457.hifi_reads.filt.kraken.nobac.fasta \
       --taxid 2 \
       --exclude \
       --include-children \
       --report kraken.filter.report
```

Step 3 Extract only the bacterial reads
```
python ./extract_kraken_reads.py \
        -k kraken2_61553288.out \
        -s m64408e_230106_060914.hifi_reads.filt.fasta \
        --output m64408e_230106_060914.hifi_onlybac.fasta \
        --taxid 2 \
        --include-children \
        --report only.bac.report
```
## Hi-C QC
The output from step 2 is clean Pacbio reads to be used in assembly with Hifiasm. We also have Hi-C reads generated with the Arima Hi-C hi-coverage pipeline.
To process these reads, we will carryout QC with Trimmomatic.

`Trimmomatic` `V.0.39` was used.

```
export PATH=/n/holylabs/LABS/extavour_lab/Lab/softwares/Trimmomatic-0.39:$PATH
module load jdk
mkdir trim
mkdir unpaired

java -jar /n/holylabs/LABS/extavour_lab/Lab/softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 16 \
SUB14149_A1_TdomHiC_S23_L004_R1_001.fastq.gz SUB14149_A1_TdomHiC_S23_L004_R2_001.fastq.gz \
"trim/SUB14149_A1_TdomHiC_S23_L004_R1_001.fastq.gz" "unpaired/SUB14149_A1_TdomHiC_S23_L004_R1_001.fastq.gz" \
"trim/SUB14149_A1_TdomHiC_S23_L004_R2_001.fastq.gz" "unpaired/SUB14149_A1_TdomHiC_S23_L004_R2_001.fastq.gz" \
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 HEADCROP:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35
```
The adapter file provided was [TruSeq3-PE-2.fa](TruSeq3-PE-2.fa)
Output:
```
Input Read Pairs: 840098291 Both Surviving: 752751661 (89.60%) Forward Only Surviving: 64729274 (7.70%) Reverse Only Surviving: 11043465 (1.31%) Dropped: 11573891 (1.38%)
```
Find the [MultiQC report here](multiqc_report_hic_trimmomatic.html)

## Assembly
We used the output from step 2 of the `Kraken` pipeline for HiFi and QC'ed Hi-C data to assemble the genome using Hifiasm.
`Hifiasm` `V0.19.5-r587` was used for the assembly using both HiFi and Hi-C reads.

```
# Step 1 assembly
export PATH="/n/home12/upendrabhattarai/bin/hifiasm.new/hifiasm-0.19.5:$PATH"
hifiasm -o Tdom_HiFC_s0.2_asm -s 0.1 -l 3 -t36 --n-weight 4 --n-perturb 85000 --f-perturb 0.10 \
--h1 SUB14149_A1_TdomHiC_S23_L004_R1_001.trim.fastq.gz \
--h2 SUB14149_A1_TdomHiC_S23_L004_R2_001.trim.fastq.gz \
*reads.filt.kraken.nobac.fasta

# Step 2 To produce fasta file from gfa use the following code

awk '/^S/{print ">"$2;print $3}' TdomHiFiC.s_0.2.asm.hic.hap1.p_ctg.gfa > TdomHiFiC.s_0.2.asm.hic.hap1.p_ctg.fasta

awk '/^S/{print ">"$2;print $3}' TdomHiFiC.s_0.2.asm.hic.hap2.p_ctg.gfa > TdomHiFiC.s_0.2.asm.hic.hap2.p_ctg.fasta

awk '/^S/{print ">"$2;print $3}' Tdom_HiFC_s0.2_asm.hic.p_ctg.gfa > Tdom_HiFC_s0.2_asm.hic.p_ctg.fasta
```
<details>
<summary markdown='span'><b>hifiasm.out ▶️ </b></summary>

Version: 0.19.5-r587

CMD: hifiasm -o Tdom_HiFC_s0.2_asm -s 0.1 -l 3 -t36 --n-weight 4 --n-perturb 85000 --f-perturb 0.10 --h1 SUB14149_A1_TdomHiC_S23_L004_R1_001.trim.fastq.gz --h2 SUB14149_A1_TdomHiC_S23_L004_R2_001.trim.fastq.gz m64408e_230104_191438.hifi_reads.filt.kraken.nobac.fasta m64408e_230106_060914.hifi_reads.filt.kraken.nobac.fasta m64408e_230107_170457.hifi_reads.filt.kraken.nobac.fasta m64408e_230109_021420.hifi_reads.filt.kraken.nobac.fasta

Real time: 53589.960 sec; CPU: 581079.470 sec; Peak RSS: 265.478 GB
</details>

## Hi-C mapping
We mapped Hi-C reads back to the assembly.
We have primary assembly, haplotype-1, and haplotype-2 from hifiasm output.
Arima pipeline was used to map the hic data back to the assembly. The pipeline is as follows:

```
#for softwares installed with conda
export PATH=~/bin/miniconda3/bin:$PATH

SRA='SUB14149_A1_TdomHiC_S23_L004'
LABEL='Tdom.asm.H1'
IN_DIR='./Hi-C.trim'
REF='./Tdom_HiFC_s0.1_nw4_np85k_fp0.1_asm.hic.hap1.p_ctg.fasta'
FAIDX='$REF.fai'
PREFIX='bwa_index'
RAW_DIR='./bam.output.raw'
FILT_DIR='./bam.output.filtered'
FILTER='/n/holyscratch01/extavour_lab/Everyone/upendrabhattarai/pacbio_Thermobia.genome/230104_PREP0300_UBhattarai/r64408e_20230104_134539/11.hicmapping/hap2/scripts/filter_five_end.pl'
COMBINER='/n/holyscratch01/extavour_lab/Everyone/upendrabhattarai/pacbio_Thermobia.genome/230104_PREP0300_UBhattarai/r64408e_20230104_134539/11.hicmapping/hap2/scripts/two_read_bam_combiner.pl'
STATS='/n/holyscratch01/extavour_lab/Everyone/upendrabhattarai/pacbio_Thermobia.genome/230104_PREP0300_UBhattarai/r64408e_20230104_134539/11.hicmapping/hap2/scripts/get_stats.pl'
PICARD='~/bin/picard.jar'
TMP_DIR='./temp'
PAIR_DIR='./bams.paired'
REP_DIR='./deduplicated'
REP_LABEL=$LABEL\_rep1
MERGE_DIR='./merged.alignment'
MAPQ_FILTER=10
CPU=36

#echo "### Step 1: Index reference" # Run only once! Skip this step if you have already generated BWA index files
bwa index -a bwtsw -p $PREFIX $REF

#echo "### Step 2.A: FASTQ to BAM (1st)"
mkdir bam.output.raw
bwa mem -t $CPU bwa_index Hi-C.trim/SUB14149_A1_TdomHiC_S23_L004_R1_001.trim.fastq.gz | samtools view -@ $CPU -Sb - > bam.output.raw/SUB14149_A1_TdomHiC_S23_L004_1.bam

#echo "### Step 2.B: FASTQ to BAM (2nd)"
bwa mem -t $CPU bwa_index Hi-C.trim/SUB14149_A1_TdomHiC_S23_L004_R2_001.trim.fastq.gz | samtools view -@ $CPU -Sb - > bam.output.raw/SUB14149_A1_TdomHiC_S23_L004_2.bam

#echo "### Step 3.A: Filter 5' end (1st)"
mkdir bam.output.filtered
samtools view -@ $CPU -h bam.output.raw/SUB14149_A1_TdomHiC_S23_L004_1.bam | perl $FILTER | samtools view -@ $CPU -Sb - > bam.output.filtered/SUB14149_A1_TdomHiC_S23_L004_1.bam

#echo "### Step 3.B: Filter 5' end (2nd)"
samtools view -@ $CPU -h bam.output.raw/SUB14149_A1_TdomHiC_S23_L004_2.bam | perl $FILTER | samtools view -@ $CPU -Sb - > bam.output.filtered/SUB14149_A1_TdomHiC_S23_L004_2.bam

#echo "### Step 4A: Pair reads & mapping quality filter"
mkdir temp
perl $COMBINER bam.output.filtered/SUB14149_A1_TdomHiC_S23_L004_1.bam bam.output.filtered/SUB14149_A1_TdomHiC_S23_L004_2.bam samtools $MAPQ_FILTER | samtools view -@ $CPU -bS -t $FAIDX - | samtools sort -@ $CPU -o temp/$SRA.bam -

#echo "### Step 4.B: Add read group"
module load jdk
mkdir bams.paired
java -Xmx4G -Djava.io.tmpdir=temp/ -jar /n/home12/upendrabhattarai/bin/picard.jar AddOrReplaceReadGroups INPUT=temp/$SRA.bam OUTPUT=bams.paired/$SRA.bam ID=$SRA LB=$SRA SM=$LABEL PL=ILLUMINA PU=none


#echo "### Step 5: Mark duplicates"
module load jdk
mkdir deduplicated
java -Xmx100G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ \
       -jar /n/home12/upendrabhattarai/bin/picard.jar MarkDuplicates \
       INPUT=$PAIR_DIR/$SRA.bam \
       OUTPUT=$REP_DIR/$REP_LABEL.bam \
       METRICS_FILE=$REP_DIR/metrics.$REP_LABEL.txt \
       TMP_DIR=$TMP_DIR \
       ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

#echo "### Step 6: Indexing"
samtools index $REP_DIR/$REP_LABEL.bam -@ 20

#echo "### step 7: stats"
perl $STATS $REP_DIR/$REP_LABEL.bam > $REP_DIR/$REP_LABEL.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"

```
The final output file will be in the `deduplicated` folder.
Mapping statistics for 
Hap-1
```
All     174640407
All intra       83663096
All intra 1kb   38694447
All intra 10kb  30965119
All intra 15kb  29058249
All intra 20kb  27739395
All inter       90977311
```
Hap-2
```
All     136652400
All intra       68640181
All intra 1kb   30701886
All intra 10kb  24208907
All intra 15kb  22642139
All intra 20kb  21581273
All inter       68012219
```
Primary assembly
```
All     161695594
All intra       82212781
All intra 1kb   39663910
All intra 10kb  32496829
All intra 15kb  30733566
All intra 20kb  29517735
All inter       79482813
```
## Hi-C scaffolding
We scaffolded all three assemblies using the `Yahs` `V1.2` pipeline.
For Hap-1
```
#for samtools to index
export PATH=~/bin/miniconda3/bin:$PATH
#samtools faidx T.domestica.M.filtered.fasta

#for softwares installed with conda
export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/yahs/:$PATH
#yahs ../Tdom_HiFC_s0.1_nw4_np85k_fp0.1_asm.hic.hap1.p_ctg.fasta \
#../deduplicated/Tdom.asm.H1_rep1.bam

### creating a .hic contact map for visualization
OUT="yahs.out" #name of output file from yahs
CONTIGS="../Tdom_HiFC_s0.1_nw4_np85k_fp0.1_asm.hic.hap1.p_ctg.fasta"

#juicer pre -a -o out_Pri_JBAT ${OUT}.bin ${OUT}_scaffolds_final.agp ${CONTIGS}.fai > out_JBAT.log 2>&1

module load jdk
JUICER=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/juicer_tools.2.20.00.jar
java -Xmx80G -jar ${JUICER} pre out_Pri_JBAT.txt out_Pri_JBAT.hic <(echo "assembly 1412685812")
```
For Hap-2
```
#for samtools to index
#export PATH=~/bin/miniconda3/bin:$PATH
#samtools faidx T.domestica.M.filtered.fasta

#for softwares installed with conda
#export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/yahs/:$PATH
#yahs ../Tdom_HiFC_s0.1_nw4_np85k_fp0.1_asm.hic.hap2.p_ctg.fasta \
#../deduplicated/Tdom.asm.Pri_rep1.bam

### creating a .hic contact map for visualization
#OUT="yahs.out" #name of output file from yahs
#CONTIGS="../Tdom_HiFC_s0.1_nw4_np85k_fp0.1_asm.hic.hap2.p_ctg.fasta"

#juicer pre -a -o out_Pri_JBAT ${OUT}.bin ${OUT}_scaffolds_final.agp ${CONTIGS}.fai > out_JBAT.log 2>&1

module load jdk/20.0.1-fasrc01
JUICER=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/juicer_tools.2.20.00.jar
java -Xmx80G -jar ${JUICER} pre out_Pri_JBAT.txt out_Pri_JBAT.hic <(echo "assembly 1197285486")
```
For Primary assembly
```
#for samtools to index
export PATH=~/bin/miniconda3/bin:$PATH
#samtools faidx T.domestica.M.filtered.fasta

#for softwares installed with conda
export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/yahs/:$PATH
#yahs ../Tdom_HiFC_s0.1_nw4_np85k_fp0.1_asm.hic.p_ctg.fasta \
#../deduplicated/Tdom.asm.Pri_rep1.bam

### creating a .hic contact map for visualization
OUT="yahs.out" #name of output file from yahs
CONTIGS="../Tdom_HiFC_s0.1_nw4_np85k_fp0.1_asm.hic.p_ctg.fasta"

#juicer pre -a -o out_Pri_JBAT ${OUT}.bin ${OUT}_scaffolds_final.agp ${CONTIGS}.fai > out_JBAT.log 2>&1

module load jdk/20.0.1-fasrc01

JUICER=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/juicer_tools.2.20.00.jar
java -Xmx100G -jar ${JUICER} pre out_Pri_JBAT.txt out_Pri_JBAT.hic <(echo "assembly 1493420936")
```
Statistics after scaffolding with Yahs 
for Hap-1, there was no need for manual editing of the Hi-C contact map. All the signals looked very good and clear.
```
Main genome scaffold total: 2531
Main genome contig total:   5129
Main genome scaffold sequence total: 5651.0 MB
Main genome contig sequence total:   5650.7 MB (->  0.0% gap)
Main genome scaffold N/L50: 9/288.8 MB
Main genome contig N/L50:   397/3.5 MB
Number of scaffolds > 150000 KB: 18
% main genome in scaffolds > 150000 KB: 87.3%

 Minimum    Number    Number     Total        Total     Scaffold
Scaffold      of        of      Scaffold      Contig     Contig
 Length   Scaffolds  Contigs     Length       Length    Coverage
--------  ---------  -------  -----------  -----------  --------
    All     2,531      5,129  5,651,003,049  5,650,743,249   100.00%
   1 kb     2,531      5,129  5,651,003,049  5,650,743,249   100.00%
 2.5 kb     2,520      5,118  5,650,985,049  5,650,725,249   100.00%
   5 kb     2,515      5,113  5,650,968,049  5,650,708,249   100.00%
  10 kb     2,500      5,098  5,650,845,506  5,650,585,706   100.00%
  25 kb     2,081      4,679  5,642,883,583  5,642,623,783   100.00%
  50 kb     1,583      4,180  5,625,723,319  5,625,463,619   100.00%
 100 kb     1,060      3,647  5,588,754,951  5,588,496,251   100.00%
 250 kb       620      3,195  5,517,352,936  5,517,095,436   100.00%
 500 kb       380      2,938  5,433,188,187  5,432,932,387   100.00%
   1 mb       178      2,724  5,288,333,668  5,288,079,068   100.00%
 2.5 mb        47      2,570  5,088,252,866  5,088,000,566   100.00%
   5 mb        24      2,538  5,011,456,855  5,011,205,455    99.99%
# BUSCO version is: 5.1.3 
# The lineage dataset is: arthropoda_odb10 (Creation date: 2020-09-10, number of genomes: 90, number of BUSCOs: 1013)
# BUSCO was run in mode: genome
# Gene predictor used: metaeuk
***** Results: *****

        C:99.5%[S:94.7%,D:4.8%],F:0.2%,M:0.3%,n:1013       
        1008    Complete BUSCOs (C)                        
        959     Complete and single-copy BUSCOs (S)        
        49      Complete and duplicated BUSCOs (D)         
        2       Fragmented BUSCOs (F)                      
        3       Missing BUSCOs (M)                         
        1013    Total BUSCO groups searched                

Dependencies and versions:
        hmmsearch: 3.2
        metaeuk: 9dee7a78db0f2a8d6aafe7dbf18ac06bb6e23bf0-MPI

```
For the Hap-2 and Primary assembly, we made a few manual edits to the Hi-C contact map in `Juicebox` `V.2.15`. It produced a *review.assembly file and post-processed to produce the final assembly with the juicer in the command line again.

Juicer post-processing after manual edits:
For Hap-2
```
export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/yahs/:$PATH
#juicer post -o Tdom.final.Hap2.asm \
#       out_hap2_JBAT.review.assembly \
#       ../out_hap2_JBAT.liftover.agp \
#       ../../Tdom_HiFC_s0.1_nw4_np85k_fp0.1_asm.hic.hap2.p_ctg.fasta

## for samtools to index the reviewed .fa
export PATH=~/bin/miniconda3/bin:$PATH
samtools faidx Tdom.final.Hap2.asm.FINAL.fa

### creating a .hic contact map for visualization

juicer pre ../yahs.out.bin \
        Tdom.final.Hap2.asm.FINAL.agp \
        ../../Tdom_HiFC_s0.1_nw4_np85k_fp0.1_asm.hic.hap2.p_ctg.fasta.fai \
        | sort -k2,2d -k6,6d -T ./ --parallel=36 -S80G \
        | awk 'NF' > alignments_sorted.txt.part && mv alignments_sorted.txt.part alignments_sorted.txt


#cat juice_post_10834303.err |grep "PRE_C_SIZE:" | cut -d' ' -f 2,3 > scaffolds_final.chrom.sizes


#module load jdk
#JUICER=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/juicer_tools.2.20.00.jar
#java -Xmx80G -jar ${JUICER} pre alignments_sorted.txt out.hic.part scaffolds_final.chrom.sizes && mv out.hic.part out.hic
```

For Primary assembly
```
export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/yahs/:$PATH
#juicer post -o Tdom.final.Pri.asm \
#       out_Pri_JBAT.review.assembly \
#       ../out_Pri_JBAT.liftover.agp \
#       ../../Tdom_HiFC_s0.1_nw4_np85k_fp0.1_asm.hic.p_ctg.fasta

## for samtools to index the reviewed .fa
export PATH=~/bin/miniconda3/bin:$PATH
#samtools faidx Tdom.final.Pri.asm.FINAL.fa

### creating a .hic contact map for visualization

#juicer pre ../yahs.out.bin Tdom.final.Pri.asm.FINAL.agp \
#       ../../Tdom_HiFC_s0.1_nw4_np85k_fp0.1_asm.hic.p_ctg.fasta.fai \
#       | sort -k2,2d -k6,6d -T ./ --parallel=36 -S80G \
#       | awk 'NF' > alignments_sorted.txt.part && 
#       mv alignments_sorted.txt.part alignments_sorted.txt


#cat juice_post_10834303.err |grep "PRE_C_SIZE:" | cut -d' ' -f 2,3 > scaffolds_final.chrom.sizes


module load jdk
JUICER=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/juicer_tools.2.20.00.jar
java -Xmx80G -jar ${JUICER} pre alignments_sorted.txt out.hic.part scaffolds_final.chrom.sizes && mv out.hic.part out.hic
```

Assembly statistics for Hap-2 final
```
Main genome scaffold total: 1187
Main genome contig total:   3647
Main genome scaffold sequence total: 4789.4 MB
Main genome contig sequence total:   4789.1 MB (->  0.0% gap)
Main genome scaffold N/L50: 8/278.5 MB
Main genome contig N/L50:   339/3.4 MB
Number of scaffolds > 150000 KB: 16
% main genome in scaffolds > 150000 KB: 87.5%

 Minimum    Number    Number     Total        Total     Scaffold
Scaffold      of        of      Scaffold      Contig     Contig
 Length   Scaffolds  Contigs     Length       Length    Coverage
--------  ---------  -------  -----------  -----------  --------
    All     1,187      3,647  4,789,387,944  4,789,141,944    99.99%
   1 kb     1,187      3,647  4,789,387,944  4,789,141,944    99.99%
 2.5 kb     1,171      3,631  4,789,366,944  4,789,120,944    99.99%
   5 kb     1,164      3,624  4,789,342,944  4,789,096,944    99.99%
  10 kb     1,160      3,620  4,789,313,438  4,789,067,438    99.99%
  25 kb     1,052      3,512  4,787,283,513  4,787,037,513    99.99%
  50 kb       792      3,252  4,778,169,251  4,777,923,251    99.99%
 100 kb       504      2,958  4,757,839,855  4,757,594,455    99.99%
 250 kb       264      2,712  4,719,475,814  4,719,231,014    99.99%
 500 kb       161      2,607  4,683,321,651  4,683,077,051    99.99%
   1 mb        80      2,516  4,625,307,012  4,625,063,412    99.99%
 2.5 mb        32      2,440  4,555,452,844  4,555,212,044    99.99%
   5 mb        22      2,405  4,521,540,488  4,521,302,188    99.99%

# BUSCO version is: 5.1.3 
# The lineage dataset is: arthropoda_odb10 (Creation date: 2020-09-10, number of genomes: 90, number of BUSCOs: 1013)
# BUSCO was run in mode: genome
# Gene predictor used: metaeuk

        ***** Results: *****

        C:96.0%[S:94.4%,D:1.6%],F:0.6%,M:3.4%,n:1013       
        972     Complete BUSCOs (C)                        
        956     Complete and single-copy BUSCOs (S)        
        16      Complete and duplicated BUSCOs (D)         
        6       Fragmented BUSCOs (F)                      
        35      Missing BUSCOs (M)                         
        1013    Total BUSCO groups searched                

Dependencies and versions:
        hmmsearch: 3.2
        metaeuk: 9dee7a78db0f2a8d6aafe7dbf18ac06bb6e23bf0-MPI
```
Statistics for final primary assembly
```
Main genome scaffold total: 2192
Main genome contig total:   3922
Main genome scaffold sequence total: 5973.9 MB
Main genome contig sequence total:   5973.7 MB (->  0.0% gap)
Main genome scaffold N/L50: 9/294.1 MB
Main genome contig N/L50:   281/5.4 MB
Number of scaffolds > 150000 KB: 18
% main genome in scaffolds > 150000 KB: 85.5%

 Minimum    Number    Number     Total        Total     Scaffold
Scaffold      of        of      Scaffold      Contig     Contig
 Length   Scaffolds  Contigs     Length       Length    Coverage
--------  ---------  -------  -----------  -----------  --------
    All     2,192      3,922  5,973,856,745  5,973,683,745   100.00%
   1 kb     2,192      3,922  5,973,856,745  5,973,683,745   100.00%
 2.5 kb     2,178      3,908  5,973,839,745  5,973,666,745   100.00%
   5 kb     2,175      3,905  5,973,828,745  5,973,655,745   100.00%
  10 kb     2,170      3,900  5,973,791,661  5,973,618,661   100.00%
  25 kb     1,770      3,500  5,966,427,671  5,966,254,671   100.00%
  50 kb     1,455      3,185  5,955,996,124  5,955,823,124   100.00%
 100 kb     1,105      2,828  5,931,181,573  5,931,009,273   100.00%
 250 kb       715      2,422  5,867,367,260  5,867,196,560   100.00%
 500 kb       449      2,142  5,772,446,206  5,772,276,906   100.00%
   1 mb       231      1,904  5,615,910,202  5,615,742,902   100.00%
 2.5 mb        70      1,716  5,367,014,557  5,366,849,957   100.00%
   5 mb        27      1,664  5,224,069,084  5,223,905,384   100.00%





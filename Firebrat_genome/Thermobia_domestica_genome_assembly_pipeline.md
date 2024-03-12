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

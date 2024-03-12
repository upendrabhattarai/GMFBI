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
`Kraken2` was used with the standard Kraken database to filter out contaminants.
```
export PATH="/n/holyscratch01/extavour_lab/Everyone/upendrabhattarai/bin/kraken2/kraken2:$PATH"
export KRAKEN2_DB_PATH=/n/holyscratch01/extavour_lab/Everyone/upendrabhattarai/bin/kraken2/kraken-std-database

## Step 1
kraken2 --db kraken2-std --threads 20 --classified-out classified.out --unclassified-out unclassified.out --report reports.txt --use-names \
m64408e_230106_060914.hifi_reads.filt.fasta \
m64408e_230109_021420.hifi_reads.filt.fasta \
m64408e_230104_191438.hifi_reads.filt.fasta \
m64408e_230107_170457.hifi_reads.filt.fasta

## Step 2
# Extract the clean reads with no bacterial contaminants. The extract_kraken_reads.py script is a part of kraken tools. you can also find it (here)[extract_kraken_reads.py]
export PATH="/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/mambaforge/bin/:$PATH"
python ./extract_kraken_reads.py \
       -k kraken2_61553288.out \
       -s m64408e_230107_170457.hifi_reads.filt.fasta \
       --output m64408e_230107_170457.hifi_reads.filt.kraken.nobac.fasta \
       --taxid 2 \
       --exclude \
       --include-children \
       --report kraken.filter.report

## Step 3
# Extract only the bacterial reads
python ./extract_kraken_reads.py \
        -k kraken2_61553288.out \
        -s m64408e_230106_060914.hifi_reads.filt.fasta \
        --output m64408e_230106_060914.hifi_onlybac.fasta \
        --taxid 2 \
        --include-children \
        --report only.bac.report
```


For the mitogenome assembly with pacbio long read data the reference genome used is the assembled mitogenome for thermobia in ncbi.
```
#!/bin/bash
#SBATCH -p serial_requeue # Partition to submit to (comma separated)
#SBATCH -J mitohifi1
#SBATCH -c 36 
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=50gb
#SBATCH --open-mode=append
#SBATCH -t 08:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out


singularity exec --cleanenv /n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/mitohifi_master.sif mitohifi.py \
-r "m64408e_230104_191438.hifi_reads.filt.fastq.gz \
m64408e_230106_060914.hifi_reads.filt.fastq.gz \
m64408e_230107_170457.hifi_reads.filt.fastq.gz \
m64408e_230109_021420.hifi_reads.filt.fastq.gz" \
-f T.dom.mito.NC_006080.1.fasta -g Thermobia_domestica_mitochondrion_complete_genome_NC_006080.1.gb \
-t 36 -o 5

```

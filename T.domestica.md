Assembly process for firebrat's genome.
This is just the collection of all the scripts I have used various times to process the assembly. This is not a pipeline.

Received the data from Pacbio hifi as bam files.
To convert these bamfiles into fastq, one can use samtools
'''
#!/bin/bash
#SBATCH -p test  # Partition to submit to (comma separated)
#SBATCH -J samtools_1
#SBATCH -c 10 
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=2gb
#SBATCH --open-mode=append
#SBATCH -t 00:20:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

module load samtools/1.10-fasrc01

for i in $(find ../1.raw_reads_pacbio_Tdom/ -type f -name "*.bam"); do
        filename=$(basename "$i")
        samtools fastq -@ 10 "$i" > "${filename%.*}.fastq"
done
'''
Check the quality of data with fastqc
```
#!/bin/bash
#SBATCH -p test #serial_requeue # Partition to submit to (comma separated)
#SBATCH -J fastqc_1
#SBATCH -c 10 
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=1gb
#SBATCH --open-mode=append
#SBATCH -t 01:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

module load fastqc/0.11.8-fasrc01

fastqc ../m64408e_230104_191438.hifi_reads.fastq -o ./
```
But converting bam files to fastq with samtools was not a part of my pipeline.
I used longqc to check the data quality in bamfiles
```
#!/bin/bash
#SBATCH -p test #serial_requeue # Partition to submit to (comma separated)
#SBATCH -J longqc_1
#SBATCH -c 10 
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=15gb
#SBATCH --open-mode=append
#SBATCH -t 05:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

export PATH="/n/home12/upendrabhattarai/bin/long-qc/LongQC:$PATH"
export PATH="/n/home12/upendrabhattarai/bin/miniconda3/bin:$PATH"
python /n/home12/upendrabhattarai/bin/long-qc/LongQC/longQC.py sampleqc --index 200M -p 10 -d -x pb-sequel -o 1_A01_longqc ../1_A01/m64408e_230104_191438.hifi_reads.bam
```

I then used Hifiadapterfilter to trimout the adapters
```
#!/bin/bash
#SBATCH -p test #serial_requeue # Partition to submit to (comma separated)
#SBATCH -J hifi.adapt.filt.1
#SBATCH -c 6 
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=10gb
#SBATCH --open-mode=append
#SBATCH -t 05:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

module load bamtools/2.3.0-fasrc01
module load GCC/8.2.0-2.31.1 OpenMPI/3.1.3 BLAST+/2.9.0
module load GCCcore/8.2.0 pigz/2.4

export PATH=$PATH:/n/home12/upendrabhattarai/bin/HiFiAdapterFilt
export PATH=$PATH:/n/home12/upendrabhattarai/bin/HiFiAdapterFilt/DB

bash pbadapterfilt.sh -p ../1_A01/m64408e_230104_191438 -t 6 -o 1_A_filt
```


I used hifiasm to assemble the reads with adapter removed reads

```
#!/bin/bash

#SBATCH -p serial_requeue # Partition to submit to (comma separated)
#SBATCH -J hifiasm_1
#SBATCH -c 32 
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=250gb
#SBATCH --open-mode=append
#SBATCH -t 3-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

export PATH="/n/home12/upendrabhattarai/bin/hifiasm:$PATH"
hifiasm -o TdomV1.asm -t32 ../reads/*.fastq.gz
```

I needed to remove mitogenome from the assembly. Mitogenome has been already assembled for T.domestica. I downloaded that from the NCBI in genebank format.
I installed MitoFinder from Github, the singularity and docker installation of MitoHiFi didn't work properly so had to install clonning the github repo.

`MitoFinder.sh`
```
#!/bin/bash
#SBATCH -p test #serial_requeue # Partition to submit to (comma separated)
#SBATCH -J mitohifi1
#SBATCH -c 16 
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=50gb
#SBATCH --open-mode=append
#SBATCH -t 08:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

/n/holyscratch01/extavour_lab/Everyone/upendrabhattarai/bin/MitoFinder/mitofinder \
-j TdomMito \
-a ../../adapter.removed/hifiasm.v1/TdomV1.asm.bp.hap1.p_ctg.fasta \
-r ../Thermobia_domestica_mitochondrion_complete_genome_NC_006080.1.gb \
-p 16 -m 50Gb -o 1
```

Mapping illumina paired end reads to the genome.

```
#!/bin/bash

#SBATCH -p test #shared  # Partition to submit to (comma separated)
#SBATCH -J bwa.mem
#SBATCH -c 16
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=20gb
#SBATCH --open-mode=append
#SBATCH -t 04:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out


#module load bwa/0.7.17-fasrc01

#bwa index TdomV1.asm.bp.hap1.p_ctg.fasta

#bwa mem -t 20 TdomV1.asm.bp.hap1.p_ctg.fasta SRR7665252_trim_1.fastq SRR7665252_trim_2.fastq > TdomV1.SRR.aln-pe.sam
module load intel/2017c impi/2017.4.239 SAMtools/1.9
#samtools sort -n -O sam TdomV1.SRR.aln-pe.sam | samtools fixmate -m -O bam - TdomV1.SRR.aln-pe.sorted.fixmate.bam
#samtools sort -O bam -o TdomV1.SRR.aln-pe.fixmate.sorted.bam TdomV1.SRR.aln-pe.sorted.fixmate.bam --threads 16
#samtools markdup -r -S TdomV1.SRR.aln-pe.fixmate.sorted.bam TdomV1.SRR.aln-pe.fixmate.sorted.dedup.bam --threads 16
samtools flagstat TdomV1.SRR.aln-pe.fixmate.sorted.dedup.bam --threads 16
```
Estimating the genome size with illumina reads
```
#!/bin/bash

#SBATCH -p serial_requeue # Partition to submit to (comma separated)
#SBATCH -J jellyfish
#SBATCH -c 6
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=100gb
#SBATCH --open-mode=append
#SBATCH -t 03:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out


module load intel/2017c impi/2017.4.239 Jellyfish/2.2.10

jellyfish count -C -m 21 -s 1000000000 -t 10 *.fastq -o reads.jf
#jellyfish count -t 6 -C -m 21 -s 20G -o 21mer_out --min-quality=30 SRR7665252_trim_1.fastq SRR7665252_trim_2.fastq
```
Pbmm2 alignment
```
#!/bin/bash

#SBATCH -p serial_requeue # Partition to submit to (comma separated)
#SBATCH -J pbmm2.1
#SBATCH -c 20
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=80gb
#SBATCH --open-mode=append
#SBATCH -t 3-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

export PATH="/n/home12/upendrabhattarai/bin/miniconda3/bin:$PATH"
pbmm2 align --preset CCS \
        ../../adapter.removed/hifiasm.v1/TdomV1.asm.bp.hap1.p_ctg.fasta \
        ../../1_A01/m64408e_230104_191438.hifi_reads.bam \
        pbmm2.1_A01.aligned.all.bam \
        --sort -j 20 -J 20
```

merge all the bam files and index it
```

#SBATCH -p serial_requeue # Partition to submit to (comma separated)
#SBATCH -J pbtkmerge
#SBATCH -c 10
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=50gb
#SBATCH --open-mode=append
#SBATCH -t 1-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

#export PATH="/n/home12/upendrabhattarai/bin/miniconda3/bin:$PATH"
#pbmerge -o merged.bam \
#../1/pbmm2.1_A01.aligned.all.bam \
#../2/pbmm2.2_B01.aligned.all.bam \
#../3/pbmm2.3_C01.aligned.all.bam \
#../4/pbmm2.4_D01.aligned.all.bam

#pbindex merged.bam -j 10
# pbindex produced .pbi file not .bai file as index so using samtool for this

module load GCC/8.2.0-2.31.1 SAMtools/1.9
#samtools sort merged.bam -o merged.sorted.bam -@ 10

samtools index merged.sorted.bam -@ 10
```
cpg calling
```
#!/bin/bash

#SBATCH -p serial_requeue # Partition to submit to (comma separated)
#SBATCH -J cpg.call
#SBATCH -c 16
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=100gb
#SBATCH --open-mode=append
#SBATCH -t 3-00:00:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out


export PATH="/n/home12/upendrabhattarai/bin/miniconda3/envs/cpg/bin:$PATH"
#export PATH="/n/home12/upendrabhattarai/bin/pb-CpG-tools:$PATH"

python /n/home12/upendrabhattarai/bin/pb-CpG-tools/aligned_bam_to_cpg_scores.py \
-b merged.sorted.bam \
-f ../../adapter.removed/hifiasm.v1/TdomV1.asm.bp.hap1.p_ctg.fasta \
-o 1_A01 -t 16 -d /n/home12/upendrabhattarai/bin/pb-CpG-tools/pileup_calling_model/
```

### To separate individual chromosomes
```
grep -w -A 1 ">TdomScaff_2" TdomHap2.FINAL.length.sorted.renamed.oneliner.fa > individual.chromosomes/TdomScaff_2.fa
```
### split each chromosome scaffolds in to contigs (break it at "N")
```
cat ../TdomH1Scaff_18.fa | perl -ne 'if(/^\>/){$scafnum++;}else{my $len=length($_);my @scaftigs=split(/N+/i,$_);my $scaftignum=0;foreach my $scaftig(@scaftigs){ my $len=length($scaftig);$scaftignum++;  print ">TdomH1Scaff_18:";print ""; print "$scaftignum-$len\n$scaftig\n";}}' > contig_chr_18.fa
```

Script to run purgehaplotigs
```
#!/bin/bash

#SBATCH -J purge.3
#SBATCH -p test #serial_requeue #test #shared #serial_requeue # Partition to submit to (comma separated)
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem 15gb
#SBATCH --open-mode=append
#SBATCH -t 00-01:00
#SBATCH -e %x_%j.err
#SBATCH -o %x_%j.out

#for softwares installed with conda
export PATH=/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/mambaforge/envs/purge_haplotigs/bin:$PATH

GENOME=Eann.GCA_023213315.1.filtered.fa
READS=../sRA/SRR178575*

#minimap2 -t 20 -ax map-ont $GENOME $READS --secondary=no \
#       | samtools sort -@ 20 -o aligned.bam -T tmp.ali

#purge_haplotigs hist -b aligned.bam -g $GENOME -t 20

#purge_haplotigs cov -i aligned.bam.gencov -l 1 -m 25 -h 65 -o coverage_stats.csv -j 80 -s 80

purge_haplotigs purge -g $GENOME -c coverage_stats.csv -t 20
```


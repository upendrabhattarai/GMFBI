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
We followed the script here to run transdecoder after pasa
https://github.com/TransDecoder/TransDecoder/blob/master/sample_data/pasa_example/runMe.sh

```
singularity exec --cleanenv /n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/transdecoder_latest.sif \
        TransDecoder.LongOrfs -t ../T.dom.pasa.sqlite.assemblies.fasta \
        --gene_trans_map pasa.gene_trnas_map.txt -O pasa.transdecoder_workdir

```
## 6. Braker2


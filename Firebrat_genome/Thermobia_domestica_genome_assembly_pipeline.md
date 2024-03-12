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
The Multiqc report after adapter filtering is [here](GMFBI/Firebrat_genome
/multiqc_report_hifiadapterfilt.html)


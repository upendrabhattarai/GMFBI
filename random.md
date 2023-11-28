To load R and 
```
module load R/4.3.1-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.3.1-fasrc01:$R_LIBS_USER
```

Get reads length of fasta file
```
cat input.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > read.length.txt
```

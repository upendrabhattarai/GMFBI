Installing transposonPSI
need to install blastall or blast legacy to run formatdb and other blast commands used by older version of blast. Those commands are modified in newer blast releases like formatdb to makeblastdb.
Conda installation of blast legacy is an easy solution.
activate conda environment and install
```
conda install -c bioconda blast-legacy
```

To download GyDB collection for ltrharvest
```
wget https://gydb.org/extensions/Collection/collection/db/GyDB_collection.zip
```

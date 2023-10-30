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

I updated the Dfam.h5 library that comes with Repeatmasker with the version from 10-Jan-20223 20:07 This is version `Dfam_3.7`
Its a big file 83.5 Gb so takes time to download
download link `https://www.dfam.org/releases/Dfam_3.7/families/`

```
wget https://www.dfam.org/releases/Dfam_3.7/families/Dfam.h5.gz
```

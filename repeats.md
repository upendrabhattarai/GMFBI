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

I wanted to updated the Dfam.h5 library that comes with Repeatmasker with the version from 10-Jan-20223 20:07 This is version `Dfam_3.7`
Its a big file 83.5 Gb compressed. I couldn't wait for the download to finish, but if needed update we can download it from here.

```
wget https://www.dfam.org/releases/Dfam_3.7/families/Dfam.h5.gz
```

I updated the Repeatmasker library with Repbase library
First make copy of original Repeatmasker library
```
/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/tetools_latest.sif cp -r /opt/RepeatMasker/Libraries/ ./
```
Download the Repeatmasker edition of Repbase library in the same directory and extract it there it will add into the RepeatMasker library you copied there with Repbase file.
Then run the following code to merge the database:
```
/n/holylabs/LABS/extavour_lab/Users/upendrabhattarai/software/tetools_latest.sif addRepBase.pl -libdir Libraries/
```
Now you can specify this Libraries/ directory by setting the LIBDIR environment variable,


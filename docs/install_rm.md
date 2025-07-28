# Install RepeatMasker (RM)

1. Download Sequence Search Engine (e.g. HMMER)
```
wget http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz
tar -xvf hmmer-3.4.tar.gz
cd hmmer-3.4 
./configure
make -j 8
cd easel
make -j 8
cd ../..
rm hmmer-3.4.tar.gz

HMMER_DIR=$(realpath hmmer-3.4/src)

```
2. Download DFAM (8.9 GB - takes about 1 hour to download) and extract (67 GB - takes about 15 mins to extract).
Skip if you already have downloaded a copy of the Dfam database.
Then just set `DFAM_FILE` variable to the correct full path of the database file.
```
wget https://www.dfam.org/releases/Dfam_3.8/families/FamDB/dfam38-1_full.0.h5.gz
gunzip dfam38-1_full.0.h5.gz

DFAM_FILE=$(realpath dfam38-1_full.0.h5)
```
3. Download RepeatMasker.
   Discard the error `Can't open DateRepeats: No such file or directory.` in the configuration step below.

   
```
wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.6.tar.gz
tar -xvzf RepeatMasker-4.1.6.tar.gz
cd RepeatMasker
mv $DFAM_FILE Libraries/famdb
./configure -trf_prgm $TRF_BINARY -hmmer_dir $HMMER_DIR -default_search_engine hmmer
cd ..
rm RepeatMasker-4.1.6.tar.gz
```

4. Make sure to run RepeatMasker on `test/sample.fa` to build libs and indices that are required for future workflows.
```
REPEAT_MASKER=$(realpath RepeatMasker/RepeatMasker)
$REPEAT_MASKER test/sample.fa
```
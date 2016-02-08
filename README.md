# pulse
predicting protein stability

Change `preprocess/preprocess_settings.json` to your desired settings before use.

### DEPENDENCIES
* preprocessing requires cufflinks and samtools
* BLAST is required in the preprocessing step after alternative splice site extraction
* last step of the preprocessing step requires the [biopython package](https://github.com/biopython/biopython)
* disorderome scoring step of feature extraction step requires iupred
* domain scoring requires pfam
* perl is required for some scripts in feature extraction and in last step
* sable is required in the last step of feature extraction

### SETUP

To create the necessary directories for output from pipeline:

```
python setup.py
```

Use `readlink -e .` when you're in root, and export that variable to PULSE_PATH using:
```
export PULSE_PATH=<result from readlink -e .>
```

### TODO:
* Need to download the uniprot using API and insert into pulse/input. After the uniprot_sprot file is downloaded into input/, need to run:
```
formatdb -i ./input/uniprot_sprot.fasta -p T
```
This creates a hashed version of uniprot_sprot.fasta database
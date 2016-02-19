# pulse
predicting protein stability.
This repo is a refactoring of @ccorbi's main pulse pipeline for my needs.

Change `preprocess/preprocess_settings.json` and `features/features_settings.json` to your desired settings before use.

### DEPENDENCIES
* python2.7
* preprocessing requires cufflinks and samtools
* BLAST is required in the preprocessing step after alternative splice site extraction
* last step of the preprocessing step requires the [biopython package](https://github.com/biopython/biopython)
* disorderome scoring step of feature extraction step requires iupred
* domain scoring requires [pfam] (http://www.vcru.wisc.edu/simonlab/bioinformatics/programs/install/pfamscan.htm) 
* perl is required for some scripts in feature extraction and in last step
* perl needs the following modules (download with `cpan <module name>`:
    * Getopt::Long
    * LWP::UserAgent
    * HTTP::Request::Common qw(POST), qw(GET)
    * HTTP::Headers
    * XML::XPath
    * XML::XPath::XMLParser
    * JSON
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

### USAGE

First, insert all of your .bam files into ```input/cell_lines```, and then run:  

```
python app.py
```

### PIPELINE

_Preprocessing_

1. Cufflinks generates transcriptome from sequence alignment data (.bam file).
2. Find exon-skipping alternative splicing events from transcriptome.
3. Identify anchors to alternative splicing events using BLAST.
4. Filter out alternative splicing events without a good anchor.
5. Generate an index file to locations of alternative splcing events.

_Feature Extraction_

1. Transmembrane scoring.
2. Post transcriptional modification scoring.
3. Eukaryotic linear motif scoring.
4. Disorderome scoring.
5. Domain scoring using PFAM.
6. Sable scoring.
7. Mutation scoring.
8. Conservation/network features generation.

_Machine Learning Step_

Five iterations of random forest classifier - takes the average of all iterations to generate final votes from ensemble and importance of each feature.

### TODO:
* Need to download the uniprot using API and insert into pulse/input. After the uniprot_sprot file is downloaded into input/, need to run:
```
formatdb -i ./input/uniprot_sprot.fasta -p T            # creates a hashed version of uniprot_sprot.fasta database
```

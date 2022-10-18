# Atlas baseline meta-analysis

## Metadata

To merge the metadata of datasets, we use within a Snakemake workflow the [MAGE-Tab-merger](https://pypi.org/project/MAGE-Tab-merger/) package to produce:

- merged condensed SDRF file
- merged SDRF file
- merged XML configuration file

from an starting list of Atlas baseline experiment accessions. The process of dealing with the metadata will filter out datasets that cannot be merged due to experimental design or metadata limitations.

To run the metadata part you will need to specify the following environment variables:

```
# The accessions of experiments that you want to merge
export ACCESSIONS="E-MTAB-4395,E-MTAB-4342,E-MTAB-4128"
# The accession to be used for the new experiment
export NEW_ACCESSION="E-CORN-1"
# The path to cache where all data will be downloaded (if INPUT_PATH is not specified)
export CACHE_PATH="/absolute/path/to/cache"
# The name used to the batch category, usually "study"
export BATCH=study
# The covariate used to merge metadata, as written in the SDRF/Condensed
export COVARIATE='organism part'
# The type of SDRF field that the covariate is, usually characteristic, but could be a factor as well (less likely).
export COVARIATE_TYPE=characteristic
```

Optionally, you can set:

```
# If you need to work with data that is not available on the FTP
# or if you have manually tweaked files, setting this to where you have
# the data will override the retrieval of data from the public FTP.
export INPUT_PATH="/some/absolute/path/for/data"
# On occassions the covariate choosen might have values that are too
# general and that should be skipping for merging.
export COVARIATE_SKIP_VALUES='whole organism'
# If SEND_TO_ANALYSIS is set, then all data outputs will be generated as well. Takes more time to download and merge things, although not a lot more.
export SEND_TO_ANALYSIS=true
# If REPLACE_DOWNLOADS is set, then all downloaded files will be
# overwritten
export REPLACE_DOWNLOADS=true
```

### Directory structure when using INPUT_PATH

If setting `INPUT_PATH` then the pipeline won't attempt to download data from the public FTP, and will assume that all the input data is available on that directory. So, when setting this var, this is how data needs to be available (here it is entirely up to to operator to provision it), following the example for `$ACCESSIONS` above:

```
$ cd $INPUT_PATH
$ ls *
E-MTAB-4395:
E-MTAB-4395-configuration.xml                       E-MTAB-4395.condensed-sdrf.tsv
E-MTAB-4395-raw-counts.tsv.undecorated              E-MTAB-4395.idf.txt
E-MTAB-4395-transcripts-raw-counts.tsv.undecorated  E-MTAB-4395.sdrf.txt

E-MTAB-4342:
E-MTAB-4342-configuration.xml                       E-MTAB-4342.condensed-sdrf.tsv
E-MTAB-4342-raw-counts.tsv.undecorated              E-MTAB-4342.idf.txt
E-MTAB-4342-transcripts-raw-counts.tsv.undecorated  E-MTAB-4342.sdrf.txt
```

So, as a summary, each accesion in `$ACCESSIONS` is expected to have its own directory named after it, and each directory contain the raw-counts (genes and transcripts), IDF, SDRF and condensed SDRF, with the name structure shown above.

### Running the metadata merge workflow

```
cd desired/working/directory
# start snakemake environment - possibly conda activate snakemake
snakemake --snakefile <PATH_TO_THIS_REPO>/meta-data/Snakefile --cores 2 --use-conda --conda-frontend mamba --config accessions=<accession list> new_accession=<new accession> cache_path=<path for temp files> batch=<batch> covariate=<covariate> covariate_type=characteristic species=<species name> retrieve_data=True
```

Example run:

```
cd desired/working/directory
# start snakemake environment - possibly conda activate snakemake
snakemake --snakefile <PATH_TO_THIS_REPO>/meta-data/Snakefile --cores 2 --use-conda --conda-frontend mamba --config accessions=E-GEOD-53197,E-GEOD-55482,E-CURD-31,E-GEOD-52806,E-GEOD-64740 new_accession=E-SUPR-1 cache_path=cache batch=study covariate="organism part" covariate_type=characteristic species="arabidopsis_thaliana" retrieve_data=True
```

### Results

Results will be present on the working directory, inside `tmp_results/data/<new accession>`. In the maximal execution setting (with data), you should see given a `$NEW_ACCESSION` set to `E-SUPR-1`:

```
E-SUPR-1-configuration.xml                      E-SUPR-1.condensed.sdrf.tsv
E-SUPR-1-raw-counts.tsv.undecorated             E-SUPR-1.sdrf.tsv
E-SUPR-1-transcripts-raw-counts.tsv.undecorated	E-SUPR-1.selected_studies.txt
E-SUPR-1.gtf                                    E-SUPR-1-analysis-methods.tsv
```

Most files are self explanatory, the `selected_studies.txt` files contains a comma separate list of the original accessions used.


## Data Analysis

The current analysis setup enables to merge Atlas RNA-Seq baseline datasets that have the required compatibility in terms of metadata (done with the metadata workflow above).

This analysis requires:

- RAW counts for genes: at `data/{accesion}/{accesion}-raw-counts.tsv.undecorated`
- RAW counts for transcripts: at `data/{accesion}/{accesion}-transcripts-raw-counts.tsv.undecorated`
- GTF file for the organism, which matches what was used to generate counts: at `data/{accesion}/{accesion}.gtf`
- Merged configuration XML file with assays: at `data/{accesion}/{accesion}-configuration.xml`
- Merged SDRF file: at `data/{accesion}/{accesion}.sdrf.txt`
- Methods file, containing analysis methods used so far: at `data/{accesion}/{accesion}-analysis-methods.tsv`

## Running the analysis

Once you have the following files in place in a `data` directory (TODO generalise this) with the `ACCESSION` being the newly minted accession for the merged dataset:

```
data/<ACCESSION>/<ACCESSION>-analysis-methods.tsv
data/<ACCESSION>/<ACCESSION>-configuration.xml
data/<ACCESSION>/<ACCESSION>-raw-counts.tsv.undecorated
data/<ACCESSION>/<ACCESSION>-transcripts-raw-counts.tsv.undecorated
data/<ACCESSION>/<ACCESSION>.sdrf.tsv
data/<ACCESSION>/<ACCESSION>.gtf
```

then run:

```
snakemake -p --snakefile <PATH_TO_THIS_REPO>/data-analysis/Snakefile --cores 2 --use-conda --conda-frontend mamba \
  --config accession=<ACCESSION> cache_path=<path for temp files> batch=<batch>
```

example run:

```
snakemake -p --snakefile <PATH_TO_THIS_REPO>/data-analysis/Snakefile --cores 2 --use-conda --conda-frontend conda --config accession=E-SUPR-1 cache_path=tmp_results batch=Study
```

This will create results in two directories:

- tmp_results: main results, with all calculations for the merged dataset.
- lengths: where lengths per gene, exon and transcripts are produced, based on the reference GTF file provided. These are used to normalise to TPMs and FPKMs.

# Atlas baseline meta-analysis
The repository contains workflows for metadata generation and data analysis.
## Metadata

To merge the metadata of datasets, we use within a Snakemake workflow the [MAGE-Tab-merger](https://pypi.org/project/MAGE-Tab-merger/) package to produce:

- merged condensed SDRF file
- merged SDRF file
- merged XML configuration file

from an starting list of Atlas baseline experiment accessions. The process of dealing with the metadata will filter out datasets that cannot be merged due to experimental design or metadata limitations.


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
snakemake --snakefile <PATH_TO_THIS_REPO>/meta-data/Snakefile --cores 2 --use-conda --conda-frontend mamba --config accessions=E-GEOD-53197,E-GEOD-55482,E-CURD-31,E-GEOD-52806,E-GEOD-64740 new_accession=E-SUPR-1 cache_path=cache batch=study covariate="organism part" covariate_type=characteristic species="arabidopsis_thaliana" retrieve_data=True output_path=output
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

### Jenkins job

Jenkins job for this part of workflow is at http://gene-expression.ebi.ac.uk/jenkins/job/bulk_baseline_meta_analysis_dataprep/ 


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
  --config accession=<ACCESSION> output_path=<path for temp files> batch=<batch> covariate=<covariate>
```

example run:

```
snakemake -p --snakefile <PATH_TO_THIS_REPO>/data-analysis/Snakefile --cores 2 --use-conda --conda-frontend conda --config accession=E-SUPR-1 output_path=tmp_results batch=Study covariate=organism_part
```

This will create results in two directories:

- tmp_results: main results, with all calculations for the merged dataset.
- lengths: where lengths per gene, exon and transcripts are produced, based on the reference GTF file provided. These are used to normalise to TPMs and FPKMs.

### Jenkins job

Jenkins job for this part of workflow is at http://gene-expression.ebi.ac.uk/jenkins/job/bulk_baseline_meta_analysis/

# Known bug

## Covariate

Meta-data pipeline needs covariate without underscore, i.e. organism part
Data analysis pipeline needs covariate with underscore, i.e. organism_part 

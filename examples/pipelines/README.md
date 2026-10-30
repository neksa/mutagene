# Run MutaGene in a workflow

These examples share the labelled SBS96 matrix interface.
Each workflow runs the offline demonstration and fits COSMICv2 signatures.
The fixture contains 24 mutations across three synthetic samples.
It demonstrates execution, not reliable signature attribution.

## Prepare the environment

These examples require the updated source checkout.
The older PyPI release does not contain the new commands.
From the repository root, install the locked environment and activate it:

```sh
uv sync --frozen
source .venv/bin/activate
mutagene --version
```

Install the workflow engine separately.
Use a fresh working directory for each engine.
Record the repository commit, `uv.lock`, engine version, and input hashes with each analysis.
Keep engine dependencies outside the core MutaGene installation.

## DVC

Copy `dvc/dvc.yaml` into a separate, empty Git repository.
Initialize DVC and execute the pipeline:

```sh
git init
dvc init
dvc repro
```

Commit `dvc.yaml` and `dvc.lock` with your environment specification.
Configure your own DVC remote before sharing cached outputs.
The example does not configure or upload to a remote.

## Nextflow

From a fresh working directory, use the absolute path to `nextflow/main.nf`:

```sh
nextflow run /path/to/mutagene/examples/pipelines/nextflow/main.nf --outdir results
```

Nextflow copies `signatures.tsv` and `diagnostics.json` into `results`.
Use `-resume` to reuse completed tasks.
The example uses the local executor and the active `mutagene` executable.

## Snakemake

Copy `snakemake/Snakefile` into a fresh directory and execute:

```sh
snakemake --cores 1
```

## Replace the demonstration

Generate the same matrix contract from your mutation files:

```sh
mutagene profile --per-sample -i cohort.maf.gz -g hg38 -o matrix.tsv
mutagene signature -f MATRIX -i matrix.tsv -o signatures.tsv --diagnostics-out diagnostics.json
```

Each row identifies one SBS96 channel, such as `A[C>A]A`.
The first column has the name `Type`.
Each remaining column contains integer counts for one uniquely named sample.
The reader maps channel labels; it does not assume row order.
Split multisample VCF files before analysis.
Provide the assembly that matches the input coordinates.

## Public repository preparation

These examples provide a starting point for separate public workflow repositories.
Before publication, pin a released MutaGene version and add engine-specific CI.
Add small fixtures, expected outputs, a license, and citation metadata.
A future release can supply a container digest after image validation.
No public workflow repository accompanies this draft.

## References

- [DVC pipeline definition](https://doc.dvc.org/user-guide/project-structure/dvcyaml-files)
- [Nextflow processes](https://docs.seqera.io/nextflow/process)
- [Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html)

## Validation records

`validation.json` records installed-wheel execution with the preceding MLE default.
`validation-sparse-default.json` records execution with the sparse-BIC default proposal.
Each record identifies the method, engine versions, wheel hash, and output hashes.

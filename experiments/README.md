# Experiment setup

This directory holds the scripts that implement the MetaEdit2 workflow. For a description of the numbered experiments see the repository README.

## Software requirements

- `eggnog-mapper`
- Python packages: `pandas`
- Singularity 3+

Create a conda environment and install the dependencies:

```bash
conda create -n metaedit2 -c conda-forge -c bioconda eggnog-mapper pandas
conda activate metaedit2
```

## Singularity container

The `env/` folder contains definition files for building a Singularity container with all remaining tools. Build and launch it with:

```bash
module load singularity
bash ../env/build_singularity.sh  # creates metaEdit.sif
bash ../env/run_singularity.sh   # enter the container
```

All `submit_*.sh` scripts assume they are executed inside this container and use SLURM for job scheduling.

## Running the workflow

Provide sample identifiers in `data/samples.txt` and execute the numbered scripts in order (`2-download_dna.sh`, `3-download_rna.sh`, ...). The submission scripts produce log files in the `logs/` directory.

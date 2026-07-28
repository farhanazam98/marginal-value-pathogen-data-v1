# marginal-value-pathogen-data-v1

**Objective:** Measuring whether removing sequences >90% identical to the SARS-CoV-2 RBD changes PSSM's ability to predict the Starr 2020 ACE2 binding assay, across cumulative database snapshots from 2010 to 2025, with a depth-matched control.

## Setup

All scripts are pure-Python stdlib; the one external dependency is the
[`ncbi-datasets-cli`](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)
(`datasets` / `dataformat` binaries), distributed via conda. On a fresh machine (e.g. a new
EC2 instance):

```bash
./setup.sh
conda activate marginal-value-pathogen-data
```

`setup.sh` installs Miniconda if it isn't already present, accepts the conda channel
Terms of Service (required non-interactively), and creates/updates the
`marginal-value-pathogen-data` environment from `environment.yml`. Re-running it is safe
and just updates the existing environment.

## Milestones:
### Monday
Come up with objective. Find and download, and the Starr 2020 DMS data for SARS-COV-2 in the [priority-viruses](https://github.com/debbiemarkslab/priority-viruses) repo and parse the N501Y mutation. 

*Deliverable: CSV file with DMS data from Starr 2020 experiment for SARS-COV-2 under data/ and a script that parses the N501Y mutation*

### Tuesday
Find viral protein database for SARS-COV-2 sequences, segmented by creation date and host field (if available). 

*Deliverable: All SARS-COV-2 protein sequences (or at least a subset) downloaded, segmented by creation date, and a table filled out with actual numbers. We should have an idea how long it will take to download all the data if possible, or how to fetch it otherwise.* 


### Wednesday
For N501Y, build MSA, calculate mutation effect prediction, and compare against DMS score. 

### Thursday






# PacBio Pipeline

## Overview

Pipeline for processing PacBio data (*subreads.bam) with PacBio SMRT command line tools.

## Installation

Download and install SMRT Link from PacBio.

```
wget -P <INSTALL_DIR> https://downloads.pacbcloud.com/public/software/installers/smrtlink_10.1.0.119588.zip
cd <INSTALL_DIR>
unzip smrtlink_*.zip
smrtlink_*.run --rootdir smrtlink --smrttools-only
```

Make sure you run in an environment with snakemake version >= 5 < 6 and conda installed.

```
module load snakemake/v5.4.4
```

## Usage

Export the path to SMRT Link binary.

```
export PATH=<INSTALL_DIR>/smrtlink/smrtcmds/bin/:$PATH
```

Edit the config.yml and enter the path to your reference genome and samples. All files of each sample must the within a separate folder, the name of the folder is assumed to be equal to the the name of the sample.

```
snakemake all --cores $(nproc)
```





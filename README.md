
# [Spt6 2018](https://doi.org/10.1101/347575) small data

## description

Small data sets used in [our preprint](https://doi.org/10.1101/347575). Includes the following:

- histone H3 ChIP-qPCR (`h3-chip-qpcr`)
- RT-qPCR during Spt6 depletion (`spt6-depletion-rtpqcr`)
- Spt6 Western blots (`spt6-western`)
- TFIIB ChIP-qPCR (`tfiib-chip-qpcr`)
- TFIIB Western blots (`tfiib-western`)

## requirements

- Unix-like operating system (tested on CentOS 7.2.1511)
- Git
- [conda](https://conda.io/docs/user-guide/install/index.html)

## instructions
**0**. To create figures for any of the ChIP-qPCR datasets or TFIIB Western blots, create and activate the `small_data` virtual environment using conda. This can take a while, be patient.

```bash

# create the small_data environment
conda env create -v -f envs/small_data.yaml

# activate the environment
source activate small_data

# to deactivate the environment
# source deactivate
```

**1**. With the `small_data` environment activated, navigate to the relevant directory (`h3-chip-qpcr`, `tfiib-chip-qpcr`, or `tfiib-western`) and run the R script in the directory.

```bash
cd tfiib-chip-qpcr

Rscript tfiib-qpcr.R
```

**2.** To generate png versions of the svg files generated, run `mogrify.sh`.

```bash
./mogrify.sh
```



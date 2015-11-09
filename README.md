# bam-toolbox

## Installation

```
pip install https://github.com/AndersenLab/bam-toolbox/archive/0.0.1.tar.gz
```

## Usage

    bam-toolbox 0.1

    usage:
      bam <command> [<args>...]
      bam -h | --help
      bam --version

    commands:
      coverage
      readgroups
      fastq

## Commands

### Coverage

Calculate depth of coverage! There are several ways to do so:

1. Calculate coverage for the whole genome, nuclear genome, and individual chromosomes. Also calculated the mtDNA:nuclear ratio which acts as a proxy for mitochondrial content.
2.

```
usage:
  bam coverage <bam> [options] [--mtchr=<mtchr>]
  bam coverage <bam> [options] <chrom:start-end>...
  bam coverage <bam> [options] --window=<size>
  bam coverage <bam> [options] --regions=<gff/bed>
```

###### Parameters

* `--header` - Output the header line.
* `--tsv` - Output in `tsv` format rather than `eav` format.
* `--mtchr` - Specifying a mitochondrial chromosome will output the mtDNA:nuclear DNA ratio, offering an estimate of mitochondrial content. 
* `--window` - Calculate depth of coverage for given windows across the genome.
* `--regions` - Specify a bed or gff file and calculate coverage in those regions.

###### Example

```
bam coverage test.bam --tsv --header
```

__Result:__

Output is in `--tsv` format.

```
ENTITY  chrom   start   end ATTR    VALUE   DATE
CB4853  I   1   15072434    bases_mapped    1688210 2015-11-08 23:22:39.498133
CB4853  I   1   15072434    depth_of_coverage   0.112006461597  2015-11-08 23:22:39.498242
CB4853  I   1   15072434    breadth_of_coverage 0.0918075342045 2015-11-08 23:22:39.498328
CB4853  I   1   15072434    length  15072434    2015-11-08 23:22:39.498372
CB4853  I   1   15072434    pos_mapped  1383763 2015-11-08 23:22:39.498413
CB4853  II  1   15279421    bases_mapped    1605336 2015-11-08 23:22:39.923802
CB4853  II  1   15279421    depth_of_coverage   0.105065237747  2015-11-08 23:22:39.923880
...
CB4853  mt_nuclear_ratio    51.7503772943   2015-11-08 23:22:41.949819
```

__Formatted Result:__

| ENTITY   | chrom   |   start |      end | ATTR                |       VALUE | DATE                       |
|:---------|:--------|--------:|---------:|:--------------------|------------:|:---------------------------|
| CB4853   | I       |       1 | 15072434 | bases_mapped        | 1.68821e+06 | 2015-11-08 23:26:41.891988 |
| CB4853   | I       |       1 | 15072434 | depth_of_coverage   | 0.112006    | 2015-11-08 23:26:41.892087 |
| CB4853   | I       |       1 | 15072434 | breadth_of_coverage | 0.0918075   | 2015-11-08 23:26:41.892155 |
| CB4853   | I       |       1 | 15072434 | length              | 1.50724e+07 | 2015-11-08 23:26:41.892198 |
| CB4853   | I       |       1 | 15072434 | pos_mapped          | 1.38376e+06 | 2015-11-08 23:26:41.892239 |
| CB4853   | II      |       1 | 15279421 | bases_mapped        | 1.60534e+06 | 2015-11-08 23:26:42.319937 |

### fastq

Generate stats from fastq files. This makes use of a great [bioinformatics one liner](https://github.com/stephenturner/oneliners)!

__Command__

```
bam fastq --header --tsv test.fq.gz
```

__Result__

| ENTITY     | ATTR                    | VALUE                                                                                                |
|:-----------|:------------------------|:-----------------------------------------------------------------------------------------------------|
| test.fq.gz | Total_Reads             | 2500                                                                                                 |
| test.fq.gz | Unique_Reads            | 2500                                                                                                 |
| test.fq.gz | Percent_Unique          | 100                                                                                                  |
| test.fq.gz | Most_Abundant_Sequence  | CCCACCCCTAAATCATATGCCTAATCAAACTTTTGAATTAAACATTTAGGTAAACACAACTTGAATACGTTTTGAAATTCTGCTCTAAATTTGTTTTTGT |
| test.fq.gz | Most_Abundant_Frequency | 1                                                                                                    |
| test.fq.gz | Percentage_Unique_fq    | 0.04                                                                                                 |
| test.fq.gz | cum_length              | 237908                                                                                               |
| test.fq.gz | A_count                 | 75456                                                                                                |
| test.fq.gz | T_count                 | 74248                                                                                                |
| test.fq.gz | C_count                 | 43765                                                                                                |
| test.fq.gz | G_count                 | 44439                                                                                                |
| test.fq.gz | N_count                 | 0                                                                                                    |
| test.fq.gz | bases                   | 237908                                                                                               |
| test.fq.gz | GC_count                | 0.370748357                                                                                          |
| test.fq.gz | min_length              | 30                                                                                                   |
| test.fq.gz | avg_length              | 95.1632                                                                                              |
| test.fq.gz | max_length              | 100                                                                                                  |

## To Do

###### Coverage

* [X] Read from stdin / rework function.
* [X] Coverage across intervals/sliding windows
* [ ] Coverage plot

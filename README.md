# bam-toolbox



## Usage

    bam-toolbox 0.1

    usage:
      bam <command> [<args>...]
      bam -h | --help
      bam --version

    commands:
      coverage

## Commands

### Coverage

Calculate depth of coverage 

```
usage:
  bam coverage header [--eav]
  bam coverage <bam> [--mtchr=<mtchr>] [--eav]
```

###### Parameters

* `--header` - Output the header line. Use `--eav` to output entity-attribute-value header.
* `--eav` - Entity-Attributes-Value output format. Useful for combining results from multiple bams for comparison.
* `--mtchr` - Specifying a mitochondrial chromosome will output the mtDNA:nuclear DNA ratio, offering an estimate of mitochondrial content.

###### Example

```
bam coverage header --eav > bam_output.txt
bam coverage test.bam --eav >> bam_output.txt
```

__Result:__

Output is in `--eav` format.

```
ENTITY  ATTRIBUTES  VALUE   TIMESTAMP
test.bam    CONTIG=I;ATTR=Bases_Mapped  9890    2015-10-29 14:22:31.571553
test.bam    CONTIG=I;ATTR=Sum_of_Depths 401490  2015-10-29 14:22:31.571599
test.bam    CONTIG=I;ATTR=Length    15072434    2015-10-29 14:22:31.571612
test.bam    CONTIG=I;ATTR=Breadth_of_Coverage   0.000656164757464   2015-10-29 14:22:31.571621
test.bam    CONTIG=I;ATTR=Depth_of_Coverage 0.0266373699165 2015-10-29 14:22:31.571639
test.bam    CONTIG=II;ATTR=Bases_Mapped 10092   2015-10-29 14:22:31.571648
test.bam    CONTIG=II;ATTR=Sum_of_Depths    262301  2015-10-29 14:22:31.571656
test.bam    CONTIG=II;ATTR=Length   15279421    2015-10-29 14:22:31.571664
test.bam    CONTIG=II;ATTR=Breadth_of_Coverage  0.000660496232154   2015-10-29 14:22:31.571672
...
test.bam    CONTIG=genome;ATTR=mt_ratio 16788.2757429   2015-10-29 14:22:31.571884
```

__Formatted Result:__

| ENTITY   | ATTRIBUTES                         |            VALUE | TIMESTAMP                  |
|:---------|:-----------------------------------|-----------------:|:---------------------------|
| test.bam | CONTIG=I;ATTR=Bases_Mapped         |   9890           | 2015-10-29... |
| test.bam | CONTIG=I;ATTR=Sum_of_Depths        | 401490           | 2015-10-29... |
| test.bam | CONTIG=I;ATTR=Length               |      1.50724e+07 | 2015-10-29... |
| test.bam | CONTIG=I;ATTR=Breadth_of_Coverage  |      0.000656165 | 2015-10-29... |
| test.bam | CONTIG=I;ATTR=Depth_of_Coverage    |      0.0266374   | 2015-10-29... |
| test.bam | CONTIG=II;ATTR=Bases_Mapped        |  10092           | 2015-10-29... |
| test.bam | CONTIG=II;ATTR=Sum_of_Depths       | 262301           | 2015-10-29... |
| test.bam | CONTIG=II;ATTR=Length              |      1.52794e+07 | 2015-10-29... |
| test.bam | CONTIG=II;ATTR=Breadth_of_Coverage |      0.000660496 | 2015-10-29... |
| test.bam | CONTIG=genome;ATTR=mt_ratio        |  16788.3         | 2015-10-29... |

The result may look a little strange but it lends itself to performing analysis in programs like R. A simpler output is available byy ommitting the `--eav` flag.



## To Do

###### Coverage

* [ ] Read from stdin / rework function.
* [ ] Coverage across intervals/sliding windows
* [ ] Coverage plot

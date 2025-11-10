# Add structural annotations to dotplots

## The example:

Recently, Franco et al. identified clusters with several similar transposable elements (belonging to the *sSaTar* families 1-3) within several grasses. On the sorgum chromosome 1, sSaTar transposons are arranged in a linear tandem-manner, just separated by short microsatellites.
Combination of structural annotation with a dotplot, as possible with FlexiDot, allows visualisation of this peculiar region.

[Franco et al. (2018)](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4653-6) "Modular assembly of transposable element arrays by microsatellite targeting in the guayule and rice genomes". BMC Genomics 19:271

## FlexiDot illustration of this region:

<img src="https://github.com/flexidot-bio/flexidot/blob/master/tests/test-data/sSaTar_example/sSaTar_cluster_flexi_300b.png" width="600">

## Input files:

- [*sSaTar.fas*: sSaTar cluster on Sorghum chromosome 1](https://github.com/flexidot-bio/flexidot/blob/master/tests/test-data/sSaTar_example/sSaTar.fas)*
- [*sSaTar.gff3*: sSaTar annotations as gff3](https://github.com/flexidot-bio/flexidot/blob/master/tests/test-data/sSaTar_example/sSaTar.gff3)*
- [*sSaTar.config*: sSatar config file to define colors for FlexiDot](https://github.com/flexidot-bio/flexidot/blob/master/tests/test-data/sSaTar_example/sSaTar.config)

\* *fasta* and *gff3 files* have been deduced from Franco et al.'s [Supplemental File 14](https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-018-4653-6/MediaObjects/12864_2018_4653_MOESM14_ESM.pdf), showing the annotated sequence of this region.

The *config file* defines color, alpha and zoom of each sequence type. Please note, that the small microsatellite region is visualized with an additive zoom of `10`.

## Command:

```bash
# Set path to data in flexidot repo.
DATA="tests/test-data/sSaTar_example"

flexidot -i $DATA/sSaTar.fas -g $DATA/sSaTar.gff3 -G $DATA/sSaTar.config -k 10 -S 1 -T 30 -E 15 -A 2 -C black -f pdf
```

---

For additional application use cases, please see the [FlexiDot in-depth documentation (pdf)](https://github.com/flexidot-bio/flexidot/blob/master/docs/SupplementaryData.pdf).

Back to [FlexiDot home](https://github.com/flexidot-bio/flexidot).

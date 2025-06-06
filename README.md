[![Anaconda-Server Badge](https://anaconda.org/bioconda/flexidot/badges/version.svg)](https://anaconda.org/bioconda/flexidot)
[![License Badge](https://anaconda.org/bioconda/flexidot/badges/license.svg)](https://anaconda.org/bioconda/flexidot)

# FlexiDot: Highly customizable, ambiguity-aware dotplots for visual sequence analyses

![alt text](https://github.com/flexidot-bio/flexidot/blob/master/docs/images/Selfdotplots_banner4.png "FlexiDot self dotplots")

FlexiDot is a cross-platform dotplot suite generating high quality self, pairwise and all-against-all visualizations. To improve dotplot suitability for comparison of consensus and error-prone sequences, FlexiDot harbors routines for strict and relaxed handling of mismatches and ambiguous residues. The custom shading modules facilitate dotplot interpretation and motif identification by adding information on sequence annotations and sequence similarities to the images. Combined with collage-like outputs, FlexiDot supports simultaneous visual screening of a large sequence sets, allowing dotplot use for routine screening.

## Citation

If you use FlexiDot in your research, please cite us:

**Kathrin M. Seibt, Thomas Schmidt, and Tony Heitkam** (2018) "FlexiDot: Highly customizable, ambiguity-aware dotplots for visual sequence analyses". *Bioinformatics* 34 (20), 3575–3577, doi: 10.1093/bioinformatics/bty395  -  [**Read article**](https://doi.org/10.1093/bioinformatics/bty395)

## FlexiDot versions and updates

<img align="right" width="100" height="100" src="https://github.com/flexidot-bio/flexidot/blob/master/docs/images/FlexiLogo.png">

**Current version (May 2025): FlexiDot v2.0.0**

For an overview of FlexiDot version updates please see the [code history](https://github.com/flexidot-bio/flexidot/blob/master/CHANGELOG.md).

Corresponding [parameter cheat sheets](https://github.com/flexidot-bio/flexidot/tree/master/docs) are available as well.

## Documentation

* [In depth documentation](https://github.com/flexidot-bio/flexidot/blob/master/docs/SupplementaryData.pdf) (This readme gives an overview, and more detail is in the documentation.)
* [Parameter cheat sheet v1.0.6](https://github.com/flexidot-bio/flexidot/blob/master/docs/usage_v1.06.pdf)
* [Artificial test sequences used for the examples](https://github.com/flexidot-bio/flexidot/tree/master/tests/test-data)
* [Example: Adding annotation-based shading to a dotplot](https://github.com/flexidot-bio/flexidot/blob/master/docs/tutorial_add_annotation.md)
* [Presentation slides introducing dotplots and our FlexiDot tool](https://zenodo.org/record/2558556)

## Implementation

FlexiDot is implemented in [Python 3](https://www.python.org/), with dependencies:

* [numpy](https://pypi.python.org/pypi/numpy)
* [matplotlib](https://pypi.python.org/pypi/matplotlib)
* [biopython](https://pypi.python.org/pypi/biopython)
* [regex](https://pypi.python.org/pypi/regex)
* [colormap](https://pypi.python.org/pypi/colormap)
* [easydev](https://pypi.python.org/pypi/easydev) (required for colormap)
* [colour](https://pypi.python.org/pypi/colour)

You can create a Conda environment with these dependencies using the YAML file in this repo.

```bash
conda env create -f environment.yml

conda activate flexidot
```

After activating the flexidot environment you can use pip to install the latest version of Flexidot.

## Installing Flexidot

Installation options:

Install from PyPI (recommended):

```bash
pip install flexidot
```

Install from bioconda:

```bash
conda install -c bioconda flexidot
```

pip install the latest development version directly from this repo.

```bash
% pip install git+https://github.com/flexidot-bio/flexidot.git
```

Test installation.

```bash
# Print version number and exit.
flexidot --version

# Get usage information
flexidot --help
```

### Setup Development Environment

If you want to contribute to the project or run the latest development version, you can clone the repository and install the package in editable mode.

```bash
# Clone repository
git clone https://github.com/flexidot-bio/flexidot.git && cd flexidot

# Create virtual environment
conda env create -f environment.yml

# Activate environment
conda activate flexidot

# Install package in editable mode
pip install -e '.[dev]'

# Optional: Install pre-commit hooks
pre-commit install
```

## Use FlexiDot

Flexidot accepts one or more uncompressed fasta files as input. The files can contain multiple sequences.

```bash
# Use individual fasta file (can contain multiple sequences)
flexidot -i input.fasta [optional arguments]

# Use multiple fasta files
flexidot -i input1.fasta input2.fasta [optional arguments]

# Use all fasta files in current directory
flexidot -i *.fasta [optional arguments]
```

Optional arguments are explained below and in detail with the `--help` option.

Importantly, `-k` defines the word size (e.g. `-k 10`) and `-t` specifies the sequence type (`-t nuc` for DNA [default]; `-t aa` for proteins). The plotting mode is chosen via `-m` and described below.

## Plotting modes

FlexiDot allows sequence investigation in three run modes via the option `-m/--mode`:

`-m 0`    self sequence comparison
`-m 1`    pairwise sequence comparison
`-m 2`    all-to-all sequence comparison

To run multiple plotting modes, call the option multiple times i.e. `-m 0 -m 1 -m 2`.

### Self dotplots

with `-m/--mode 0`

In **self** dotplot mode, each sequence is compared with itself. The resulting dotplots can be combined to form a **collage** (with `--collage`) or written to separate files.

![alt text](https://github.com/flexidot-bio/flexidot/blob/master/docs/images/Selfdotplots_banner.png "FlexiDot self dotplots")

```bash
# A single sequence compared to itself
flexidot -i Seq2.fasta -m 0 -k 10 -P 15

# Single sequence with annotations
flexidot -i Seq2.fasta -m 0 -k 10 -P 15 -g example.gff3 -G gff_color.config

# Collage of 6 sequences each compared to themselves with Seq2 annotated (shown above)
flexidot -i test-seqs.fasta -m 0 -k 10 --n_col 6 -P 15 -g example.gff3 -G gff_color.config --collage
```

### Pairwise comparisons

with `-m/--mode 1`

For **pairwise** dotplots, the collage output is recommended for larger numbers of sequences. The collage output of the 15 pairwise dotplots for the test sequences is shown below. By default, dotplot images are in square format (panel A). This maximizes the visibility of matches, if the compared sequences differ drastically in length. To enable scaling according to the respective sequence lengths, the FlexiDot scaling feature is callable via option `-L/--length_scaling` (panel B). If scaling is enabled, a red line indicates the end of the shorter sequence in the collage output.

Pairwise comparisons can be limited to only pairs that contain the first sequence in a fasta file using `--only_vs_first_seq`.

<img src="https://github.com/flexidot-bio/flexidot/blob/master/docs/images/pairwise_low_res.png" width="600">

```bash
# Panel A
flexidot -i test-seqs.fasta -m 1 -k 10 --n_col 3 -c
# Panel B (with length scaling)
flexidot -i test-seqs.fasta -m 1 -k 10 --n_col 3 -c -L
```

### All-against-all comparisons

with `-m/--mode 2`

In **all-against-all** mode, FlexiDot compares each pair from a set of input sequences. To enable the identification of long shared subsequences at a glance, FlexiDot offers similarity shading (switched on/off via option `-x/--lcs_shading`) based on the LCS length in all-against-all comparisons (see below).

<img src="https://github.com/flexidot-bio/flexidot/blob/master/docs/images/all_against_all.png" width="500">

```bash
# All-by-all plot, LCS shading using maximal LCS length
# -y/--lcs_shading_ref: 0 = maximal LCS length
# -x/--lcs_shading
flexidot -i test-seqs.fasta -m 2 -k 10 -y 0 -x
```

## Major features

### Mismatch and ambiguity handling

In diverged or distantly related sequences matches may be interrupted by mismatches, or residues might be represented as ambiguities to refer to frequent variants or mutations. Similarly, relaxed matching is helpful when analyzing error-prone sequences like SMRT reads. Relaxation of the matching conditions thus increases sensitivity, while decreasing specificity.

Firstly, FlexiDot handles **ambiguous residues**, often found in consensus sequences. This allows the comparison of species-specific representations of multigene or repeat families as well as common variants or sequence subfamilies. The ambiguity handling is controlled via`-w/--wobble_conversion`.

Secondly, a defined number of **mismatches** within the window can be allowed with `-S/--substitution_count [number of allowed mismatches (substitutions)]`. This is even less stringent than the ambiguity handling. Please note, that only substitution mutations are allowed but not indels.

Lastly, both mismatch and ambiguity handling can be combined for the analysis.

<img src="https://github.com/flexidot-bio/flexidot/blob/master/docs/images/Fig-Suppl-MismatchesWobbles.png" width="600">

```bash
# Mismatch tolerance -S
#Panel tl
flexidot -i Seq1.fasta Seq4.fasta -m 1 -k 10
#Panel tm
flexidot -i Seq1.fasta Seq4.fasta -m 1 -k 10 -S 1
#Panel tr
flexidot -i Seq1.fasta Seq4.fasta -m 1 -k 10 -S 2

# Wobble -w (tolerate ambiguities)
#Panel bl
flexidot -i Seq1.fasta Seq4.fasta -m 1 -k 10 -w
#Panel bm
flexidot -i Seq1.fasta Seq4.fasta -m 1 -k 10 -w -S 1
#Panel br
flexidot -i Seq1.fasta Seq4.fasta -m 1 -k 10 -w -S 2
```

### Annotation-based shading

Note: See also [**our tutorial**](https://github.com/flexidot-bio/flexidot/blob/master/docs/tutorial_add_annotation.md) on how to integrate annotation shadings with a real-life example.

In FlexiDot self dotplots, annotated sequence regions can be highlighted by **shading** to allow clear assignment of dotplot regions to specific sequence contexts (see Seq2 in self dotplots). The underlying **annotation** information must be provided in general feature format (**gff3**), either as individual file or file list via the `-g/--input_gff_files` option. To customize GFF-based shading, a user-defined configuration file can be provided via the `-G/--gff_color_config option`. Example files are provided in the test-data directory. Please note, that a legend is generated in a separate file.

If you wish to find out more on the gff3 file format used here, Ensembl provides a [good overview](https://www.ensembl.org/info/website/upload/gff3.html).

<img src="https://github.com/flexidot-bio/flexidot/blob/master/docs/images/Selfdotplot_shaded.png" width="500">

```bash
flexidot -i Seq2.fasta -m 0 -k 10 -w -P 5 -g example.gff3 -G gff_color.config
```

### [since FlexiDot_v1.03] Annotation-based shading also available for all-against-all dotplots

Previously only available for self dotplots, we added annotation-based shading to all-against-all dotplots, allowing for many new visualizations. As before, annotation information is provided as general feature file (GFF3). These features are added to the middle diagonal (see our example below).

<img src="https://github.com/flexidot-bio/flexidot/blob/master/docs/images/all_against_all_annotation_based_shading_cool.png" width="700">

Basic command:

```bash
flexidot -i test-seqs.fasta -g example2.gff3 -G gff_color.config -m 2
```

Command plus aesthetics as shown here (+ LCS shading, wordsize 10, change of subplot spacing and line width):

```bash
flexidot -i test-seqs.fasta -g example2.gff3 -G gff_color.config -m 2 -x -k 10 -F 0.06 -A 1.5
```

The test files used here are [provided](https://github.com/flexidot-bio/flexidot/tree/master/test-data):

* [test-seqs.fasta](https://github.com/flexidot-bio/flexidot/blob/master/tests/test-data/test-seqs.fasta)
* [example2.gff3](https://github.com/flexidot-bio/flexidot/blob/master/tests/test-data/example2.gff3)
* [gff_color.config](https://github.com/flexidot-bio/flexidot/blob/master/tests/test-data/gff_color.config)

### Similarity shading

In all-against-all mode, FlexiDot compares each pair from a set of input sequences. To enable the identification of long shared subsequences at a glance, FlexiDot offers similarity shading (switched on/off via option `-x/--lcs_shading`) based on the **LCS length** (longest common subsequence, or longest match if mismatches are considered) in all-against-all comparisons. Longer matches are represented by darker background shading. A separate shading **legend** output file is created written according to mathematical interval notation, where interval boundaries are represented by a pair of numbers. Consequently, the symbols “(” or “)” represent exclusion, whereas “[” or “]” represent inclusion of the respective number.

FlexiDot similarity shading is highly customizable with the following parameters, explained in depth in the documentation:

* Reference for shading (option `-y/--lcs_shading_ref`)
* Number of shading intervals (option `-X/--lcs_shading_num`)
* Shading based on sequence orientation (option `-z/--lcs_shading_ori`)

Shading examples based on sequence orientation (forward, panel A; reverse, panel B; both, panel C) are shown:

![alt text](https://github.com/flexidot-bio/flexidot/blob/master/docs/images/all_against_all_shaded_orientation2.png "FlexiDot shaded dotplots")

```bash
#Panel A - lcs_shading_ori: 0 = forward
flexidot -i test-seqs.fasta -m 2 -k 10 -x -y 0 -z 0
#Panel B - lcs_shading_ori: 1 = reverse
flexidot -i test-seqs.fasta -m 2 -k 10 -x -y 0 -z 1
#Panel C - lcs_shading_ori: 2 = both
flexidot -i test-seqs.fasta -m 2 -k 10 -x -y 0 -z 2
```

### Custom matrix shading

When comparing related sequences, multiple sequence alignments are frequently applied. The resulting pairwise **sequence similarities** can be integrated in the FlexiDot images by providing a **matrix file** via `-u/--user_matrix_file <matrix.txt>`. This allows a shading of the upper right triangle according to the matrix (here orange). With `-U/--user_matrix_print` the matrix values can be printed into the respective fields. Besides, also **text** information can be provided in the matrix, but then shading is suppressed.

In the example, LCS and matrix shading are combined to visualize the relationships between different members of a repeat family.

<img src="https://github.com/flexidot-bio/flexidot/blob/master/docs/images/Beetle_matrix_shading.png" width="750">

```bash
# Beetle TE plot
flexidot -i Beetle.fas -m 2 -k 10 -S 1 -r -x -u custom_matrix.txt -U

# Example with test dataset
flexidot -i test-seqs.fasta -m 2 -k 10 -S 1 -x -u custom_matrix.txt -U
```

# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.2] 2025-11-10

### Changed

- Bump third-party dependency versions in GitHub actions.
- Minor formatting fixes.

### Added

- Add CFF maintenance actions.

## [2.0.0] - 2025-05

This release is a major refactor of the Flexidot codebase that migrates to Python 3 and a modern package structure.

### Added

- Flexidot and its dependencies are now pip installable - uses Hatch and pyproject.toml
- Versioning is now managed dynamically using git tags
- Some basic tests for the core `find_match_pos_diag()` function have been added
- Command line options are now managed with argparse
- Repo includes env yaml to set up conda env for flexidot
- Check that input files exist
- Auto cleanup temp files
- Add action to run pytests
- Add action to format code with Ruff
- Uses logging module to manage status logging (removed time logging)

### Changed

- When comparing a sequence to an identical sequence `find_match_pos_diag()` will recycle kmer counts from the first seq. Saves 33% runtime.
- Several command line options have been renamed or have changes to their expected input formatting. See --help.
- If not using the `--wobble_conversion` option then kmers containing any Ns will be skipped by default.
- If `--wobble_conversion` is set then `--max_n` determines the max percentage of Ns that will be tolerated in a kmer. Default changed to 10% from hard coded 49%.

### Fixed

- Fix deprecation issue with numpy creating an ndarray from ragged nested sequences in `find_match_pos_diag()`. Closes issue #15
- Read files with r instead of rb
- Fix unicode issue referenced in #10

## [1.06] - 2019-04-14

- [new parameter cheat sheet v1.06](https://github.com/molbio-dresden/flexidot/blob/master/documentation/usage_v1.06.pdf)
- [new FlexiDot script v1.06](https://github.com/molbio-dresden/flexidot/blob/master/code/flexidot_v1.06.py)

### Fixed

- Corrected a few bugs, including a bug introduced in version 1.05, affecting dotplots with substitutions allowed
- Reverted (for now) to the pattern matching algorithm of version 1.04

## [1.05] - 2018-12-14

- [new parameter cheat sheet v1.05](https://github.com/molbio-dresden/flexidot/blob/master/documentation/usage_v1.05.pdf)
- [new FlexiDot script v1.05](https://github.com/molbio-dresden/flexidot/blob/master/code/flexidot_v1.05.py)

### Changed

- Modified word match recognition, speeding up FlexiDot's runtime
- The default wordsize has been changed from 7 to 10 in order to prevent people from running FlexiDot with small word sizes on large datasets, as this presumably takes a very long time

### Added

- New `-O, --only_vs_first_seq` option to limit the output of pairwise dotplots. Instead of printing all possible pairwise combinations from a multi-fasta-sequence, only the pairwise comparisons against the first sequence are generated, if switched on (`-O y`). We use this feature to compare a new/unknown sequence against a batch of references.

### Fixed

- Fixed a few bugs with the dotplot shading legends

## [1.04] - 2018-06-29

- [new parameter cheat sheet v1.04](https://github.com/molbio-dresden/flexidot/blob/master/documentation/usage_v1.04.pdf)
- [new FlexiDot script v1.04](https://github.com/molbio-dresden/flexidot/blob/master/code/flexidot_v1.04.py)

### Added

- New graphic formatting options for all-against-all dotplots
- With `-M/--mirror` it is now possible to mirror the middle diagonal

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/all_against_all_Flexi1.04_Para_Mirror.png" width="700">

Basic commands:

```text
python flexidot.py -i test-seqs.fas -p 2 -M n
python flexidot.py -i test-seqs.fas -p 2 -M y
```

Command plus aesthetics as shown here (as described in version update 1.03):

```text
python flexidot.py -i test-seqs.fas -p 2 -M n -g example2.gff3 -G gff_color.config -x y -k 10 -F 0.06 -A 1.5
python flexidot.py -i test-seqs.fas -p 2 -M y -g example2.gff3 -G gff_color.config -x y -k 10 -F 0.06 -A 1.5
```

- The `-R/--representation` parameter allows partial dotplotting, either printing the complete `-R 0`, the top `-R 1`, or the bottom dotplot `-R 2`

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/all_against_all_Flexi1.04_Para_Representation.png" width="900">

Basic commands:

```text
python flexidot.py -i test-seqs.fas -p 2 -R 0
python flexidot.py -i test-seqs.fas -p 2 -R 1
python flexidot.py -i test-seqs.fas -p 2 -R 2
```

Command plus aesthetics as shown here (as described in version update 1.03):

```text
python flexidot.py -i test-seqs.fas -p 2 -R 0 -g example2.gff3 -G gff_color.config -x y -k 10 -F 0.06 -A 1.5
python flexidot.py -i test-seqs.fas -p 2 -R 1 -g example2.gff3 -G gff_color.config -x y -k 10 -F 0.06 -A 1.5
python flexidot.py -i test-seqs.fas -p 2 -R 2 -g example2.gff3 -G gff_color.config -x y -k 10 -F 0.06 -A 1.5
```

### Fixed

- Fixed a distortion issue in `-p/--plotting_mode 0` (self dotplots)

## [1.03] - 2018-06-17

- [new parameter cheat sheet v1.03](https://github.com/molbio-dresden/flexidot/blob/master/documentation/usage_v1.03.pdf)
- [new FlexiDot script v1.03](https://github.com/molbio-dresden/flexidot/blob/master/code/flexidot_v1.03.py)

### Added

- Annotation-based shading now available for all-against-all dotplots. Previously only available for self dotplots, we added annotation-based shading to all-against-all dotplots, allowing for many new visualizations. As before, annotation information is provided as general feature file (GFF3). These features are added to the middle diagonal (see our example below).

<img src="https://github.com/molbio-dresden/flexidot/blob/master/images/all_against_all_annotation_based_shading_cool.png" width="700">

Basic command:

```text
python flexidot.py -i test-seqs.fas -g example2.gff3 -G gff_color.config -p 2
```

Command plus aesthetics as shown here (+ LCS shading, wordsize 10, change of subplot spacing and line width):

```text
python flexidot.py -i test-seqs.fas -g example2.gff3 -G gff_color.config -p 2 -x y -k 10 -F 0.06 -A 1.5
```

The test files used here are provided:

- [test-seqs.fas](https://github.com/molbio-dresden/flexidot/blob/master/test-data/test-seqs.fas)
- [example2.gff3](https://github.com/molbio-dresden/flexidot/blob/master/test-data/example2.gff3)
- [gff_color.config](https://github.com/molbio-dresden/flexidot/blob/master/test-data/gff_color.config)

## [1.02] - 2018-05-09

- [new parameter cheat sheet v1.02](https://github.com/molbio-dresden/flexidot/blob/master/documentation/usage_v1.02.pdf)
- [new FlexiDot script v1.02](https://github.com/molbio-dresden/flexidot/blob/master/code/flexidot_v1.02.py)

### Changed

- Changed handling of `-T` parameter: The character count of the sequence titles has been limited to `20` by default. This limit can be changed with `-T`. If an `E` (end) is added to the limit, the last characters are chosen instead of the first.

```text
-T 20  (the first 20 characters)
-T 20E (the last 20 characters)
```

## [1.01] - 2018-04-21

- [parameter cheat sheet v1.01](https://github.com/molbio-dresden/flexidot/blob/master/documentation/usage_v1.01.pdf)
- [FlexiDot script v1.01](https://github.com/molbio-dresden/flexidot/blob/master/code/flexidot_v1.01.py)

### Fixed

- Minor bugfixing

## [1.00] - 2018-03-21

### Added

- First FlexiDot release

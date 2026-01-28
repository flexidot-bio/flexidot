# Plotting Pre-calculated Alignments

FlexiDot now supports plotting pre-calculated alignments from external alignment tools like blastn, Nucmer, and Minimap2. This feature allows you to visualize alignments that have been generated using more sensitive or specialized alignment algorithms, rather than relying solely on k-mer matching.

## Supported Alignment Formats

FlexiDot supports two popular alignment output formats:

### BLAST6 Format (Tabular)

BLAST6 is the tabular output format from BLAST (output format 6). It contains 12 tab-separated columns:

| Column | Description |
|--------|-------------|
| 1 | Query sequence ID |
| 2 | Subject sequence ID |
| 3 | Percent identity |
| 4 | Alignment length |
| 5 | Number of mismatches |
| 6 | Number of gap openings |
| 7 | Query start |
| 8 | Query end |
| 9 | Subject start |
| 10 | Subject end |
| 11 | E-value |
| 12 | Bit score |

Example BLAST6 file:

```text
seq1	seq2	95.5	100	4	1	1	100	1	100	1e-50	180
seq1	seq2	90.0	50	5	0	150	199	200	249	1e-20	90
```

### PAF Format (Pairwise mApping Format)

PAF is the output format used by Minimap2 and other modern aligners. It contains at least 12 tab-separated columns:

| Column | Description |
|--------|-------------|
| 1 | Query sequence name |
| 2 | Query sequence length |
| 3 | Query start (0-based) |
| 4 | Query end (0-based, open) |
| 5 | Strand (+/-) |
| 6 | Target sequence name |
| 7 | Target sequence length |
| 8 | Target start (0-based) |
| 9 | Target end (0-based, open) |
| 10 | Number of matching bases |
| 11 | Alignment block length |
| 12 | Mapping quality |

Example PAF file:

```text
seq1	1000	0	100	+	seq2	1200	0	100	95	100	60
seq1	1000	149	199	+	seq2	1200	199	249	45	50	30
```

## Command-Line Usage

### Basic Usage

To plot pre-calculated alignments, use the `-a` or `--alignment_file` option:

```bash
# Plot alignments from a BLAST6 file
flexidot -i sequences.fasta -a alignments.blast6 -m 1

# Plot alignments from a PAF file
flexidot -i sequences.fasta -a alignments.paf -m 2
```

### Specifying Alignment Format

FlexiDot auto-detects the alignment format from the file extension:
- `.blast6`, `.b6`, `.blastn`, `.blast`, `.m8` → BLAST6 format
- `.paf` → PAF format

If your file has a different extension, specify the format explicitly:

```bash
flexidot -i sequences.fasta -a alignments.txt --alignment_format blast6 -m 1
```

### Filtering Alignments

You can filter alignments by minimum percent identity or minimum length:

```bash
# Only plot alignments with ≥95% identity
flexidot -i sequences.fasta -a alignments.paf --min_identity 95 -m 1

# Only plot alignments ≥100 bp long
flexidot -i sequences.fasta -a alignments.paf --min_length 100 -m 1

# Combine filters
flexidot -i sequences.fasta -a alignments.paf --min_identity 90 --min_length 50 -m 2
```

## Generating Alignment Files

### Using BLASTN

```bash
# Create blast database
mkdir db
makeblastdb -in sequences.fasta -dbtype nucl -out db/sequences_db -parse_seqids

# Run BLASTN with output format 6
blastn -query sequences.fasta -db db/sequences_db -outfmt 6 -out alignments.blast6 \
-word_size 4 -evalue 1e-3 -perc_identity 60.0 -max_target_seqs 10000 -num_threads 8
```

### Using Minimap2

```bash
# Run minimap2 for nucleotide sequences (5% divergence)
minimap2 -x asm5 -t 8 sequences.fasta sequences.fasta > alignments.paf

# For more sensitive alignments (20% divergence)
minimap2 -x asm20 -t 8 sequences.fasta sequences.fasta > alignments.paf
```

### Using Nucmer

```bash
# Self-alignment with nucmer (use --nosimplify for repeats in self alignments)
nucmer --maxmatch --nosimplify --minmatch 15 --mincluster 20 --diagfactor 0.3 \
--prefix self_align sequences.fasta sequences.fasta

# Convert directly using paftools (if installed with minimap2)
paftools.js delta2paf self_align.delta > self_align.paf
```

## Example Workflow

Here's a complete workflow comparing k-mer matching with pre-calculated alignments:

### 1. Test case data

```bash
SEQ="tests/test-data/sSaTar_example/sSaTar.fas"
ANNOTATION="tests/test-data/sSaTar_example/sSaTar.gff3"
COLOURS="tests/test-data/sSaTar_example/sSaTar.config"
COLORS=$COLOURS
```


### 2. Standard K-mer Matching

```bash
# Use FlexiDot's built-in k-mer matching
flexidot -i $SEQ -m 2 -k 15 -o kmer_dotplot --gff $ANNOTATION --gff_color_config $COLOURS
```

### 3. Using BLAST Alignments

```bash
# Generate BLAST alignments
blastn -query $SEQ -subject $SEQ -outfmt 6 -word_size 4 -perc_identity 60.0 -max_target_seqs 10000 -evalue 0.001 -out alignments.blast6

# Plot alignments
flexidot -i $SEQ -m 2 -a alignments.blast6 -m 2 -o blast_dotplot --gff $ANNOTATION --gff_color_config $COLOURS --min_identity 80 --min_length 20
```

### 3. Using PAF alignments from Nucmer

```bash
# Self-alignment with nucmer (use --nosimplify for repeats in self alignments)
nucmer --maxmatch --nosimplify --minmatch 15 --mincluster 20 --diagfactor 0.3 \
--prefix self_align $SEQ $SEQ

# Convert directly using paftools (if installed with minimap2)
paftools.js delta2paf self_align.delta > self_align.paf

# Plot alignments
flexidot -i $SEQ -a self_align.paf -m 2 -o nucmer_dotplot --gff $ANNOTATION --gff_color_config $COLOURS
```

### 4. Using Minimap2 Alignments

```bash
# Generate minimap2 alignments
# Note: Minimap2 is not particularly well suited to detecting small secondary alignments in small sequences.
# Try tinkering with settings: -k 10 -N 1000 -p 0.05 -r 2k, or use Nucmer.
minimap2 -x asm20 $SEQ $SEQ > alignments.paf

# Plot alignments
flexidot -i $SEQ -a alignments.paf -m 2 -o minimap_dotplot --gff $ANNOTATION --gff_color_config $COLOURS
```

## Tips and Best Practices

1. **Redundant Alignment Filtering**: FlexiDot automatically filters redundant alignments where the same sequence pair appears in both directions (e.g., SeqA vs SeqB and SeqB vs SeqA). Only one copy is kept.

2. **Sequence Names**: Ensure the sequence names in your FASTA file match exactly the names in your alignment file. FlexiDot uses these names to associate alignments with the correct sequences.

3. **Self-Alignments**: Self-alignments (sequence aligned to itself) are preserved and can be useful for identifying repeats within sequences.

4. **Strand Information**:
   - In BLAST6 format, strand is determined by the subject coordinates (start > end indicates reverse strand).
   - In PAF format, strand is explicitly provided in column 5 (+/-).

5. **Performance**: Using pre-calculated alignments can be significantly faster than k-mer matching for large datasets, especially when alignments have already been computed for other purposes.

## Comparison: K-mer Matching vs Pre-calculated Alignments

| Aspect | K-mer Matching | Pre-calculated Alignments |
|--------|---------------|---------------------------|
| Speed | Fast for small datasets | Very fast (alignments already computed) |
| Sensitivity | Limited by k-mer size | Depends on alignment tool |
| Gap handling | No gap tolerance | Handles gaps (depending on aligner) |
| Mismatch tolerance | Limited (with `-S` option) | Full flexibility |
| Setup | Built-in | Requires external tool |
| Use case | Quick visualization | Sensitive comparisons |

## Troubleshooting

### Common Issues

1. **No alignments plotted**:
   - Check that sequence names in the alignment file match the FASTA headers
   - Verify the alignment file format is correct
   - Try relaxing the `--min_identity` or `--min_length` filters

2. **Format detection fails**:
   - Explicitly specify the format with `--alignment_format`

3. **Some sequences missing from plot**:
   - Ensure all sequences in your FASTA file have at least one alignment in the alignment file

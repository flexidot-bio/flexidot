"""Functions for parsing and processing pre-calculated alignments.

This module provides functions to parse alignment files in BLAST6 (tabular) and
PAF (Pairwise mApping Format) formats, and convert them to coordinates suitable
for plotting in FlexiDot.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

import numpy as np


def parse_blast6(
    filepath: Union[str, Path],
    min_identity: float = 0.0,
    min_length: int = 0,
) -> List[Dict]:
    """
    Parse a BLAST6 (tabular) format alignment file.

    BLAST6 format columns (tab-separated):
    0: qseqid   - Query sequence ID
    1: sseqid   - Subject sequence ID
    2: pident   - Percentage identity
    3: length   - Alignment length
    4: mismatch - Number of mismatches
    5: gapopen  - Number of gap openings
    6: qstart   - Start of alignment in query
    7: qend     - End of alignment in query
    8: sstart   - Start of alignment in subject
    9: send     - End of alignment in subject
    10: evalue  - E-value
    11: bitscore - Bit score

    Parameters
    ----------
    filepath : str or Path
        Path to the BLAST6 format file.
    min_identity : float, optional
        Minimum percent identity to include alignment (0-100). Default is 0.0.
    min_length : int, optional
        Minimum alignment length to include. Default is 0.

    Returns
    -------
    list of dict
        List of alignment dictionaries, each containing:
        - query_id: str
        - subject_id: str
        - query_start: int
        - query_end: int
        - subject_start: int
        - subject_end: int
        - identity: float
        - length: int
        - strand: str ('+' or '-')

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    ValueError
        If the file format is invalid or min_identity is out of range.

    Examples
    --------
    >>> alignments = parse_blast6('alignments.blast6')
    >>> for aln in alignments:
    ...     print(f"{aln['query_id']} vs {aln['subject_id']}")
    """
    # Validate min_identity
    if not 0 <= min_identity <= 100:
        raise ValueError(f'min_identity must be between 0 and 100, got {min_identity}')

    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f'Alignment file not found: {filepath}')

    alignments = []

    with open(filepath) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < 12:
                logging.warning(
                    f'Skipping malformed line {line_num}: expected 12 columns, got {len(fields)}'
                )
                continue

            try:
                query_id = fields[0]
                subject_id = fields[1]
                identity = float(fields[2])
                length = int(fields[3])
                query_start = int(fields[6])
                query_end = int(fields[7])
                subject_start = int(fields[8])
                subject_end = int(fields[9])
            except (ValueError, IndexError) as e:
                logging.warning(f'Error parsing line {line_num}: {e}')
                continue

            # Apply filters
            if identity < min_identity or length < min_length:
                continue

            # Determine strand based on subject coordinates
            # In BLAST6, if sstart > send, the alignment is on the minus strand
            if subject_start <= subject_end:
                strand = '+'
            else:
                strand = '-'
                # Swap to ensure start < end for consistent processing
                subject_start, subject_end = subject_end, subject_start

            alignments.append(
                {
                    'query_id': query_id,
                    'subject_id': subject_id,
                    'query_start': query_start,
                    'query_end': query_end,
                    'subject_start': subject_start,
                    'subject_end': subject_end,
                    'identity': identity,
                    'length': length,
                    'strand': strand,
                }
            )

    if len(alignments) == 0:
        logging.warning(
            f'No alignments found in BLAST6 file: {filepath}. '
            'Plot will be generated without alignment overlays.'
        )
    else:
        logging.info(f'Parsed {len(alignments)} alignments from BLAST6 file: {filepath}')
    return alignments


def parse_paf(
    filepath: Union[str, Path],
    min_identity: float = 0.0,
    min_length: int = 0,
) -> List[Dict]:
    """
    Parse a PAF (Pairwise mApping Format) alignment file.

    PAF format columns (tab-separated):
    0: query_name    - Query sequence name
    1: query_length  - Query sequence length
    2: query_start   - Query start coordinate (0-based)
    3: query_end     - Query end coordinate (0-based, open)
    4: strand        - '+' or '-'
    5: target_name   - Target sequence name
    6: target_length - Target sequence length
    7: target_start  - Target start coordinate (0-based)
    8: target_end    - Target end coordinate (0-based, open)
    9: matches       - Number of matching bases
    10: block_length - Total number of bases in alignment
    11: mapq         - Mapping quality (0-255, 255 = unavailable)

    Parameters
    ----------
    filepath : str or Path
        Path to the PAF format file.
    min_identity : float, optional
        Minimum percent identity to include alignment (0-100). Default is 0.0.
    min_length : int, optional
        Minimum alignment length to include. Default is 0.

    Returns
    -------
    list of dict
        List of alignment dictionaries, each containing:
        - query_id: str
        - subject_id: str
        - query_start: int (1-based)
        - query_end: int (1-based, inclusive)
        - subject_start: int (1-based)
        - subject_end: int (1-based, inclusive)
        - identity: float
        - length: int
        - strand: str ('+' or '-')

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    ValueError
        If the file format is invalid or min_identity is out of range.

    Examples
    --------
    >>> alignments = parse_paf('alignments.paf')
    >>> for aln in alignments:
    ...     print(f"{aln['query_id']} vs {aln['subject_id']}: {aln['strand']}")
    """
    # Validate min_identity
    if not 0 <= min_identity <= 100:
        raise ValueError(f'min_identity must be between 0 and 100, got {min_identity}')

    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f'Alignment file not found: {filepath}')

    alignments = []

    with open(filepath) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < 12:
                logging.warning(
                    f'Skipping malformed line {line_num}: expected at least 12 columns, got {len(fields)}'
                )
                continue

            try:
                query_id = fields[0]
                query_start = int(fields[2])  # 0-based
                query_end = int(fields[3])  # 0-based, open
                strand = fields[4]
                subject_id = fields[5]
                subject_start = int(fields[7])  # 0-based
                subject_end = int(fields[8])  # 0-based, open
                matches = int(fields[9])
                block_length = int(fields[10])
            except (ValueError, IndexError) as e:
                logging.warning(f'Error parsing line {line_num}: {e}')
                continue

            # Calculate identity
            if block_length > 0:
                identity = (matches / block_length) * 100
            else:
                identity = 0.0

            length = block_length

            # Apply filters
            if identity < min_identity or length < min_length:
                continue

            # Convert PAF 0-based half-open coordinates to 1-based inclusive
            # PAF: [start, end) 0-based -> 1-based inclusive: start+1, end
            alignments.append(
                {
                    'query_id': query_id,
                    'subject_id': subject_id,
                    'query_start': query_start + 1,  # Convert to 1-based
                    'query_end': query_end,  # Already correct (half-open to inclusive)
                    'subject_start': subject_start + 1,  # Convert to 1-based
                    'subject_end': subject_end,  # Already correct
                    'identity': identity,
                    'length': length,
                    'strand': strand,
                }
            )

    if len(alignments) == 0:
        logging.warning(
            f'No alignments found in PAF file: {filepath}. '
            'Plot will be generated without alignment overlays.'
        )
    else:
        logging.info(f'Parsed {len(alignments)} alignments from PAF file: {filepath}')
    return alignments


def filter_redundant_alignments(
    alignments: List[Dict],
) -> List[Dict]:
    """
    Filter redundant alignments where sequences are duplicated between query and target.

    When the same sequences appear in both query and target sets, alignments may be
    present in both directions (SeqA vs SeqB and SeqB vs SeqA). This function keeps
    only one copy of each alignment pair, preferring the alignment where the query
    ID is alphabetically before the subject ID.

    Parameters
    ----------
    alignments : list of dict
        List of alignment dictionaries from parse_blast6 or parse_paf.

    Returns
    -------
    list of dict
        Filtered list of alignments with redundant pairs removed.

    Examples
    --------
    >>> alignments = [
    ...     {'query_id': 'seqA', 'subject_id': 'seqB', 'query_start': 1, 'query_end': 100,
    ...      'subject_start': 1, 'subject_end': 100, 'strand': '+', 'identity': 95, 'length': 100},
    ...     {'query_id': 'seqB', 'subject_id': 'seqA', 'query_start': 1, 'query_end': 100,
    ...      'subject_start': 1, 'subject_end': 100, 'strand': '+', 'identity': 95, 'length': 100},
    ... ]
    >>> filtered = filter_redundant_alignments(alignments)
    >>> len(filtered)
    1
    """
    seen_pairs: Set[Tuple[str, str, int, int, int, int, str]] = set()
    filtered = []

    for aln in alignments:
        query_id = aln['query_id']
        subject_id = aln['subject_id']

        # Skip self-alignments where query and subject are the same sequence
        if query_id == subject_id:
            # Keep self-alignments (they're useful for self dotplots)
            filtered.append(aln)
            continue

        # Create a canonical key for this alignment pair
        # Sort IDs alphabetically to create consistent key regardless of direction
        if query_id < subject_id:
            pair_key = (
                query_id,
                subject_id,
                aln['query_start'],
                aln['query_end'],
                aln['subject_start'],
                aln['subject_end'],
                aln['strand'],
            )
        else:
            # Reverse the pair and swap coordinates
            pair_key = (
                subject_id,
                query_id,
                aln['subject_start'],
                aln['subject_end'],
                aln['query_start'],
                aln['query_end'],
                aln['strand'],
            )

        if pair_key not in seen_pairs:
            seen_pairs.add(pair_key)
            filtered.append(aln)

    removed_count = len(alignments) - len(filtered)
    if removed_count > 0:
        logging.info(f'Filtered {removed_count} redundant alignments')

    return filtered


def alignments_to_coordinates(
    alignments: List[Dict],
    query_id: str,
    subject_id: str,
    subject_length: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert alignments for a specific sequence pair to plotting coordinates.

    This function extracts alignments between the specified query and subject
    sequences and converts them to the coordinate format expected by FlexiDot's
    plotting functions.

    Parameters
    ----------
    alignments : list of dict
        List of alignment dictionaries from parse_blast6 or parse_paf.
    query_id : str
        ID of the query sequence (plotted on x-axis).
    subject_id : str
        ID of the subject sequence (plotted on y-axis).
    subject_length : int
        Length of the subject sequence (needed for reverse complement coordinates).

    Returns
    -------
    tuple of numpy.ndarray
        Four arrays containing:
        - x1: x coordinates for forward strand alignments
        - y1: y coordinates for forward strand alignments
        - x2: x coordinates for reverse strand alignments
        - y2: y coordinates for reverse strand alignments

        Each array contains sub-arrays where each sub-array represents
        the start and end coordinates of a single alignment line.

    Examples
    --------
    >>> alignments = [
    ...     {'query_id': 'seq1', 'subject_id': 'seq2', 'query_start': 10,
    ...      'query_end': 50, 'subject_start': 20, 'subject_end': 60,
    ...      'strand': '+', 'identity': 95, 'length': 41},
    ... ]
    >>> x1, y1, x2, y2 = alignments_to_coordinates(alignments, 'seq1', 'seq2', 100)
    >>> len(x1)  # One forward alignment
    1
    """
    x1_list = []  # x coords for forward
    y1_list = []  # y coords for forward
    x2_list = []  # x coords for reverse
    y2_list = []  # y coords for reverse

    for aln in alignments:
        # Check if this alignment involves the requested sequences
        # Handle both directions: query->subject and subject->query
        if aln['query_id'] == query_id and aln['subject_id'] == subject_id:
            q_start = aln['query_start']
            q_end = aln['query_end']
            s_start = aln['subject_start']
            s_end = aln['subject_end']
            strand = aln['strand']
        elif aln['query_id'] == subject_id and aln['subject_id'] == query_id:
            # Swap roles: what was query becomes subject and vice versa
            q_start = aln['subject_start']
            q_end = aln['subject_end']
            s_start = aln['query_start']
            s_end = aln['query_end']
            strand = aln['strand']
        else:
            continue

        # Create coordinate arrays for the alignment line
        # FlexiDot plots x (query) on horizontal, y (subject) on vertical
        if strand == '+':
            x1_list.append(np.array([q_start, q_end]))
            y1_list.append(np.array([s_start, s_end]))
        else:
            # For reverse strand, we need to flip the y coordinates
            # The line goes from (q_start, s_end) to (q_end, s_start)
            x2_list.append(np.array([q_start, q_end]))
            y2_list.append(np.array([s_end, s_start]))

    return (
        np.array(x1_list, dtype=object),
        np.array(y1_list, dtype=object),
        np.array(x2_list, dtype=object),
        np.array(y2_list, dtype=object),
    )


def load_alignments(
    filepath: Union[str, Path],
    file_format: Optional[str] = None,
    min_identity: float = 0.0,
    min_length: int = 0,
    filter_redundant: bool = True,
) -> List[Dict]:
    """
    Load alignments from a file, auto-detecting format if not specified.

    Parameters
    ----------
    filepath : str or Path
        Path to the alignment file.
    file_format : str, optional
        Format of the alignment file: 'blast6' or 'paf'.
        If None, format is auto-detected from file extension.
    min_identity : float, optional
        Minimum percent identity to include alignment (0-100). Default is 0.0.
    min_length : int, optional
        Minimum alignment length to include. Default is 0.
    filter_redundant : bool, optional
        If True, filter redundant alignments. Default is True.

    Returns
    -------
    list of dict
        List of alignment dictionaries.

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    ValueError
        If the file format cannot be determined.

    Examples
    --------
    >>> alignments = load_alignments('alignments.paf')
    >>> alignments = load_alignments('alignments.txt', file_format='blast6')
    """
    filepath = Path(filepath)

    # Auto-detect format from extension if not specified
    if file_format is None:
        suffix = filepath.suffix.lower()
        if suffix in ['.paf']:
            file_format = 'paf'
        elif suffix in ['.blast6', '.b6', '.blastn', '.blast', '.m8']:
            file_format = 'blast6'
        else:
            # Try to detect from content
            file_format = _detect_alignment_format(filepath)

    file_format = file_format.lower()

    if file_format == 'blast6':
        alignments = parse_blast6(
            filepath, min_identity=min_identity, min_length=min_length
        )
    elif file_format == 'paf':
        alignments = parse_paf(
            filepath, min_identity=min_identity, min_length=min_length
        )
    else:
        raise ValueError(
            f"Unknown alignment format: {file_format}. Use 'blast6' or 'paf'."
        )

    if filter_redundant:
        alignments = filter_redundant_alignments(alignments)

    return alignments


def _detect_alignment_format(filepath: Path) -> str:
    """
    Attempt to detect alignment file format from content.

    Parameters
    ----------
    filepath : Path
        Path to the alignment file.

    Returns
    -------
    str
        Detected format: 'blast6' or 'paf'.

    Raises
    ------
    ValueError
        If the format cannot be determined.
    """
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) >= 12:
                # Check if field 4 (index 4) is '+' or '-' (PAF strand)
                if fields[4] in ['+', '-']:
                    return 'paf'
                # Check if field 2 looks like percent identity (BLAST6)
                try:
                    pident = float(fields[2])
                    if 0 <= pident <= 100:
                        return 'blast6'
                except ValueError:
                    pass

            break

    raise ValueError(
        f'Could not auto-detect alignment format for {filepath}. '
        'Please specify format using --alignment_format option.'
    )


def get_sequence_ids_from_alignments(alignments: List[Dict]) -> Set[str]:
    """
    Extract all unique sequence IDs from alignments.

    Parameters
    ----------
    alignments : list of dict
        List of alignment dictionaries.

    Returns
    -------
    set of str
        Set of unique sequence IDs found in alignments.

    Examples
    --------
    >>> alignments = [
    ...     {'query_id': 'seq1', 'subject_id': 'seq2', ...},
    ...     {'query_id': 'seq1', 'subject_id': 'seq3', ...},
    ... ]
    >>> ids = get_sequence_ids_from_alignments(alignments)
    >>> sorted(ids)
    ['seq1', 'seq2', 'seq3']
    """
    ids = set()
    for aln in alignments:
        ids.add(aln['query_id'])
        ids.add(aln['subject_id'])
    return ids

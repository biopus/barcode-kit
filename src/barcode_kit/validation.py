from __future__ import annotations

from Bio.Seq import Seq

from barcode_kit.models import SequenceQuality


def sequence_quality(sequence: Seq | str) -> SequenceQuality:
    seq = str(sequence).upper().replace(" ", "").replace("\n", "")
    length = len(seq)
    canonical_count = sum(seq.count(base) for base in "ACGT")
    gc_content = ((seq.count("G") + seq.count("C")) / canonical_count) if canonical_count else 0.0
    ambiguous_content = sum(1 for base in seq if base not in "ACGT") / length if length else 0.0
    return SequenceQuality(
        length=length,
        gc_content=gc_content,
        ambiguous_content=ambiguous_content,
    )

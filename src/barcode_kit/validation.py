from __future__ import annotations

from Bio.Seq import Seq

from barcode_kit.models import Marker, SequenceQuality


STOP_CODONS = {"TAA", "TAG", "TGA"}


def sequence_quality(sequence: Seq | str, marker: Marker) -> SequenceQuality:
    seq = str(sequence).upper().replace(" ", "").replace("\n", "")
    length = len(seq)
    canonical_count = sum(seq.count(base) for base in "ACGT")
    gc_content = ((seq.count("G") + seq.count("C")) / canonical_count) if canonical_count else 0.0
    n_content = seq.count("N") / length if length else 0.0
    has_stop_codon: bool | None = None
    has_frameshift: bool | None = None
    if marker.is_coding:
        has_frameshift = length % 3 != 0
        codons = {seq[index : index + 3] for index in range(0, max(length - 2, 0), 3)}
        has_stop_codon = bool(codons & STOP_CODONS)
    return SequenceQuality(
        length=length,
        gc_content=gc_content,
        n_content=n_content,
        has_stop_codon=has_stop_codon,
        has_frameshift=has_frameshift,
    )

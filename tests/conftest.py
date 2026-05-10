from __future__ import annotations

import io

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


def genbank_text(
    accession: str = "PP476489",
    version: int = 4,
    organism: str = "Iris japonica",
    taxon_id: int = 12345,
    gene: str = "rbcL",
    sequence: str | None = None,
) -> str:
    sequence = sequence or ("ATG" + "GCT" * 199)
    record = SeqRecord(
        Seq(sequence),
        id=f"{accession}.{version}",
        name=accession,
        description=f"{organism} voucher test {gene}",
    )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["accessions"] = [accession]
    record.annotations["sequence_version"] = version
    record.annotations["organism"] = organism
    record.features = [
        SeqFeature(
            FeatureLocation(0, len(sequence), strand=1),
            type="source",
            qualifiers={"organism": [organism], "db_xref": [f"taxon:{taxon_id}"]},
        ),
        SeqFeature(
            FeatureLocation(0, len(sequence), strand=1),
            type="gene",
            qualifiers={"gene": [gene]},
        ),
    ]
    handle = io.StringIO()
    SeqIO.write(record, handle, "genbank")
    return handle.getvalue()


@pytest.fixture
def genbank_text_factory():
    return genbank_text


@pytest.fixture(name="genbank_text")
def genbank_text_fixture():
    return genbank_text

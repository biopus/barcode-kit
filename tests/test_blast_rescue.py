from __future__ import annotations

from Bio.Seq import Seq

from barcode_kit import blast as blast_module
from barcode_kit.blast import BlastRecord, BlastRunner
from barcode_kit.config import BlastRescueConfig, BlastRescueMarkerConfig
from barcode_kit.models import Marker


def test_blast_module_exports_only_runner_contract():
    assert blast_module.__all__ == [
        "BlastRecord",
        "BlastRescueResult",
        "BlastRunner",
    ]


def test_blast_rescue_accepts_full_subject_coverage_for_its():
    hit = blast_module.BlastHit(
        query_id="failed",
        subject_id="seed1",
        identity=90.0,
        alignment_length=576,
        qstart=101,
        qend=676,
        sstart=5,
        send=580,
        evalue=1e-80,
        bitscore=500.0,
    )

    result = _rescue_with_hits([hit], Marker.ITS, {"seed1": 600})

    assert result.fallback_reason is None
    assert result.sequence == Seq("A" * 576)
    assert result.metadata["blast_seed_accession"] == "seed1"
    assert result.metadata["blast_subject_coverage"] == 0.96
    assert result.metadata["blast_query_start"] == 101
    assert result.metadata["blast_query_end"] == 676
    assert result.metadata["blast_strand"] == "+"


def test_blast_rescue_rejects_low_subject_coverage():
    hit = blast_module.BlastHit(
        query_id="failed",
        subject_id="seed1",
        identity=95.0,
        alignment_length=401,
        qstart=101,
        qend=501,
        sstart=50,
        send=450,
        evalue=1e-40,
        bitscore=300.0,
    )

    result = _rescue_with_hits([hit], Marker.ITS, {"seed1": 600})

    assert result.sequence is None
    assert result.fallback_reason == "blast_subject_coverage_low"


def test_blast_rescue_rejects_low_identity():
    hit = blast_module.BlastHit(
        query_id="failed",
        subject_id="seed1",
        identity=79.0,
        alignment_length=576,
        qstart=101,
        qend=676,
        sstart=5,
        send=580,
        evalue=1e-40,
        bitscore=300.0,
    )

    result = _rescue_with_hits([hit], Marker.ITS, {"seed1": 600})

    assert result.sequence is None
    assert result.fallback_reason == "blast_identity_low"


def test_blast_rescue_rejects_seed_endpoint_miss():
    hit = blast_module.BlastHit(
        query_id="failed",
        subject_id="seed1",
        identity=95.0,
        alignment_length=521,
        qstart=101,
        qend=621,
        sstart=70,
        send=590,
        evalue=1e-40,
        bitscore=300.0,
    )

    result = _rescue_with_hits([hit], Marker.ITS, {"seed1": 600})

    assert result.sequence is None
    assert result.fallback_reason == "blast_seed_endpoint_miss"


def test_blast_rescue_rejects_implausible_query_span_length():
    hit = blast_module.BlastHit(
        query_id="failed",
        subject_id="seed1",
        identity=95.0,
        alignment_length=576,
        qstart=101,
        qend=1050,
        sstart=5,
        send=580,
        evalue=1e-40,
        bitscore=300.0,
    )

    result = _rescue_with_hits([hit], Marker.ITS, {"seed1": 600})

    assert result.sequence is None
    assert result.fallback_reason == "blast_query_length_ratio"


def test_blast_rescue_rejects_ambiguous_distant_high_scoring_hit():
    best = blast_module.BlastHit(
        query_id="failed",
        subject_id="seed1",
        identity=95.0,
        alignment_length=576,
        qstart=101,
        qend=676,
        sstart=5,
        send=580,
        evalue=1e-80,
        bitscore=500.0,
    )
    second = blast_module.BlastHit(
        query_id="failed",
        subject_id="seed2",
        identity=94.0,
        alignment_length=576,
        qstart=1201,
        qend=1776,
        sstart=5,
        send=580,
        evalue=1e-75,
        bitscore=480.0,
    )

    result = _rescue_with_hits([best, second], Marker.ITS, {"seed1": 600, "seed2": 600})

    assert result.sequence is None
    assert result.fallback_reason == "ambiguous_blast_hit"


def test_blast_rescue_reverse_complements_minus_orientation_hit():
    sequence = Seq("AAAACCCCGGGGTTTT")
    hit = blast_module.BlastHit(
        query_id="failed",
        subject_id="seed1",
        identity=99.0,
        alignment_length=8,
        qstart=12,
        qend=5,
        sstart=1,
        send=8,
        evalue=1e-20,
        bitscore=100.0,
    )

    result = blast_module._rescue_record(
        BlastRecord(accession_version="failed", sequence=sequence),
        [hit],
        Marker.ITS,
        {"seed1": 8},
        BlastRescueConfig(),
    )

    assert result.sequence == Seq("CCCCGGGG").reverse_complement()


def test_blast_rescue_uses_configured_marker_thresholds():
    hit = blast_module.BlastHit(
        query_id="failed",
        subject_id="seed1",
        identity=90.0,
        alignment_length=576,
        qstart=101,
        qend=676,
        sstart=5,
        send=580,
        evalue=1e-80,
        bitscore=500.0,
    )
    config = BlastRescueConfig(
        its=BlastRescueMarkerConfig(
            min_subject_coverage=0.85,
            min_identity=0.92,
            min_query_length_ratio=0.75,
            max_query_length_ratio=1.30,
        )
    )

    result = _rescue_with_hits([hit], Marker.ITS, {"seed1": 600}, config)

    assert result.sequence is None
    assert result.fallback_reason == "blast_identity_low"


def test_subprocess_blast_runner_uses_hardcoded_commands_and_configured_search_params(monkeypatch):
    commands: list[list[str]] = []

    def fake_which(command: str) -> str:
        return f"/usr/bin/{command}"

    def fake_run(command, **kwargs):
        commands.append(command)

        class Completed:
            returncode = 0
            stdout = (
                "failed\tseed1\t99.0\t8\t0\t0\t1\t8\t1\t8\t1e-20\t100.0\n"
                if command[0] == "blastn"
                else ""
            )

        return Completed()

    monkeypatch.setattr("barcode_kit.blast.shutil.which", fake_which)
    monkeypatch.setattr("barcode_kit.blast.subprocess.run", fake_run)
    config = BlastRescueConfig(
        blastn_dust="yes",
        word_size=15,
        evalue=1e-5,
    )

    result = BlastRunner(config).rescue(
        [BlastRecord(accession_version="failed", sequence="AAAACCCC")],
        [BlastRecord(accession_version="seed1", sequence="AAAACCCC")],
        Marker.ITS,
    )

    assert commands[0][0] == "makeblastdb"
    assert commands[1][0] == "blastn"
    assert commands[1][commands[1].index("-dust") + 1] == "yes"
    assert commands[1][commands[1].index("-word_size") + 1] == "15"
    assert commands[1][commands[1].index("-evalue") + 1] == "1e-05"
    assert commands[1][commands[1].index("-outfmt") + 1].startswith("6 qseqid")
    assert result["failed"].sequence == Seq("AAAACCCC")


def _rescue_with_hits(
    hits: list[blast_module.BlastHit],
    marker: Marker,
    seed_lengths: dict[str, int],
    config: BlastRescueConfig | None = None,
):
    return blast_module._rescue_record(
        BlastRecord(accession_version="failed", sequence="N" * 100 + "A" * 2000),
        hits,
        marker,
        seed_lengths,
        config or BlastRescueConfig(),
    )

from __future__ import annotations

from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import barcode_kit.validation as validation_module
from barcode_kit.config import TreeShrinkConfig
from barcode_kit.models import BuildReportEntry, Marker, SequenceQuality
from barcode_kit import parser as parser_module
from barcode_kit.parser import (
    annotation_marker_evidence,
    detect_markers,
    extract_marker,
    parse_accession_version,
    parse_genbank_text,
)
from barcode_kit.phylogeny import TreeShrinkResult
from barcode_kit.validation import sequence_quality, tree_shrink_qc


def test_parser_exports_only_external_contract():
    assert parser_module.__all__ == [
        "AccessionVersion",
        "AnnotationMarkerEvidence",
        "ParsedGenBankRecord",
        "annotation_marker_evidence",
        "detect_markers",
        "extract_marker",
        "format_fasta_record",
        "parse_accession_version",
        "parse_genbank_file",
        "parse_genbank_text",
        "read_single_genbank",
    ]


def test_parse_accession_version():
    parsed = parse_accession_version("PP476489.4")
    assert parsed.root == "PP476489"
    assert parsed.version == 4
    assert parsed.value == "PP476489.4"


def test_parse_genbank_and_extract_rbcl(genbank_text):
    parsed = parse_genbank_text(genbank_text())
    assert parsed.accession.value == "PP476489.4"
    assert parsed.organism == "Iris japonica"
    assert parsed.taxon_id == 12345
    assert parsed.marker_flags[Marker.RBCL] is True
    sequence = extract_marker(parsed.record, Marker.RBCL)
    assert sequence is not None
    quality = sequence_quality(sequence)
    assert quality.length == 600
    assert quality.ambiguous_content == 0.0


def test_sequence_quality_counts_non_canonical_bases_as_ambiguous():
    quality = sequence_quality(Seq("ACGTNRYS"))

    assert quality.length == 8
    assert quality.gc_content == 0.5
    assert quality.ambiguous_content == 0.5
    assert not hasattr(quality, "has_stop_codon")
    assert not hasattr(quality, "has_frameshift")


def test_tree_shrink_qc_validation_filters_fasta_and_updates_report(
    tmp_path: Path,
    monkeypatch,
):
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "rbcl.fasta").write_text(
        ">KEEP.1|Iris_a\nACGT\n>DROP.1|Iris_b\nACGT\n",
        encoding="utf-8",
    )
    report = [
        _build_report_entry("KEEP.1", "Iris a", True, "KEEP.1|Iris_a"),
        _build_report_entry("DROP.1", "Iris b", True, "DROP.1|Iris_b"),
    ]
    removed_taxa = {"DROP.1|Iris_b"}

    class FakeAlignmentRunner:
        def align(self, input_path, output_path, *, threads):
            self.call = (input_path, output_path, threads)
            output_path.write_text(Path(input_path).read_text(encoding="utf-8"), encoding="utf-8")
            return output_path

    class FakeTreeRunner:
        def build_tree(self, input_path, output_path, *, bootstrap):
            self.call = (input_path, output_path, bootstrap)
            output_path.write_text("(KEEP:0.1,DROP:1.5);\n", encoding="utf-8")
            return output_path

    class FakeTreeShrinkRunner:
        def detect_outliers(
            self,
            tree_path,
            output_dir,
            *,
            output_prefix="output",
            quantile=0.05,
            max_removed=None,
        ):
            self.call = (tree_path, output_dir, output_prefix, quantile, max_removed)
            output_dir.mkdir(parents=True, exist_ok=True)
            removed_path = output_dir / f"{output_prefix}.txt"
            removed_path.write_text("\n".join(sorted(removed_taxa)), encoding="utf-8")
            return TreeShrinkResult(
                removed_taxa=removed_taxa,
                output_dir=output_dir,
                removed_taxa_path=removed_path,
            )

    alignment_runner = FakeAlignmentRunner()
    tree_runner = FakeTreeRunner()
    tree_shrink_runner = FakeTreeShrinkRunner()
    monkeypatch.setattr(validation_module, "AlignmentRunner", lambda: alignment_runner)
    monkeypatch.setattr(validation_module, "TreeRunner", lambda: tree_runner)
    monkeypatch.setattr(validation_module, "TreeShrinkRunner", lambda: tree_shrink_runner)

    updated_report = tree_shrink_qc(
        outdir,
        Marker.RBCL,
        report,
        TreeShrinkConfig(quantile=0.05, bootstrap=1000, max_removed=6),
    )

    assert updated_report[0].included is True
    assert updated_report[1].included is False
    assert updated_report[1].reason == "TreeShrink long-branch outlier"
    assert updated_report[1].metadata["tree_shrink_removed_taxa"] == str(
        outdir / "treeshrink_qc" / "treeshrink" / "output.txt"
    )
    fasta_text = (outdir / "rbcl.fasta").read_text(encoding="utf-8")
    assert "KEEP.1|Iris_a" in fasta_text
    assert "DROP.1|Iris_b" not in fasta_text
    assert alignment_runner.call == (
        outdir / "treeshrink_qc" / "input.fasta",
        outdir / "treeshrink_qc" / "mafft.fasta",
        1,
    )
    assert tree_runner.call == (
        outdir / "treeshrink_qc" / "mafft.fasta",
        outdir / "treeshrink_qc" / "iqtree.tree",
        1000,
    )
    assert tree_shrink_runner.call == (
        outdir / "treeshrink_qc" / "iqtree.tree",
        outdir / "treeshrink_qc" / "treeshrink",
        "output",
        0.05,
        6,
    )


def test_complete_its_single_feature_sets_flags_and_extracts_full_its():
    record = _record_with_features(
        "AAAACCCGGGG",
        [
            SeqFeature(
                FeatureLocation(0, 11, strand=1),
                type="misc_RNA",
                qualifiers={
                    "product": [
                        "internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2"
                    ]
                },
            )
        ],
    )

    flags = detect_markers(record)
    evidence = annotation_marker_evidence(record, Marker.ITS)

    assert flags[Marker.ITS] is True
    assert flags[Marker.ITS2] is True
    assert extract_marker(record, Marker.ITS) == Seq("AAAACCCGGGG")
    assert evidence.contains_marker is True
    assert evidence.extractable_marker is True
    assert evidence.annotation_pattern == "complete_single_feature"


def test_complete_its_separate_features_extract_from_its1_start_through_its2_end():
    record = _record_with_features(
        "TTAAAACCCGGGGTT",
        [
            SeqFeature(
                FeatureLocation(2, 6, strand=1),
                type="misc_feature",
                qualifiers={"note": ["ITS1"]},
            ),
            SeqFeature(
                FeatureLocation(6, 9, strand=1),
                type="rRNA",
                qualifiers={"product": ["5.8S ribosomal RNA"]},
            ),
            SeqFeature(
                FeatureLocation(9, 13, strand=1),
                type="misc_feature",
                qualifiers={"note": ["internal transcribed spacer 2"]},
            ),
        ],
    )

    flags = detect_markers(record)
    evidence = annotation_marker_evidence(record, Marker.ITS)

    assert flags[Marker.ITS] is True
    assert flags[Marker.ITS2] is True
    assert extract_marker(record, Marker.ITS) == Seq("AAAACCCGGGG")
    assert extract_marker(record, Marker.ITS2) == Seq("GGGG")
    assert evidence.contains_marker is True
    assert evidence.extractable_marker is True
    assert evidence.annotation_pattern == "complete_separate_features"


def test_broad_ssu_its_lsu_feature_extracts_annotation_span_for_its():
    record = _record_with_features(
        "TTTTAAAACCCGGGGTTTT",
        [
            SeqFeature(
                FeatureLocation(0, 19, strand=1),
                type="misc_RNA",
                qualifiers={
                    "product": [
                        "18S ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, "
                        "internal transcribed spacer 2, 28S ribosomal RNA"
                    ]
                },
            )
        ],
    )

    flags = detect_markers(record)
    evidence = annotation_marker_evidence(record, Marker.ITS)

    assert flags[Marker.ITS] is True
    assert flags[Marker.ITS2] is True
    assert extract_marker(record, Marker.ITS) == Seq("TTTTAAAACCCGGGGTTTT")
    assert evidence.contains_marker is True
    assert evidence.extractable_marker is True
    assert evidence.annotation_pattern == "broad_ssu_its_lsu"


def test_complete_its_prefers_clean_feature_when_broad_feature_is_also_present():
    record = _record_with_features(
        "TTTTAAAACCCGGGGTTTT",
        [
            SeqFeature(
                FeatureLocation(0, 19, strand=1),
                type="misc_RNA",
                qualifiers={
                    "product": [
                        "18S ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, "
                        "internal transcribed spacer 2, 28S ribosomal RNA"
                    ]
                },
            ),
            SeqFeature(
                FeatureLocation(4, 15, strand=1),
                type="misc_RNA",
                qualifiers={
                    "product": [
                        "internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2"
                    ]
                },
            ),
        ],
    )

    evidence = annotation_marker_evidence(record, Marker.ITS)

    assert extract_marker(record, Marker.ITS) == Seq("AAAACCCGGGG")
    assert evidence.annotation_pattern == "complete_single_feature"


def test_complete_its_can_use_its1_feature_and_combined_5_8s_its2_feature():
    record = _record_with_features(
        "TTAAAACCCGGGGTT",
        [
            SeqFeature(
                FeatureLocation(2, 6, strand=1),
                type="misc_feature",
                qualifiers={"note": ["ITS1"]},
            ),
            SeqFeature(
                FeatureLocation(6, 13, strand=1),
                type="misc_feature",
                qualifiers={"note": ["5.8S ribosomal RNA and internal transcribed spacer 2"]},
            ),
        ],
    )

    evidence = annotation_marker_evidence(record, Marker.ITS)

    assert extract_marker(record, Marker.ITS) == Seq("AAAACCCGGGG")
    assert evidence.contains_marker is True
    assert evidence.extractable_marker is True
    assert evidence.annotation_pattern == "complete_separate_features"


def test_its2_only_annotation_detects_its2_without_complete_its():
    record = _record_with_features(
        "GGGG",
        [
            SeqFeature(
                FeatureLocation(0, 4, strand=1),
                type="misc_feature",
                qualifiers={"product": ["internal transcribed spacer 2"]},
            )
        ],
    )

    flags = detect_markers(record)
    evidence = annotation_marker_evidence(record, Marker.ITS2)

    assert flags[Marker.ITS] is False
    assert flags[Marker.ITS2] is True
    assert extract_marker(record, Marker.ITS) is None
    assert extract_marker(record, Marker.ITS2) == Seq("GGGG")
    assert evidence.contains_marker is True
    assert evidence.extractable_marker is True
    assert evidence.annotation_pattern == "its2_only"


def test_partial_5_8s_its2_lsu_annotation_extracts_annotation_span_for_its2():
    record = _record_with_features(
        "CCCGGGGTTTT",
        [
            SeqFeature(
                FeatureLocation(0, 11, strand=1),
                type="misc_feature",
                qualifiers={
                    "note": [
                        "5.8S ribosomal RNA, internal transcribed spacer 2, large subunit ribosomal RNA"
                    ]
                },
            )
        ],
    )

    flags = detect_markers(record)
    evidence = annotation_marker_evidence(record, Marker.ITS2)

    assert flags[Marker.ITS] is False
    assert flags[Marker.ITS2] is True
    assert extract_marker(record, Marker.ITS2) == Seq("CCCGGGGTTTT")
    assert evidence.contains_marker is True
    assert evidence.extractable_marker is True
    assert evidence.annotation_pattern == "partial_5.8s_its2_lsu"


def test_its2_extracts_from_separate_subfeature_when_broad_partial_feature_is_present():
    record = _record_with_features(
        "CCCGGGGTTTT",
        [
            SeqFeature(
                FeatureLocation(0, 11, strand=1),
                type="misc_feature",
                qualifiers={
                    "note": [
                        "5.8S ribosomal RNA, internal transcribed spacer 2, large subunit ribosomal RNA"
                    ]
                },
            ),
            SeqFeature(
                FeatureLocation(3, 7, strand=1),
                type="misc_RNA",
                qualifiers={"product": ["ITS2"]},
            ),
        ],
    )

    evidence = annotation_marker_evidence(record, Marker.ITS2)

    assert extract_marker(record, Marker.ITS2) == Seq("GGGG")
    assert evidence.contains_marker is True
    assert evidence.extractable_marker is True


def test_source_qualifier_strings_do_not_create_its2_marker_evidence():
    record = _record_with_features("AAAA", [])
    record.features[0].qualifiers["note"] = ["BOLD identifier mentions ITS2"]

    flags = detect_markers(record)

    assert flags[Marker.ITS] is False
    assert flags[Marker.ITS2] is False
    assert extract_marker(record, Marker.ITS2) is None


def _record_with_features(sequence: str, features: list[SeqFeature]) -> SeqRecord:
    record = SeqRecord(Seq(sequence), id="TEST000001.1", name="TEST000001")
    record.annotations["molecule_type"] = "DNA"
    record.features = [
        SeqFeature(
            FeatureLocation(0, len(sequence), strand=1),
            type="source",
            qualifiers={"organism": ["Test species"], "db_xref": ["taxon:1"]},
        ),
        *features,
    ]
    return record


def _build_report_entry(
    accession_version: str,
    scientific_name: str,
    included: bool,
    output_id: str | None,
) -> BuildReportEntry:
    return BuildReportEntry(
        accession_version=accession_version,
        scientific_name=scientific_name,
        included=included,
        reason=None,
        quality=SequenceQuality(length=4, gc_content=0.5, ambiguous_content=0.0),
        output_id=output_id,
        metadata={"marker": Marker.RBCL.value},
    )

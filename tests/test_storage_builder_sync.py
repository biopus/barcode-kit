from __future__ import annotations

import io
import json
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import barcode_kit.builder as builder_module
from barcode_kit.builder import build_dataset
from barcode_kit.config import (
    AppConfig,
    BlastRescueConfig,
    CollectorConfig,
    ItsxrustConfig,
    TreeShrinkConfig,
)
from barcode_kit.genbank import DownloadItem, DownloadReport, SyncService
from barcode_kit.blast import BlastRecord, BlastRescueResult
from barcode_kit.itsxrust import ItsxrustExtractionResult
from barcode_kit.models import (
    GenBankCacheRecord,
    Marker,
    TaxonConstraint,
    TaxonQuery,
    TaxonomyRecord,
)
from barcode_kit.storage import Storage


class FakeClient:
    def __init__(self, records: dict[str, str]):
        self.records = records
        self.fetch_calls = 0
        self.download_calls = 0
        self.search_terms: list[str] = []

    def search_accessions(self, term: str) -> list[str]:
        self.search_terms.append(term)
        return list(self.records)

    def fetch_records(self, accessions: list[str], output_dir: Path) -> DownloadReport:
        self.fetch_calls += 1
        items = {}
        for accession in accessions:
            path = output_dir / f"{accession}.gb"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(self.records[accession], encoding="utf-8")
            items[accession] = DownloadItem(
                accession=accession,
                status="succeeded",
            )
        return DownloadReport(items=items)


class FakeResolver:
    def standardize(self, scientific_name: str, taxon_id_hint: int | None = None) -> TaxonomyRecord:
        return TaxonomyRecord(
            taxon_id=taxon_id_hint or 12345,
            scientific_name=scientific_name,
            kingdom="Viridiplantae",
            phylum="Streptophyta",
            class_name="Magnoliopsida",
            order="Asparagales",
            family="Iridaceae",
            genus=scientific_name.split()[0],
            species=scientific_name.split()[1],
        )


def test_sync_ingests_new_record_and_skips_current(tmp_path: Path, genbank_text):
    config = _config(tmp_path)
    storage = Storage(config.database_path)
    client = FakeClient({"PP476489.4": genbank_text()})
    service = SyncService(config, storage, FakeResolver(), client)

    first = service.sync(TaxonQuery("genus", "Iris"), Marker.RBCL)
    assert first.downloaded == 1
    assert first.skipped == 0
    assert not first.failed
    assert (config.genbank_cache_dir / "PP476489.4.gb").exists()

    second = service.sync(TaxonQuery("genus", "Iris"), Marker.RBCL)
    assert second.downloaded == 0
    assert second.skipped == 1
    assert client.fetch_calls == 1


def test_sync_reuses_valid_cache_file_when_database_is_missing(tmp_path: Path, genbank_text):
    config = _config(tmp_path)
    storage = Storage(config.database_path)
    cache_path = config.genbank_cache_dir / "PP476489.4.gb"
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(genbank_text(), encoding="utf-8")
    client = FakeClient({"PP476489.4": genbank_text()})
    service = SyncService(config, storage, FakeResolver(), client)

    result = service.sync(TaxonQuery("genus", "Iris"), Marker.RBCL)

    assert result.downloaded == 0
    assert result.reused_local == 1
    assert result.ingested == 1
    assert result.failed == []
    assert client.fetch_calls == 0


def test_sync_uses_ncbi_search_term_for_its_marker(tmp_path: Path):
    config = _config(tmp_path)
    storage = Storage(config.database_path)
    client = FakeClient({})
    service = SyncService(config, storage, FakeResolver(), client)

    result = service.sync(TaxonQuery("genus", "Iris"), Marker.ITS)

    assert result.remote_count == 0
    assert client.search_terms == ['Iris[Organism] AND "internal transcribed spacer"']


def test_sync_uses_expanded_taxid_organism_term(tmp_path: Path):
    config = _config(tmp_path)
    storage = Storage(config.database_path)
    client = FakeClient({})
    service = SyncService(config, storage, FakeResolver(), client)

    result = service.sync(TaxonQuery("taxid", "58920"), Marker.RBCL)

    assert result.remote_count == 0
    assert client.search_terms == ["txid58920[Organism:exp] AND rbcl"]


def test_sync_adds_lineage_constraints_to_ncbi_search_term(tmp_path: Path):
    config = _config(tmp_path)
    storage = Storage(config.database_path)
    client = FakeClient({})
    service = SyncService(config, storage, FakeResolver(), client)
    query = TaxonQuery(
        "genus",
        "Iris",
        constraints=(TaxonConstraint(rank="kingdom", name="Viridiplantae"),),
    )

    result = service.sync(query, Marker.RBCL)

    assert result.remote_count == 0
    assert client.search_terms == ["Iris[Organism] AND Viridiplantae[Organism] AND rbcl"]


def test_build_dataset_exports_raw_dataset_manifest_and_compact_report(
    tmp_path: Path,
    genbank_text,
):
    config = _config(tmp_path)
    storage = Storage(config.database_path)
    service = SyncService(
        config,
        storage,
        FakeResolver(),
        FakeClient({"PP476489.4": genbank_text()}),
    )
    service.sync(TaxonQuery("genus", "Iris"), Marker.RBCL)

    outdir = tmp_path / "out"
    report = build_dataset(
        config,
        storage,
        TaxonQuery("genus", "Iris"),
        Marker.RBCL,
        outdir,
    )

    assert len(report) == 1
    assert report[0].included is True
    assert "PP476489.4|Iris_japonica" in (outdir / "raw" / "rbcl.fasta").read_text()
    assert json.loads((outdir / "dataset.json").read_text(encoding="utf-8")) == {
        "format_version": 1,
        "marker": "rbcl",
        "raw_fasta": "raw/rbcl.fasta",
        "build_report": "build_report.json",
    }
    assert json.loads((outdir / "build_report.json").read_text(encoding="utf-8")) == [
        {
            "accession": "PP476489.4",
            "sequence_id": "PP476489.4|Iris_japonica",
            "scientific_name": "Iris japonica",
            "infraspecific_rank": None,
            "is_hybrid": False,
            "is_uncertain": False,
            "extraction_backend": "annotation",
            "error": None,
        }
    ]


def test_build_dataset_keeps_ambiguous_sequences_for_later_qc(
    tmp_path: Path,
    genbank_text,
):
    config = _config(tmp_path)
    storage = Storage(config.database_path)
    service = SyncService(
        config,
        storage,
        FakeResolver(),
        FakeClient({"AMB000001.1": genbank_text(accession="AMB000001", version=1, sequence="ACGTNNRY")}),
    )
    service.sync(TaxonQuery("genus", "Iris"), Marker.RBCL)

    report = build_dataset(
        config,
        storage,
        TaxonQuery("genus", "Iris"),
        Marker.RBCL,
        tmp_path / "out",
    )

    assert len(report) == 1
    assert report[0].included is True
    assert report[0].reason is None
    assert report[0].quality is None
    assert "ACGTNNRY" in (tmp_path / "out" / "raw" / "rbcl.fasta").read_text()


def test_build_dataset_keeps_infraspecific_records_for_later_qc(tmp_path: Path, genbank_text):
    config = _config(tmp_path)
    storage = Storage(config.database_path)
    storage.initialize()
    config.genbank_cache_dir.mkdir(parents=True, exist_ok=True)
    (config.genbank_cache_dir / "VAR000001.1.gb").write_text(
        genbank_text(accession="VAR000001", version=1, organism="Iris japonica var. alba", taxon_id=111),
        encoding="utf-8",
    )
    (config.genbank_cache_dir / "SPC000001.1.gb").write_text(
        genbank_text(accession="SPC000001", version=1, organism="Iris japonica", taxon_id=222),
        encoding="utf-8",
    )
    with storage.connect() as connection:
        storage.upsert_taxonomy(
            TaxonomyRecord(
                taxon_id=111,
                scientific_name="Iris japonica var. alba",
                family="Iridaceae",
                genus="Iris",
                species="japonica",
                infraspecific_rank="variety",
            ),
            connection,
        )
        storage.upsert_taxonomy(
            TaxonomyRecord(
                taxon_id=222,
                scientific_name="Iris japonica",
                family="Iridaceae",
                genus="Iris",
                species="japonica",
            ),
            connection,
        )
        storage.upsert_genbank_cache(
            GenBankCacheRecord("VAR000001", 1, "VAR000001.1", 111, has_rbcl=True),
            connection,
        )
        storage.upsert_genbank_cache(
            GenBankCacheRecord("SPC000001", 1, "SPC000001.1", 222, has_rbcl=True),
            connection,
        )

    report = build_dataset(
        config,
        storage,
        TaxonQuery("genus", "Iris"),
        Marker.RBCL,
        tmp_path / "out",
    )

    by_accession = {entry.accession_version: entry for entry in report}
    assert by_accession["VAR000001.1"].included is True
    assert by_accession["SPC000001.1"].included is True
    build_report = json.loads(
        (tmp_path / "out" / "build_report.json").read_text(encoding="utf-8")
    )
    report_by_accession = {entry["accession"]: entry for entry in build_report}
    assert report_by_accession["VAR000001.1"]["infraspecific_rank"] == "variety"


def test_build_dataset_exports_raw_records_without_treeshrink_qc(
    tmp_path: Path,
    genbank_text,
):
    config = _config(tmp_path)
    storage = Storage(config.database_path)
    service = SyncService(
        config,
        storage,
        FakeResolver(),
        FakeClient(
            {
                "PP476489.4": genbank_text(),
                "PP476490.1": genbank_text(accession="PP476490", version=1),
            }
        ),
    )
    service.sync(TaxonQuery("genus", "Iris"), Marker.RBCL)

    report = build_dataset(
        config,
        storage,
        TaxonQuery("genus", "Iris"),
        Marker.RBCL,
        tmp_path / "out",
    )

    assert len(report) == 2
    assert report[0].included is True
    assert report[1].included is True
    fasta_text = (tmp_path / "out" / "raw" / "rbcl.fasta").read_text(encoding="utf-8")
    assert "PP476489.4|Iris_japonica" in fasta_text
    assert "PP476490.1|Iris_japonica" in fasta_text


def test_build_its_uses_annotation_when_all_itsxrust_extractions_fail(
    tmp_path: Path,
    monkeypatch,
):
    config = _config(tmp_path)
    storage = Storage(config.database_path)
    storage.initialize()
    _insert_cached_record(storage, "ITS200001.1", has_its=True, has_its2=True)
    _write_genbank(
        config.genbank_cache_dir / "ITS200001.1.gb",
        "ITS200001",
        1,
        "GGGG",
        [
            SeqFeature(
                FeatureLocation(0, 4, strand=1),
                type="misc_feature",
                qualifiers={"product": ["internal transcribed spacer 2"]},
            )
        ],
    )
    itsxrust_runner = FailingItsxrustRunner()
    monkeypatch.setattr(builder_module, "ItsxrustRunner", lambda config: itsxrust_runner)
    monkeypatch.setattr(
        builder_module,
        "BlastRunner",
        lambda config: (_ for _ in ()).throw(AssertionError("BLAST should not run")),
    )

    report = build_dataset(
        config,
        storage,
        TaxonQuery("genus", "Iris"),
        Marker.ITS,
        tmp_path / "out",
    )

    assert itsxrust_runner.extract_many_calls == 1
    assert len(report) == 1
    assert report[0].included is False
    assert report[0].reason == "marker_not_extracted"
    assert report[0].metadata["extraction_backend"] == "annotation"
    assert report[0].metadata["hmm_fallback_reason"] == "hmm_failed"
    assert report[0].metadata["fallback_reason"] is None
    assert report[0].metadata["annotation_contains_marker"] is False
    assert report[0].metadata["annotation_extractable_marker"] is False


def test_build_its_rescues_failed_itsxrust_record_with_seed(
    tmp_path: Path,
    monkeypatch,
):
    config = _config(tmp_path)
    storage = Storage(config.database_path)
    storage.initialize()
    _insert_cached_record(storage, "ITS000001.1", has_its=True, has_its2=True)
    _insert_cached_record(storage, "ITS000002.1", has_its=True, has_its2=True)
    _write_genbank(
        config.genbank_cache_dir / "ITS000001.1.gb",
        "ITS000001",
        1,
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
    _write_genbank(
        config.genbank_cache_dir / "ITS000002.1.gb",
        "ITS000002",
        1,
        "TTTTAAAACCCGGGGAAAA",
        [],
    )
    blast_runner = FakeBlastRunner(
        {
            "ITS000002.1": BlastRescueResult(
                sequence=Seq("AAAACCCGGGG"),
                fallback_reason=None,
                metadata={
                    "blast_seed_accession": "ITS000001.1",
                    "blast_identity": 95.0,
                    "blast_subject_coverage": 1.0,
                    "blast_query_start": 5,
                    "blast_query_end": 15,
                    "blast_bitscore": 200.0,
                    "blast_strand": "+",
                },
            )
        }
    )
    itsxrust_runner = BatchSeedThenFailItsxrustRunner()
    monkeypatch.setattr(builder_module, "ItsxrustRunner", lambda config: itsxrust_runner)
    monkeypatch.setattr(builder_module, "BlastRunner", lambda config: blast_runner)

    report = build_dataset(
        config,
        storage,
        TaxonQuery("genus", "Iris"),
        Marker.ITS,
        tmp_path / "out",
    )

    assert itsxrust_runner.extract_many_calls == 1
    assert len(report) == 2
    assert report[0].included is True
    assert report[0].metadata["extraction_backend"] == "itsxrust"
    assert report[1].included is True
    assert report[1].metadata["extraction_backend"] == "blastn"
    assert report[1].metadata["fallback_reason"] is None
    assert report[1].metadata["blast_seed_accession"] == "ITS000001.1"
    assert "ITS000001.1|Iris_japonica" in (
        tmp_path / "out" / "raw" / "its.fasta"
    ).read_text()
    assert "ITS000002.1|Iris_japonica" in (
        tmp_path / "out" / "raw" / "its.fasta"
    ).read_text()


def test_build_its_marks_unrescued_partial_itsxrust_failures_as_extraction_errors(
    tmp_path: Path,
    monkeypatch,
):
    config = _config(tmp_path)
    storage = Storage(config.database_path)
    storage.initialize()
    _insert_cached_record(storage, "ITS000001.1", has_its=True, has_its2=True)
    _insert_cached_record(storage, "ITS000003.1", has_its=True, has_its2=True)
    _write_genbank(config.genbank_cache_dir / "ITS000001.1.gb", "ITS000001", 1, "AAAACCCCGGGG", [])
    _write_genbank(config.genbank_cache_dir / "ITS000003.1.gb", "ITS000003", 1, "TTTTAAAACCCCGGGG", [])

    itsxrust_runner = BatchSeedThenFailItsxrustRunner()
    blast_runner = FakeBlastRunner({})
    monkeypatch.setattr(builder_module, "ItsxrustRunner", lambda config: itsxrust_runner)
    monkeypatch.setattr(builder_module, "BlastRunner", lambda config: blast_runner)

    report = build_dataset(
        config,
        storage,
        TaxonQuery("genus", "Iris"),
        Marker.ITS,
        tmp_path / "out",
    )

    assert itsxrust_runner.extract_many_calls == 1
    assert len(report) == 2
    assert report[0].included is True
    assert report[1].included is False
    assert report[1].reason == "marker_not_extracted"
    assert report[1].metadata["extraction_backend"] == "blastn"
    assert report[1].metadata["fallback_reason"] == "no_blast_hit"


def test_build_its_all_itsxrust_failures_use_annotation_without_blast(
    tmp_path: Path,
    monkeypatch,
):
    config = _config(tmp_path)
    storage = Storage(config.database_path)
    storage.initialize()
    _insert_cached_record(storage, "ITS000002.1", has_its=True, has_its2=True)
    _write_genbank(
        config.genbank_cache_dir / "ITS000002.1.gb",
        "ITS000002",
        1,
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

    _insert_cached_record(storage, "ITS000004.1", has_its=True, has_its2=True)
    _write_genbank(
        config.genbank_cache_dir / "ITS000004.1.gb",
        "ITS000004",
        1,
        "CCCCAAAAGGGG",
        [],
    )
    itsxrust_runner = FailingItsxrustRunner()
    monkeypatch.setattr(builder_module, "ItsxrustRunner", lambda config: itsxrust_runner)
    monkeypatch.setattr(
        builder_module,
        "BlastRunner",
        lambda config: (_ for _ in ()).throw(AssertionError("BLAST should not run")),
    )

    report = build_dataset(
        config,
        storage,
        TaxonQuery("genus", "Iris"),
        Marker.ITS,
        tmp_path / "out",
    )

    assert itsxrust_runner.extract_many_calls == 1
    assert len(report) == 2
    assert report[0].included is True
    assert report[0].metadata["extraction_backend"] == "annotation"
    assert report[0].metadata["hmm_fallback_reason"] == "hmm_failed"
    assert report[1].included is False
    assert report[1].reason == "marker_not_extracted"
    assert report[1].metadata["extraction_backend"] == "annotation"


def test_build_uses_configured_blast_rescue_defaults_for_default_runner(
    tmp_path: Path,
    monkeypatch,
):
    captured_configs: list[BlastRescueConfig] = []

    class CapturingBlastRunner:
        def __init__(self, config: BlastRescueConfig):
            captured_configs.append(config)

        def rescue(self, failed_records, seeds, marker):
            return {
                record.accession_version: BlastRescueResult(
                    sequence=None,
                    fallback_reason="no_blast_hit",
                )
                for record in failed_records
            }

    config = _config(
        tmp_path,
        blast_rescue=BlastRescueConfig(blastn_dust="yes"),
    )
    storage = Storage(config.database_path)
    storage.initialize()
    _insert_cached_record(storage, "ITS000001.1", has_its=True, has_its2=True)
    _insert_cached_record(storage, "ITS000005.1", has_its=True, has_its2=True)
    _write_genbank(config.genbank_cache_dir / "ITS000001.1.gb", "ITS000001", 1, "AAAACCCCGGGG", [])
    _write_genbank(config.genbank_cache_dir / "ITS000005.1.gb", "ITS000005", 1, "AAAACCCCGGGG", [])
    monkeypatch.setattr(builder_module, "BlastRunner", CapturingBlastRunner)
    monkeypatch.setattr(
        builder_module,
        "ItsxrustRunner",
        lambda config: BatchSeedThenFailItsxrustRunner(),
    )

    build_dataset(
        config,
        storage,
        TaxonQuery("genus", "Iris"),
        Marker.ITS,
        tmp_path / "out",
    )

    assert captured_configs == [config.blast_rescue]


def test_build_uses_configured_itsxrust_defaults_for_default_runner(
    tmp_path: Path,
    monkeypatch,
):
    captured_configs: list[ItsxrustConfig] = []

    class CapturingItsxrustRunner:
        def __init__(self, config: ItsxrustConfig):
            captured_configs.append(config)

        def extract_many(self, records, marker, hmm_path):
            return {
                item.accession_version: ItsxrustExtractionResult(
                    sequence=None,
                    fallback_reason="hmm_failed",
                )
                for item in records
            }

    config = _config(
        tmp_path,
        itsxrust=ItsxrustConfig(max_per_anchor=25),
    )
    storage = Storage(config.database_path)
    storage.initialize()
    _insert_cached_record(storage, "ITS000006.1", has_its=True, has_its2=True)
    _write_genbank(config.genbank_cache_dir / "ITS000006.1.gb", "ITS000006", 1, "AAAACCCCGGGG", [])
    monkeypatch.setattr(builder_module, "ItsxrustRunner", CapturingItsxrustRunner)

    build_dataset(
        config,
        storage,
        TaxonQuery("genus", "Iris"),
        Marker.ITS,
        tmp_path / "out",
    )

    assert captured_configs == [config.itsxrust]


def _config(
    tmp_path: Path,
    *,
    blast_rescue: BlastRescueConfig | None = None,
    itsxrust: ItsxrustConfig | None = None,
    tree_shrink_qc: TreeShrinkConfig | None = None,
) -> AppConfig:
    return AppConfig(
        data_dir=tmp_path / "barcode-kit",
        collectors=CollectorConfig(
            download_workers=1,
            genbank_email="test@example.com",
        ),
        blast_rescue=blast_rescue or BlastRescueConfig(),
        itsxrust=itsxrust or ItsxrustConfig(),
        tree_shrink_qc=tree_shrink_qc or TreeShrinkConfig(),
    )


class FailingItsxrustRunner:
    def __init__(self):
        self.extract_many_calls = 0

    def extract_many(self, records, marker, hmm_path) -> dict[str, ItsxrustExtractionResult]:
        self.extract_many_calls += 1
        return {
            item.accession_version: ItsxrustExtractionResult(sequence=None, fallback_reason="hmm_failed")
            for item in records
        }


class BatchSeedThenFailItsxrustRunner:
    def __init__(self):
        self.extract_many_calls = 0

    def extract_many(self, records, marker, hmm_path) -> dict[str, ItsxrustExtractionResult]:
        self.extract_many_calls += 1
        return {
            item.accession_version: (
                ItsxrustExtractionResult(sequence=Seq("AAAACCCGGGG"))
                if item.accession_version == "ITS000001.1"
                else ItsxrustExtractionResult(sequence=None, fallback_reason="no_anchor_full")
            )
            for item in records
        }


class FakeBlastRunner:
    def __init__(self, results: dict[str, BlastRescueResult]):
        self.results = results
        self.seed_accessions: list[str] = []

    def rescue(
        self,
        failed_records,
        seeds: list[BlastRecord],
        marker: Marker,
    ) -> dict[str, BlastRescueResult]:
        self.seed_accessions = [seed.accession_version for seed in seeds]
        return {
            record.accession_version: self.results.get(
                record.accession_version,
                BlastRescueResult(sequence=None, fallback_reason="no_blast_hit"),
            )
            for record in failed_records
        }


class FakeAlignmentRunner:
    def __init__(self):
        self.calls = []

    def align(self, input_path, output_path, *, program, threads):
        input_path = Path(input_path)
        output_path = Path(output_path)
        self.calls.append((input_path, output_path, program, threads))
        output_path.write_text(input_path.read_text(encoding="utf-8"), encoding="utf-8")
        return output_path


class FakeTreeRunner:
    def __init__(self):
        self.calls = []

    def build_tree(self, input_path, output_path, *, bootstrap):
        input_path = Path(input_path)
        output_path = Path(output_path)
        self.calls.append((input_path, output_path, bootstrap))
        output_path.write_text("(PP476489.4|Iris_japonica:0.1,PP476490.1|Iris_japonica:1.5);\n", encoding="utf-8")
        return output_path


class FakeTreeShrinkRunner:
    def __init__(self, removed_taxa: set[str]):
        self.removed_taxa = removed_taxa
        self.calls = []

    def detect_outliers(
        self,
        tree_path,
        output_dir,
        *,
        output_prefix="output",
        quantile=0.05,
        max_removed=None,
    ):
        tree_path = Path(tree_path)
        output_dir = Path(output_dir)
        self.calls.append((tree_path, output_dir, output_prefix, quantile, max_removed))
        output_dir.mkdir(parents=True, exist_ok=True)
        removed_path = output_dir / f"{output_prefix}.txt"
        removed_path.write_text("\t".join(sorted(self.removed_taxa)) + "\n", encoding="utf-8")
        return TreeShrinkResult(
            removed_taxa=self.removed_taxa,
            output_dir=output_dir,
            removed_taxa_path=removed_path,
        )


def _insert_cached_record(
    storage: Storage,
    accession_version: str,
    *,
    has_its: bool,
    has_its2: bool,
) -> None:
    accession_root, version_text = accession_version.split(".", 1)
    taxonomy = TaxonomyRecord(
        taxon_id=12345,
        scientific_name="Iris japonica",
        kingdom="Viridiplantae",
        family="Iridaceae",
        genus="Iris",
        species="japonica",
    )
    with storage.connect() as connection:
        storage.upsert_taxonomy(taxonomy, connection)
        storage.upsert_genbank_cache(
            GenBankCacheRecord(
                accession_root=accession_root,
                version=int(version_text),
                accession_version=accession_version,
                taxon_id=taxonomy.taxon_id,
                has_its=has_its,
                has_its2=has_its2,
            ),
            connection,
        )


def _write_genbank(
    path: Path,
    accession: str,
    version: int,
    sequence: str,
    features: list[SeqFeature],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    record = SeqRecord(
        Seq(sequence),
        id=f"{accession}.{version}",
        name=accession,
        description="Iris japonica ITS test record",
    )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["accessions"] = [accession]
    record.annotations["sequence_version"] = version
    record.annotations["organism"] = "Iris japonica"
    record.features = [
        SeqFeature(
            FeatureLocation(0, len(sequence), strand=1),
            type="source",
            qualifiers={"organism": ["Iris japonica"], "db_xref": ["taxon:12345"]},
        ),
        *features,
    ]
    handle = io.StringIO()
    SeqIO.write(record, handle, "genbank")
    path.write_text(handle.getvalue(), encoding="utf-8")

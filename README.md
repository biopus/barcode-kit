# barcode-kit

`barcode-kit` 是一个本地 GenBank-backed DNA barcode 数据收集与构建工具。它把公共数据库检索、GenBank 记录缓存、taxonomy 标准化、marker 检测和 FASTA 数据集导出串成一个可重复的本地流程。

当前版本聚焦 GenBank 单数据源闭环：

- 从 NCBI GenBank 检索指定 taxon 和 marker 的 accession。
- 下载并缓存 GenBank 原始记录。
- 使用 ETE `NCBITaxa` 标准化 NCBI taxonomy。
- 在本地 SQLite 中维护 taxonomy 与 accession 级缓存索引。
- 从缓存构建指定 marker 的 FASTA 数据集和 JSON 报告。

## 安装

本项目使用 `uv` 管理环境和依赖。

```bash
uv sync --dev
```

安装后可通过项目脚本访问 CLI：

```bash
uv run barcode-kit --help
```

## 配置

首次使用前需要设置 NCBI E-utilities email。建议在开发和测试时用 `BARCODE_KIT_CONFIG` 指向临时配置，避免写入真实 home 目录。

```bash
BARCODE_KIT_CONFIG=/tmp/barcode-kit.toml \
  uv run barcode-kit config set genbank.email user@example.com
```

可选配置项：

```bash
uv run barcode-kit config set genbank.api_key YOUR_API_KEY
uv run barcode-kit config set paths.data_dir /path/to/barcode-kit-data
uv run barcode-kit config set collectors.batch_size 500
uv run barcode-kit config set collectors.timeout 30
uv run barcode-kit config set collectors.retry_attempts 3
uv run barcode-kit config set build.tree_shrink_qc.quantile 0.05
uv run barcode-kit config set build.tree_shrink_qc.bootstrap 1000
uv run barcode-kit config set build.tree_shrink_qc.max_removed auto-select
```

查看当前有效配置：

```bash
uv run barcode-kit config list
```

默认数据目录是 `~/.barcode-kit`，其中包含：

- `database.db`：SQLite 元数据数据库。
- `cache/genbank/`：下载的 GenBank 原始记录。
- `logs/`：同步日志。

## Taxonomy 数据

`barcode-kit` 使用 `ete3.NCBITaxa` 读取本地 NCBI taxonomy 数据库。首次初始化 ETE taxonomy 数据库时，ETE 可能需要下载并构建本地 taxonomy 数据；如果在离线环境中运行，需要提前准备 ETE 可用的本地数据库。

同步阶段会优先使用 GenBank 记录中的 NCBI TaxId hint；没有 TaxId 时会通过科学名解析 TaxId。解析后的标准名和 lineage 会写入本地 `taxonomy` 表。

## 常用命令

### 同步 GenBank 记录

每次同步必须且只能指定一个 taxon 条件：`--family`、`--genus` 或 `--species`。

```bash
uv run barcode-kit sync --genus Iris --marker rbcl
uv run barcode-kit sync --species "Iris japonica" --marker its
uv run barcode-kit sync --family Iridaceae --marker matk
```

支持的 marker：

- `its`
- `its2`
- `matk`
- `rbcl`

同步命令会输出 JSON，包含远端命中数、下载数、跳过数、更新数和失败条目。

### 构建 FASTA 数据集

```bash
uv run barcode-kit build --genus Iris --marker rbcl --outdir ./out
```

可在构建时应用质量和 taxonomy 过滤：

```bash
uv run barcode-kit build --genus Iris --marker rbcl \
  --min-length 500 \
  --max-ambiguous-content 0.05 \
  --exclude-hybrid \
  --exclude-uncertain \
  --outdir ./out
```

如需在构建后用 TreeShrink 去除长枝异常序列，可显式开启系统发育质控。使用前需确保 `mafft`、`iqtree` 和 `run_treeshrink.py` 已安装并可在 `PATH` 中找到。该流程固定使用 IQ-TREE `-m MFP -T AUTO` 建树，并以 TreeShrink `per-gene` 模式检测异常序列：

```bash
uv run barcode-kit build --genus Iris --marker rbcl \
  --tree-shrink-qc \
  --outdir ./out
```

TreeShrink 的 `-k` 上限通过 `build.tree_shrink_qc.max_removed` 配置。默认值为 `auto-select`，表示沿用 TreeShrink 根据数据自动选择的上限；需要更激进时可设置为正整数，例如：

```bash
uv run barcode-kit config set build.tree_shrink_qc.max_removed 6
```

输出目录中会生成：

- `<marker>.fasta`：构建出的 FASTA 数据集。
- `build_report.json`：每条候选记录的纳入状态、排除原因和质量指标。

### 系统发育分析 Python API

`barcode_kit.phylogeny` 提供外部工具封装接口，可在 Python 流程中调用多序列比对、trimAl 裁剪、系统发育树构建和 TreeShrink 长枝异常序列识别。使用前需确保对应命令已安装并可在 `PATH` 中找到，例如 `mafft`、`iqtree`、`trimal` 和 `run_treeshrink.py`。

```python
from pathlib import Path

from barcode_kit.phylogeny import (
    AlignmentProgram,
    SubprocessAlignmentRunner,
    SubprocessTreeRunner,
    SubprocessTreeShrinkRunner,
    SubprocessTrimalRunner,
    run_tree_shrink_qc,
)

aligned = SubprocessAlignmentRunner().align(
    Path("rbcl.fasta"),
    Path("rbcl.aligned.fasta"),
    program=AlignmentProgram.MAFFT,
    threads=4,
)
trimmed = SubprocessTrimalRunner().trim(aligned, Path("rbcl.trimmed.fasta"))
SubprocessTreeRunner().build_tree(
    trimmed,
    Path("rbcl.tree"),
    bootstrap=1000,
)

run_tree_shrink_qc(
    Path("rbcl.fasta"),
    Path("rbcl.filtered.fasta"),
    Path("treeshrink-qc"),
    bootstrap=1000,
    max_removed=6,
    tree_shrink_runner=SubprocessTreeShrinkRunner(),
)
```

`build_dataset()` 也可以显式开启 TreeShrink 质控。开启后会先输出候选 FASTA 到临时工作目录，依次运行 MAFFT、IQ-TREE 和 TreeShrink，再把 TreeShrink 标记的异常序列从最终 `<marker>.fasta` 中删除，并在 `build_report.json` 中记录排除原因。

### 查看本地缓存

```bash
uv run barcode-kit db status
uv run barcode-kit db info --genus Iris
uv run barcode-kit db info --family Iridaceae
uv run barcode-kit db info --rank family
uv run barcode-kit db info --rank genus --family Asparagaceae
uv run barcode-kit db info --rank species --genus Aspidistra
```

不带 `--rank` 时，`db info` 输出整体或指定 taxon 的 marker 覆盖统计。带 `--rank family|genus|species` 时，输出当前缓存中可操作的分类阶元列表，每项包含名称、缓存记录数和 marker 覆盖数，便于确认后再执行 `db remove`。

### 管理本地缓存

删除缓存是显式本地操作，不会由 `sync` 隐式触发。`remove` 删除匹配的 accession 级数据库记录和对应 GenBank 文件；`clear` 清空本地缓存；`prune` 清理数据库与文件系统之间的不一致状态。

```bash
uv run barcode-kit db remove --accession PP476489.4 --yes
uv run barcode-kit db remove --genus Iris --yes
uv run barcode-kit db clear --yes
uv run barcode-kit db prune --yes
```

这些命令会输出 JSON 报告，包含删除的数据库记录数、GenBank 文件数和孤儿 taxonomy 记录数。不加 `--yes` 时会要求确认。

## 开发

常用开发命令：

```bash
uv sync --dev
uv run pytest
uv run barcode-kit --help
uv build
```

测试中涉及配置、缓存或数据库路径时，应使用 `tmp_path`、`monkeypatch` 或显式设置 `BARCODE_KIT_CONFIG`，不要写入真实用户目录。

## 项目结构

```text
src/barcode_kit/
  cli.py          Typer CLI
  config.py       配置加载与写入
  genbank.py      NCBI ESearch/EFetch 同步流程
  taxonomy.py     ETE NCBITaxa taxonomy 标准化
  parser.py       GenBank 记录解析与 marker 检测
  storage.py      SQLite schema 与查询
  builder.py      FASTA 构建与报告生成
  models.py       领域模型

tests/            pytest 测试
concepts/         架构与设计说明
ref/              外部 API 和工具参考
imp/              实现参考，不是生产包代码
```

## 注意事项

- 不要提交本地缓存数据库、下载的 GenBank 记录、NCBI API key 或个人 email 配置。
- `sync` 会访问 NCBI E-utilities，并依赖本地 ETE taxonomy 数据。
- `build` 默认只使用本地缓存，不会隐式触发远端同步。

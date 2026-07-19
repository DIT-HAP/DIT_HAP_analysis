# 硬编码常量归纳到 config 设计

**日期：** 2026-07-19
**范围：** 中等 — 把两组会影响分析结果的"魔法数字"从脚本常量提升为可配置参数，走既有的 config→params→argparse 通路。

## 1. 背景

`config/analysis.yaml` 已经承载了 clustering / enrichment / ml 三个 section 的实验参数（k 值范围、FDR 阈值、DR 分割阈值等）。但代码里还散落着几类硬编码常量，其中两组符合"应该进 config"的判断标准（会影响分析结果、值本身有讨论空间）：

1. `workflow/src/clustering/candidates.py` 的 `DR_CAP = 1.3` / `DL_DIVISOR = 10`（`scale_features` 用于裁剪/缩放 DR、DL）。
2. `workflow/scripts/enrichment/run_network_enrichment.py` 的 `REVIGO_CUTOFFS = [0.7, 0.5]`（REVIGO 两轮语义相似度截断阈值）。

不进 config 的常量（明确排除）：
- `_GAF_HEADER_DATE`（`ontology.py`）——固定可复现性时间戳，不是分析参数。
- `STRING_API_URL` / `REVIGO_URL` / `STRING_CALLER_IDENTITY` 等——API 端点/身份标识，不是分析参数。
- `STRING_MAX_RETRIES` / `STRING_RETRY_DELAY`——运维参数，本次不动（用户选择了"中等范围"，排除了"更大范围"选项里的这一项）。

## 2. 现有模式（不引入新模式，照抄）

以 `enrichment.wt_cluster` 为例，现有链路是：

```
模块常量（默认值）
  → dataclass 字段 default 引用常量
  → argparse --flag default 引用常量
  → Snakemake rule params: xxx=config.get("section", {}).get("key", default)
  → shell 命令 --flag {params.xxx}
```

`config/analysis.yaml` 的值始终是"最终默认值来源"，脚本里的模块常量退化为"CLI 直接调用时的默认值 + 单元测试用的已知值"，两边保持一致但 config 优先。

## 3. 改动清单

### 3.1 `config/analysis.yaml`

```yaml
clustering:
  n_clusters: 64
  random_state: 42
  k_min: 2
  k_max: 20
  dr_cap: 1.3        # DR 超过此值裁剪为此值（byte-faithful quirk）
  dl_divisor: 10     # DL 在聚类前除以该值

enrichment_network:            # 新 section
  revigo_cutoffs: [0.7, 0.5]   # REVIGO 语义相似度截断：先松后紧两轮
```

### 3.2 `workflow/src/clustering/candidates.py`

- `scale_features(data_df, selected_features, dr_cap=DR_CAP, dl_divisor=DL_DIVISOR)`：新增两个可覆盖参数，模块常量保留作默认值。
- 现有测试 `test_scale_features_caps_DR_and_divides_DL` 用默认值调用，行为不变，不用改。

### 3.3 `workflow/scripts/clustering/prepare_clustering_data.py`

- `PrepareConfig` 加 `dr_cap: float = DR_CAP`、`dl_divisor: int = DL_DIVISOR` 字段。
- `run()` 里调用 `scale_features(data_df, config.selected_features, dr_cap=config.dr_cap, dl_divisor=config.dl_divisor)`。
- argparse 加 `--dr-cap`（float, default `DR_CAP`）、`--dl-divisor`（int, default `DL_DIVISOR`）。
- `main()` 里把 `args.dr_cap` / `args.dl_divisor` 传进 `PrepareConfig`。

### 3.4 `workflow/rules/clustering.smk`（`prepare_clustering_data` rule）

```python
params:
    random_state=config.get("clustering", {}).get("random_state", 42),
    k_min=config.get("clustering", {}).get("k_min", 2),
    k_max=config.get("clustering", {}).get("k_max", 20),
    dr_cap=config.get("clustering", {}).get("dr_cap", 1.3),
    dl_divisor=config.get("clustering", {}).get("dl_divisor", 10),
...
shell:
    """
    python workflow/scripts/clustering/prepare_clustering_data.py \
        ... \
        --dr-cap {params.dr_cap} \
        --dl-divisor {params.dl_divisor} &> {log}
    """
```

### 3.5 `workflow/scripts/enrichment/run_network_enrichment.py`

- `NetworkConfig` 加 `revigo_cutoffs: list[float] = field(default_factory=lambda: list(REVIGO_CUTOFFS))`。
- `annotate_go_with_revigo` 签名加 `revigo_cutoffs: list[float]` 参数（由调用处传 `config.revigo_cutoffs`），内部 `for threshold in REVIGO_CUTOFFS` 改为 `for threshold in revigo_cutoffs`。
- argparse 加 `--revigo-cutoffs`（`nargs="+"`, `type=float`, default `REVIGO_CUTOFFS`）。
- `main()` 里把 `args.revigo_cutoffs` 传进 `NetworkConfig`。

### 3.6 `workflow/rules/enrichment_network.smk`

```python
params:
    ...
    revigo_cutoffs=" ".join(str(c) for c in config.get("enrichment_network", {}).get("revigo_cutoffs", [0.7, 0.5])),
...
shell:
    """
    python workflow/scripts/enrichment/run_network_enrichment.py \
        ... \
        --revigo-cutoffs {params.revigo_cutoffs} &> {log}
    """
```

Snakemake 的 shell 块里列表参数没有先例，这里用空格 join 成字符串，`argparse` 的 `nargs="+"` 天然能接收空格分隔的多个值。

## 4. 测试

- `test_clustering.py`：新增一个用例，验证 `scale_features` 传入非默认的 `dr_cap`/`dl_divisor` 时行为按新值变化（现有默认值用例不改）。
- `test_network_enrichment.py`：新增一个用例，验证 `annotate_go_with_revigo` 遍历的是传入的 `revigo_cutoffs` 列表而非硬编码的两个值（可以传单一阈值列表断言只跑一轮）。
- 跑 `pytest` 全量回归，确认现有 79+ 个测试仍然通过。
- `snakemake -n`（dry-run）确认两条 rule 的 shell 命令能正确插值。

## 5. 不做的事

- 不改 `_GAF_HEADER_DATE`、STRING/REVIGO 的 URL/身份标识/重试参数。
- 不改其他 section（`enrichment.*`、`ml.*`）已经在 config 里的字段。
- 不新增抽象层或通用的"常量注册机制"——两组常量各自按现有模式接入即可。

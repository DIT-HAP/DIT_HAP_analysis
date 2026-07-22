#!/usr/bin/env python3
"""
Automated migration script: pickle → parquet

Systematically replaces all pickle I/O with parquet I/O across the codebase.
"""

import re
from pathlib import Path
from typing import List, Tuple

# Files to migrate (DataFrame/Series only)
SCRIPT_FILES = [
    # Features
    "workflow/scripts/features/collect_dna_features.py",
    "workflow/scripts/features/collect_rna_features.py",
    "workflow/scripts/features/collect_protein_features.py",
    "workflow/scripts/features/collect_evolutionary_features.py",
    "workflow/scripts/features/collect_network_features.py",
    "workflow/scripts/features/collect_phenotype_features.py",
    "workflow/scripts/features/merge_features.py",
    # Clustering
    "workflow/scripts/clustering/prepare_clustering_data.py",
    "workflow/scripts/clustering/cluster_one_method.py",
    "workflow/scripts/clustering/finalize_direct_clusters.py",
    "workflow/scripts/clustering/finalize_auto_merge_clusters.py",
    "workflow/scripts/clustering/finalize_grid_clusters.py",
    # ML
    "workflow/scripts/ml/prepare_ml_data.py",
    "workflow/scripts/ml/train_automl.py",
]

RULE_FILES = [
    "workflow/rules/features.smk",
    "workflow/rules/clustering.smk",
    "workflow/rules/ml.smk",
]

def migrate_script(file_path: Path) -> List[str]:
    """Migrate one Python script from pickle to parquet."""
    content = file_path.read_text()
    changes = []

    # 1. Add import for write_parquet/read_parquet if not present
    if "from workflow.src.io import" in content and "write_parquet" not in content:
        # Find existing import line and extend it
        content = re.sub(
            r'(from workflow\.src\.io import [^\n]+)',
            r'\1, read_parquet, write_parquet',
            content
        )
        changes.append("Added parquet I/O imports")
    elif "from workflow.src.io import" not in content and ("to_pickle" in content or "read_pickle" in content):
        # Add new import after local imports section
        import_section = re.search(r'(# 4\. Local Imports.*?sys\.path\.insert.*?\n)', content, re.DOTALL)
        if import_section:
            insertion_point = import_section.end()
            content = content[:insertion_point] + "from workflow.src.io import read_parquet, write_parquet\n" + content[insertion_point:]
            changes.append("Added new parquet I/O import section")

    # 2. Replace df.to_pickle() with write_parquet()
    original = content
    content = re.sub(
        r'(\w+)\.to_pickle\(([^)]+)\)',
        r'write_parquet(\1, \2)',
        content
    )
    if content != original:
        changes.append("Replaced to_pickle() with write_parquet()")

    # 3. Replace pd.read_pickle() with read_parquet()
    original = content
    content = re.sub(
        r'pd\.read_pickle\(([^)]+)\)',
        r'read_parquet(\1)',
        content
    )
    if content != original:
        changes.append("Replaced pd.read_pickle() with read_parquet()")

    # 4. Update .pkl extensions in docstrings and comments to .parquet
    original = content
    content = re.sub(
        r'\.pkl(?=[\s:\)\],])',
        r'.parquet',
        content
    )
    if content != original:
        changes.append("Updated .pkl → .parquet in documentation")

    # 5. Update "pickle" mentions in docstrings
    original = content
    content = re.sub(
        r'(\s)pickle(\s)',
        r'\1parquet\2',
        content,
        flags=re.IGNORECASE
    )
    if content != original:
        changes.append("Updated 'pickle' → 'parquet' in text")

    file_path.write_text(content)
    return changes

def migrate_rule_file(file_path: Path) -> List[str]:
    """Migrate one Snakemake rule file from .pkl to .parquet extensions."""
    content = file_path.read_text()
    changes = []

    # Replace .pkl with .parquet in file paths
    original = content
    content = re.sub(
        r'\.pkl(?=["\'  \n])',
        r'.parquet',
        content
    )
    if content != original:
        changes.append("Updated .pkl → .parquet in rule file")

    file_path.write_text(content)
    return changes

def main():
    root = Path(__file__).parent

    print("=" * 70)
    print("PICKLE → PARQUET MIGRATION")
    print("=" * 70)
    print()

    # Migrate script files
    print("Migrating Python scripts...")
    for script_rel in SCRIPT_FILES:
        script_path = root / script_rel
        if not script_path.exists():
            print(f"  ⚠️  SKIP: {script_rel} (not found)")
            continue

        changes = migrate_script(script_path)
        if changes:
            print(f"  ✓ {script_rel}")
            for change in changes:
                print(f"      - {change}")
        else:
            print(f"  · {script_rel} (no changes needed)")

    print()

    # Migrate rule files
    print("Migrating Snakemake rule files...")
    for rule_rel in RULE_FILES:
        rule_path = root / rule_rel
        if not rule_path.exists():
            print(f"  ⚠️  SKIP: {rule_rel} (not found)")
            continue

        changes = migrate_rule_file(rule_path)
        if changes:
            print(f"  ✓ {rule_rel}")
            for change in changes:
                print(f"      - {change}")
        else:
            print(f"  · {rule_rel} (no changes needed)")

    print()
    print("=" * 70)
    print("Migration complete!")
    print()
    print("Next steps:")
    print("  1. Update workflow/src/features/assembly.py (read_coding_genes function)")
    print("  2. Update all test files in tests/")
    print("  3. Run: pytest tests/")
    print("  4. Run: snakemake -n")
    print("=" * 70)

if __name__ == "__main__":
    main()

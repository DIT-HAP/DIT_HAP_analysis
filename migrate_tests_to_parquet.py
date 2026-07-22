#!/usr/bin/env python3
"""
Automated migration script for test files: pickle → parquet
"""

import re
from pathlib import Path

TEST_FILES = [
    "tests/test_clustering.py",
    "tests/test_cluster_enrichment.py",
    "tests/test_train_automl.py",
    "tests/test_features_assembly.py",
]

def migrate_test_file(file_path: Path) -> list[str]:
    """Migrate one test file from pickle to parquet."""
    content = file_path.read_text()
    changes = []

    # 1. Replace .to_pickle() with write_parquet()
    # First check if write_parquet is imported
    if ".to_pickle(" in content or ".to_parquet(" in content:
        # Add import if not present
        if "from workflow.src.io import" not in content and "write_parquet" not in content:
            # Find imports section and add
            import_match = re.search(r'(import pytest.*?\n)', content, re.DOTALL)
            if import_match:
                insert_pos = import_match.end()
                content = content[:insert_pos] + "\nfrom workflow.src.io import write_parquet, read_parquet\n" + content[insert_pos:]
                changes.append("Added parquet I/O imports")

    # 2. Replace .to_pickle() calls with write_parquet()
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

    # 4. Update .pkl file extensions to .parquet
    original = content
    # In strings and paths
    content = re.sub(
        r'(["\'])([^"\']*?)\.pkl(["\'])',
        r'\1\2.parquet\3',
        content
    )
    if content != original:
        changes.append("Updated .pkl → .parquet in file paths")

    # 5. Update variable names
    original = content
    content = re.sub(r'\bpkl\b', 'parquet_file', content)
    content = re.sub(r'\b_pkl\b', '_parquet', content)
    if content != original:
        changes.append("Updated variable names pkl → parquet")

    file_path.write_text(content)
    return changes

def main():
    root = Path(__file__).parent

    print("=" * 70)
    print("TEST FILES: PICKLE → PARQUET MIGRATION")
    print("=" * 70)
    print()

    for test_rel in TEST_FILES:
        test_path = root / test_rel
        if not test_path.exists():
            print(f"  ⚠️  SKIP: {test_rel} (not found)")
            continue

        changes = migrate_test_file(test_path)
        if changes:
            print(f"  ✓ {test_rel}")
            for change in changes:
                print(f"      - {change}")
        else:
            print(f"  · {test_rel} (no changes needed)")

    print()
    print("=" * 70)
    print("Test migration complete!")
    print()
    print("Next: Run pytest to verify")
    print("=" * 70)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Fix missing parquet imports in scripts."""

import re
from pathlib import Path

files_to_fix = [
    "workflow/scripts/features/collect_evolutionary_features.py",
    "workflow/scripts/features/collect_network_features.py",
    "workflow/scripts/features/collect_phenotype_features.py",
    "workflow/scripts/features/collect_rna_features.py",
    "workflow/scripts/ml/prepare_ml_data.py",
    "workflow/scripts/ml/train_automl.py",
]

for file_path in files_to_fix:
    path = Path(file_path)
    if not path.exists():
        print(f"SKIP: {file_path} (not found)")
        continue
    
    content = path.read_text()
    
    # Check if already has the import
    if "from workflow.src.io import" in content and "parquet" in content:
        print(f"OK: {file_path} (already has import)")
        continue
    
    # Find the local imports section
    pattern = r'(# 3\. Local Imports.*?sys\.path\.insert.*?\n)(from workflow\.src\.)'
    match = re.search(pattern, content, re.DOTALL)
    
    if match:
        # Add the import after sys.path.insert
        insertion_point = match.end(1)
        new_content = (
            content[:insertion_point] +
            "from workflow.src.io import read_parquet, write_parquet\n" +
            content[insertion_point:]
        )
        path.write_text(new_content)
        print(f"FIXED: {file_path}")
    else:
        print(f"WARN: {file_path} (pattern not found)")

print("\nAll imports fixed!")

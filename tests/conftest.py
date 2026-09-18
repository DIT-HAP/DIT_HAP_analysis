"""Pytest configuration: make the repo importable the way the scripts do.

Modules under workflow/src/ import their siblings by bare name (`from io_table
import ...`, `from release_schema import ...`) — that is how every script imports
them, and it needs workflow/src itself on sys.path. Tests, meanwhile, import the
same code as `workflow.src.<domain>.<module>`, which needs the repo root. Both go
on the path here, once, instead of in each test file — three test modules were
silently un-collectable because they only had the repo root.
"""

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]

for _path in (_REPO_ROOT, _REPO_ROOT / "workflow" / "src"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

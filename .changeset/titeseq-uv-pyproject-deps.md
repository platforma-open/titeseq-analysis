---
'@platforma-open/platforma-open.titeseq-analysis.software': patch
---

Manage Python runtime dependencies in pyproject.toml as the single source of truth and generate src/requirements.txt from it (`pnpm deps:export`, top-level pins only). No runtime change: identical pinned versions, still built/installed via the pip toolset against the scientific-slim runenv.

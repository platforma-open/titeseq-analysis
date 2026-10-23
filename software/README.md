# Titeseq Analysis Software

The Python calculation script behind the titeseq-analysis block. Sources are in `src/`, tests in `tests/`.
Run every command from this directory.

## Development

```bash
uv sync                                     # install runtime and dev dependencies
uv run pytest -m "not slow"                 # fast test suite
uv run ruff check --fix && uv run ruff format
uv run ty check                             # type check; CI requires 0 diagnostics
```

Type every function signature in `src/`: all parameters and the return type, in plain 3.12 syntax
(`list[float]`, `X | None`). Leave local variables unannotated. Use `np.ndarray` for arrays and the polars types
(`pl.DataFrame`) callers already see. In tests, annotate helpers and fixtures; test functions need no
annotations. Add a `# ty: ignore[<rule>]` only for a gap in a third-party library's types, with a short reason.

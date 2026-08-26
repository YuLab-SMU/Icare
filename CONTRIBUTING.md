# Contributing to the icare tutorial

Changes to this book are reviewed through GitHub pull requests targeting the `docs` branch.

## Before opening a pull request

- Keep prose and code in English.
- Use small R chunks with a unique `#| label:` and place an explanation immediately before the chunk.
- Let important plots display directly below the chunk that creates them.
- Do not commit `_book/`, Quarto caches, frozen execution output, temporary R objects, or analysis output directories.
- Preserve a strict development/validation boundary in modeling examples.
- Do not change a package function name or claimed capability without checking the current `main` branch of `YuLab-SMU/icare`.
- Add new static manuscript figures under the matching `image/module*/` directory and prefer vector PDF output.
- Use informative alternative text or a clear figure caption for every reader-facing figure.

## Validation

Run the structural check:

```bash
quarto inspect
```

When the required R packages and data are available, run a full render:

```bash
quarto render
```

Review the GitHub Actions result before merging. A pull request should explain the affected chapter, any changed output, new dependencies, and whether numerical results or figures were regenerated.

## Scope

Tutorial-only changes belong in this branch. Changes to the R package, exported API, package documentation, or bundled R data should be proposed against the package's `main` branch instead.

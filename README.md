# icare Quarto Book

This directory contains the source of the English tutorial for [`YuLab-SMU/icare`](https://github.com/YuLab-SMU/icare), a modular R framework for clinical statistics, machine learning, unsupervised subtyping, and prognostic modeling.

The book is source-first: Quarto pages, tutorial data, branding assets, and publication-ready static figures are versioned; generated HTML and execution caches are not.

## Contents

- `index.qmd`: package overview, authors, scope, installation, and citation
- `prerequisites.qmd`: software setup, object flow, data contracts, and reproducibility guidance
- `module1-*.qmd`: `Stat` data cleaning and exploratory analysis
- `module2-*.qmd`: `Train_Model` classification modeling
- `module3-*.qmd`: `Subtyping` unsupervised discovery and validation
- `module4-*.qmd`: `PrognosiX` survival and prognostic modeling
- `data/`: tutorial data used by the survival workflow
- `.github/workflows/quarto-pages.yml`: PR build check and deployment from `docs`

## Local rendering

Install Quarto, R, and the package dependencies described in `prerequisites.qmd`, then run:

```bash
quarto render
```

The rendered site is written to `_book/`. That directory is ignored because GitHub Actions rebuilds it from source.

## Contribution workflow

1. Branch from the current `docs` branch.
2. Edit the smallest relevant `.qmd`, asset, data file, or workflow file.
3. Run `quarto inspect` and, when dependencies are available, `quarto render`.
4. Check that plots appear below their generating chunks and that all links and figure references resolve.
5. Open a pull request targeting `docs` and describe the scientific or instructional reason for the change.

For package behavior, function signatures, or defects, use the main [`icare` issue tracker](https://github.com/YuLab-SMU/icare/issues). See `CONTRIBUTING.md` for book-specific review checks.

## License and citation

`icare` is licensed under GPL-3. Please cite the package as described on the book home page and in the official repository.

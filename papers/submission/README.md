# Submission Package (arXiv + Zenodo)

This folder describes what to include for submission.

## arXiv
Include:
- `../null_results.tex` (standalone LaTeX manuscript)
- Figures referenced (if any) from `../../results/analysis/` and `../../results/reports/`
- `../references.bib` (optional; if switching to BibTeX)

Build locally to verify:
```bash
cd papers
pdflatex null_results.tex
pdflatex null_results.tex  # repeat as needed
```

## Zenodo (Data DOI)
Archive the following for reproducibility:
- `results/analysis/*.json` (raw sweep outputs)
- `results/analysis/*.{png,pdf}` (plots)
- `results/reports/analysis_report.md` and `analysis_tables.tex`
- `results/reports/csv/*.csv`

Recommended metadata:
- Title: "Null Results and Exclusion Limits for Coherence–Gravity and Curvature Couplings"
- Creators: Project Team
- Description: Short abstract + methods, link to GitHub repo
- Keywords: non-minimal coupling, curvature–EM, exclusion limits, tabletop gravity, null results

After DOI assignment, add the DOI badge to the repository README.

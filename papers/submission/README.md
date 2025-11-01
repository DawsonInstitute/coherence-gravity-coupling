# Submission Package (arXiv + Zenodo)

This folder describes what to include for submission.

## arXiv
Include:
- `../final_preprint.tex` (or `../null_results.tex`)
- `../null_results.md` (if using the markdown package)
- Figures referenced (if any) from `../../results/analysis/` and `../../results/reports/`
- `../references.bib` (if bibliography enabled)

Build locally to verify:
```bash
cd papers
pdflatex final_preprint.tex
# bibtex references  # uncomment if using bibliography
pdflatex final_preprint.tex
pdflatex final_preprint.tex
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

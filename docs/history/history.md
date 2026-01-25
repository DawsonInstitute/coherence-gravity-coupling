bibliography moved into null_results.bib and `null_results.tex` now uses BibTeX (`\bibliographystyle{unsrt}` + `\bibliography{null_results}`). 
<!-- ------ -->
## ✅ Update complete — removed unused entries

**What I changed**
- Removed `will2014` and `schlamminger2008` from null_results.bib.

**Verification**
- Confirmed null_results.tex does **not** cite those keys.  
- Found `will2014` used as an inline `\bibitem{will2014}` in other papers (e.g., extra_dof_mode_transitions.tex and curvature_em_to_bsm.tex). **`schlamminger2008` is not referenced elsewhere.**
<!-- ------ -->
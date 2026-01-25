bibliography moved into null_results.bib and `null_results.tex` now uses BibTeX (`\bibliographystyle{unsrt}` + `\bibliography{null_results}`). 
<!-- ------ -->
## ✅ Update complete — removed unused entries

**What I changed**
- Removed `will2014` and `schlamminger2008` from null_results.bib.

**Verification**
- Confirmed null_results.tex does **not** cite those keys.  
- Found `will2014` used as an inline `\bibitem{will2014}` in other papers (e.g., extra_dof_mode_transitions.tex and curvature_em_to_bsm.tex). **`schlamminger2008` is not referenced elsewhere.**
<!-- ------ -->
Done — I converted null_results.tex to REVTeX 4.2 and recompiled.

What I changed
- Switched document class to `\documentclass[aps,prb,preprint,amsmath,amssymb]{revtex4-2}`.
- Replaced the article-style title/author/footnote block with standard REVTeX front matter:
  - `\title{...}`, `\author{\authorname}`, `\affiliation{...}`, `\email{...}`, `\date{...}`
  - Moved `\maketitle` to appear after `\end{abstract}` (REVTeX requires the abstract before `\maketitle`).
- Switched bibliography style to `apsrev4-2`.
- Left `author_config.tex` in place so the same author/affiliation/email variables are still used.

Build result
- Recompiled successfully; output null_results.pdf created (11 pages).
- LaTeX warnings are non-critical (underfull/overfull hbox and hyperref token warnings). BibTeX ran and produced the bibliography.
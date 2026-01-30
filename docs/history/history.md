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
<!-- ------ -->
## ✅ Fixed — citation updated & paper recompiled 🔧

**What I changed**
- Replaced both occurrences of `eotwash2008` → **`schlamminger2008`** in:
  - coherence_gravity_coupling.tex (lines ~154 and ~178). ✅

**What I ran**
- Recompiled the paper (pdflatex → bibtex → pdflatex ×2). The bibliography now contains `\bibitem{schlamminger2008}` in `coherence_gravity_coupling.bbl` and coherence_gravity_coupling.pdf was regenerated.
<!-- ------ -->
---
Now I'll create the annotation file with a comprehensive summary of Verlinde's paper and its relevance to the manuscript.

Based on my analysis of the Verlinde (2011) paper and the manuscript citation context, here's my annotation:

**Created**: /home/echo_/Code/asciimath/coherence-gravity-coupling/papers/coherence_gravity_coupling-bib-annotations.md

The annotation is complete. The Verlinde citation is accurate and well-placed in the manuscript. The connection between emergent gravity and the $G_{\text{eff}}(\Phi)$ formalism is strong—both frameworks treat gravity as arising from information/entropy rather than as a fundamental force.
<!-- ------ -->
---
1.  **Annotated jacobson1995**: Read the converted text from jacobson1995.md and appended a detailed annotation to coherence_gravity_coupling-bib-annotations.md. The annotation highlights the thermodynamic derivation of the Einstein equation and explicitly connects the entropy area law to the manuscript's $G_{\text{eff}}(\Phi)$ formalism.
<!-- ------ -->
---
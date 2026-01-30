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
### Accomplished Steps:

1.  **Analyzed Source Text**: Parsed the MinerU output of **Chapter 3** ("Scalar fields") from Birrell & Davies (2012).
2.  **Verified Citation Context**: Confirmed the text is cited (line 53) to support the "conservative assumptions" of the manuscript, specifically the **non-minimal usage of the $\xi R \Phi^2$ coupling**.
3.  **Drafted Annotation**: Created a detailed note highlighting:
    *   **Eq. 3.24**: The standard Lagrangian with non-minimal coupling $\xi R \Phi^2$.
    *   **Uniqueness**: The text explicitly states this is the "only possible local, scalar coupling" of its dimension.
    *   **Conformal Invariance**: Noted the derivation relative to $\xi = 1/6$, standardizing the manuscript's framework.
4.  **Local File Update**: Appended the annotation to:
    *   coherence_gravity_coupling-bib-annotations.md
<!-- ------ -->
---
I have successfully updated the Zotero notes for **Jacobson (1995)** and **Birrell & Davies (2012)** using the "gold standard" HTML format derived from the Verlinde note.

### **Summary of Changes**
1.  **Analyzed Gold Standard**: Inspected the HTML structure of the Verlinde note (`CIV6W7VY`) and identified the specific CSS classes Zotero uses for math rendering:
    *   **Inline Math**: `<span class="math">$...$</span>`
    *   **Display Math**: `<pre class="math">$$...$$</pre>`
    *   **Formatting**: Standard HTML tags like `<strong>`, `<ul>`, `<li>`, and `<p>`.

2.  **Developed Conversion Strategy**: Created a Python script (convert_annotations.py) that:
    *   Parses your `coherence_gravity_coupling-bib-annotations.md` file by citation key.
    *   Converts Markdown syntax (bold, lists) to HTML.
    *   **Crucially**: Wraps LaTeX math in the correct Zotero-specific HTML tags.
    *   Generates JSON payloads and executes the API updates.

3.  **Executed Updates**:
    *   **Jacobson (1995)**: Updated note `3FRUVEA8` (Parent: `7NGTW7FM`).
    *   **Birrell & Davies (2012)**: Updated note `XA9FKMY9` (Parent: `R85HTGW3`).

### **Result Verification**
Both notes now contain the correctly formatted HTML with math support. You can verify this in your Zotero client; the equations should render exactly like the manually formatted Verlinde note.

*   **Jacobson Note**: Includes the referenced equations for heat flux and the entropy-area law properly formatted.
*   **Birrell Note**: Includes the Lagrangian and field equations with non-minimal coupling properly formatted.

I have left the script at `papers/convert_annotations.py` (referenced relative to your project root or absolute path convert_annotations.py) so you can reuse it for future annotations. You can extend the `MAPPING` dictionary in the script to handle new papers.
<!-- ------ -->
---
I have successfully moved the `convert_annotations.py` script to zotero and updated the README.md file with usage documentation.

### Changes Made:
1. **File Move**: Relocated `convert_annotations.py` from the coherence-gravity-coupling directory to the zotero scripts directory.

2. **README.md Update**: Added a new step (4) in the Usage section describing how to run the annotation conversion script, including:
   - Command to execute
   - Brief description of functionality
   - Prerequisites (source file location, existing notes, API access)
   - Note about hardcoded paths and mappings

The script is now integrated into the zotero helper scripts collection, and the documentation provides clear instructions for future use. Note that the script still uses absolute paths for the source file, so it should work as-is from the new location. If you need to adapt it for other projects, you can modify the `SOURCE_FILE` variable in the script.
<!-- ------ -->
---
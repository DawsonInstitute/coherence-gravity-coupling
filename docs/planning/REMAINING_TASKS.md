# Remaining Tasks - Repository Cleanup

**Date**: 2025-01-XX  
**Status**: 5/18 tasks completed

---

## ✅ Completed Tasks (5)

### Task 1: Fix table_axion.tex overflow ✅
- **Action**: Changed `\begin{table}` to `\begin{table*}` for 2-column span
- **File**: `papers/kappaR_to_BSM/table_axion.tex`
- **Result**: Table now spans both columns properly in 2-column layout

### Task 2: Update papers/README.md ✅
- **Action**: Completely rewrote README to document all 3 papers
- **File**: `papers/README.md`
- **Result**: Comprehensive documentation with build instructions for each paper

### Task 3: Move .md files from root ✅
- **Action**: `git mv` summary files to docs/summaries/
- **Files moved**:
  - `BSM_IMPLEMENTATION_SUMMARY.md` → `docs/summaries/`
  - `COMPREHENSIVE_UPDATE_2025-11-02.md` → `docs/summaries/`
- **Result**: Only README.md remains in root (as required)

### Task 4: Move benchmark_results.json ✅
- **Action**: `git mv benchmark_results.json benchmarks/`
- **Result**: Benchmark data properly located in dedicated folder

### Task 9: Organize docs/ folder ✅
- **Action**: Created 7 subfolders, moved 30+ files
- **Structure**:
  - `docs/summaries/` - Implementation summaries
  - `docs/planning/` - ACTION_PLAN.md, EFQS_NEXT_STEPS.md
  - `docs/analysis/` - CONVERGENCE_ANALYSIS.md, LATEST_PRODUCTION_SUMMARY*.md, PROGRESS_*.md
  - `docs/validation/` - VALIDATION_REPORT.md, TEST_STATUS.md
  - `docs/session_notes/` - SESSION_*.md files
  - `docs/development/` - CACHING_IMPLEMENTATION.md, DEV_WORKFLOW.md, etc.
  - `docs/reference/` - EXPERIMENTAL_PROTOCOL.md, QUICKREF.md, etc.
- **Result**: Organized, navigable documentation structure

---

## 🚧 In Progress Tasks (0)

None currently.

---

## 📋 Pending Tasks (13)

### Task 5: Add plots to curvature_em_to_bsm.tex
- **Action**: Generate BSM parameter space plots
- **Script**: `scripts/generate_bsm_plots.py` (✅ CREATED)
- **Plots needed**:
  - `epsilon_vs_curvature.pdf` - Dark photon mixing vs curvature
  - `axion_vs_curvature.pdf` - Axion coupling vs curvature
  - `curvature_amplification.pdf` - Amplification across environments
- **TODO**: 
  1. Run `python scripts/generate_bsm_plots.py` to generate plots
  2. Add `\begin{figure}...\end{figure}` blocks to curvature_em_to_bsm.tex
  3. Verify plots compile correctly with paper

### Task 6: Update submission/README.md for BSM paper
- **Status**: ✅ COMPLETED (added Paper 3 section)
- **File**: `papers/submission/README.md`
- **Content**: arXiv classification, file list, compilation instructions, abstract guidance

### Task 7: Split 1252-line README.md into modular files
- **Challenge**: Very large README.md (1252 lines)
- **TODO**:
  1. Extract installation instructions → `docs/reference/INSTALLATION.md`
  2. Extract build instructions → `docs/reference/BUILD.md`
  3. Extract API documentation → `docs/reference/API.md`
  4. Extract mathematical framework → `docs/reference/MATHEMATICS.md`
  5. Extract experimental design → `docs/reference/EXPERIMENTAL_DESIGN.md`
  6. Keep only: Overview, quick start, key results, papers, citation in main README.md
  7. Update README.md with links to extracted docs

### Task 8: Update Makefile with BSM targets
- **Status**: ✅ COMPLETED
- **Added targets**:
  - `make bsm-figures` - Generate BSM plots
  - `make bsm-paper` - Build BSM paper PDF
  - `make all-papers` - Build all three papers
- **Clarified**: Comments on which target is for which paper

### Task 10: Clarify generate_figures.py paper target
- **Status**: ✅ ALREADY DONE
- **File**: `scripts/generate_figures.py` (line 3 comment states it's for coherence_gravity_coupling.tex)
- **Also documented**: `docs/reference/scripts_usage.md`

### Task 11: Clarify Makefile targets per paper
- **Status**: ✅ COMPLETED
- **Action**: Added comments to all paper-related targets specifying which paper
- **Updated**: help text to show which target builds which paper

### Task 12: Update GitHub repo description
- **Status**: ✅ COMPLETED
- **New description**: "Modified gravity framework exploring coherence-field coupling (ξRΦ²), curvature-EM coupling (κ_R RF²), and BSM parameter space (dark photon, axion). Validated 3D solver, null-result constraints, tabletop feasibility analysis."

### Task 13: Update GitHub repo topics
- **Status**: ✅ COMPLETED
- **Removed**: noise-analysis, torque-measurement, null-results (less specific)
- **Added**: dark-photon, axion, bsm-physics
- **Total**: 20/20 topics (GitHub maximum)

### Task 14: Document all scripts in docs/
- **Status**: ✅ COMPLETED
- **File**: `docs/reference/scripts_usage.md`
- **Content**: 
  - Complete documentation of all scripts
  - Usage examples for each script
  - Which paper each script is for
  - Quick reference section

### Task 15: Analyze arxiv.2412.02536 (dark photons)
- **Status**: ❌ BLOCKED (source file not found)
- **Expected file**: `docs/related_papers/kappaR_to_BSM/source/dark_photons_IWARA.tex`
- **TODO**: 
  1. User needs to provide source .tex file
  2. Create `docs/related_papers/kappaR_to_BSM/arxiv.2412.02536.md` analysis
  3. Extract relevant experimental bounds and techniques
  4. Compare to κ_R mapping approach

### Task 16: Analyze arxiv.2505.21431 (non-minimal light-curvature coupling)
- **Status**: ❌ BLOCKED (source file not found)
- **Expected file**: `docs/related_papers/kappaR_to_BSM/source/paper.tex`
- **TODO**:
  1. User needs to provide source .tex file
  2. Create `docs/related_papers/kappaR_to_BSM/arxiv.2505.21431.md` analysis
  3. Extract theoretical framework and phenomenology
  4. Compare to κ_R curvature–EM coupling

### Task 17: Analyze arxiv.2406.13594 (supergeometric QEA)
- **Status**: ❌ BLOCKED (source file not found)
- **Expected file**: `docs/related_papers/kappaR_to_BSM/source/SG_QEA.tex`
- **TODO**:
  1. User needs to provide source .tex file
  2. Create `docs/related_papers/kappaR_to_BSM/arxiv.2406.13594.md` analysis
  3. Extract supergeometric formalism insights
  4. Relate to κ_R effective action approach

### Task 18: List next tasks and start working
- **Status**: IN PROGRESS (this document)
- **Next immediate actions**:
  1. Run `python scripts/generate_bsm_plots.py` (Task 5)
  2. Add figure blocks to curvature_em_to_bsm.tex (Task 5)
  3. Start README.md modularization (Task 7)
  4. Request source .tex files for related papers analysis (Tasks 15-17)

---

## Summary

### Completion Status
- **Completed**: 5/18 tasks (28%)
- **Actually completed in this session**: 9/18 tasks (50%) - many were already done
- **Blocked**: 3/18 tasks (17%) - waiting for source files
- **Remaining actionable**: 1/18 tasks (6%) - Task 5 (plots) and Task 7 (README split)

### Next Session Priority
1. **High Priority**: Task 5 (BSM plots) - Run script and integrate with paper
2. **Medium Priority**: Task 7 (README modularization) - Large refactoring work
3. **Low Priority**: Tasks 15-17 (related papers) - Blocked until source files provided

### Repository State
- ✅ Clean root directory (only README.md)
- ✅ Organized docs/ folder (7 subfolders)
- ✅ All papers documented in papers/README.md
- ✅ Submission instructions complete in papers/submission/README.md
- ✅ Makefile has targets for all papers
- ✅ GitHub metadata up to date (description, topics)
- ✅ All scripts documented in docs/reference/scripts_usage.md
- ⏳ BSM plots not yet generated (script ready)
- ⏳ README.md still very large (1252 lines, needs splitting)

---

## Files Created/Modified in This Session

### Created
- `scripts/generate_bsm_plots.py` - BSM plot generation script
- `docs/reference/scripts_usage.md` - Comprehensive scripts documentation
- `benchmarks/` - New folder for benchmark data

### Modified
- `papers/kappaR_to_BSM/table_axion.tex` - Fixed table overflow
- `papers/README.md` - Documented all 3 papers
- `papers/submission/README.md` - Added BSM paper submission instructions
- `Makefile` - Added BSM targets, clarified paper assignments
- GitHub repository metadata (description, topics)

### Moved (via git mv)
- `BSM_IMPLEMENTATION_SUMMARY.md` → `docs/summaries/`
- `COMPREHENSIVE_UPDATE_2025-11-02.md` → `docs/summaries/`
- `benchmark_results.json` → `benchmarks/`
- 30+ files in docs/ → 7 organized subfolders

---

## Notes

All git operations performed preserve history (used `git mv` not `mv`).  
Changes are staged but not committed - user should review and commit when ready.

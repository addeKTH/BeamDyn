# Railway bridge dynamics — web version

Browser version of `Main_gui.m` and `BeamDyn2D.m` (MATLAB version of 2026-11-01): multi-span 2D
Bernoulli-Euler beam, modal dynamics with the exact Nigam-Jennings recurrence, mode-acceleration
correction, HSLM-A (EN 1991-2:2024) and the EN 15528 line categories A1-E5, dynamic factor φ′ and
φ″ (EN 1991-2:2024 Annex C). Static files, no build step, no libraries.

## Running

Serve the folder over http (Web Workers do not start from `file://`), e.g.

```
python -m http.server 8000
```

and open http://localhost:8000, or publish the folder on GitHub Pages. Opened from `file://` the page
falls back to a serial run in the page itself (slower, Stop still works).

## Interface

Same layout and workflow as the MATLAB figure (560 × 908 px):

- Symbol bar: Run (Stop after 5 s), Reset, Open and Save of the input (JSON, same fields as MATLAB),
  Export (xlsx), Print (report), progress bar, status text.
- Input: spans, mass, EI; damping from EN 1991-2 Table 8.6 (menu, L = L_Φ) or typed; N_mod and f_max
  automatic per EN 1990:2023 A.2.9.4.2.1(3) (grey) or typed; dL; train (HSLM-A, A1-E5, All); cars;
  speeds (`100:1:200`); load distribution b_Q (1-, 3-, 5-point or a typed triangle base, default
  triangular 3.0 m).
- Model: beam, supports (click to clamp), mode shapes four at a time, legend on the right.
- Results: a_max, d_max or φ_max, M_max or M_min, V_max; φ′ and φ″ tick boxes, A1-A10 separately,
  Cumulative, Clear, Print of the four plots. Click a curve for a datatip (drag along it,
  double-click to remove).

`window.inp`, `window.mdl` and `window.out` hold the input, model and results (console).

## Files

| File | Content |
|---|---|
| `index.html` | layout and styles |
| `app.js` | GUI logic (port of `Main_gui.m`): input, model update, runs, plots, envelope, φ″, Open/Save, Export, Print |
| `beamdyn2d.js` | solver (port of `BeamDyn2D.m`): model, eigenvalues, trains, load distribution, static part with exact element-end forces, speed loop, φ′, Eurocode helpers |
| `linalg.js` | band storage, banded Cholesky, dense symmetric eigensolver (tred2/tql2), IIR filter |
| `plot.js` | SVG plots (ticks, legend, datatips) and the LaTeX subset for labels |
| `solver.worker.js` | one train (or a group of speeds) per worker job |
| `xlsx.js` | minimal xlsx writer (from the earlier web version) |

## Differences to MATLAB

- Export writes xlsx only (no .mat) and without native Excel charts.
- Print and the plot print open the browser print dialog (choose "Save as PDF"); the report is HTML
  with vector SVG figures in the layout of the LaTeX report. Date and page number "N (M)" use CSS page
  margin boxes (Chrome/Edge 131 or newer).
- Eigenvalues: dense up to 600 free dofs (whole spectrum, as `eig`); above that subspace iteration
  for the modes below f_max (Sturm count) instead of `eigs`.
- Mode-shape signs may differ from MATLAB; results do not depend on them.

## Validation

`../../test/compare_web.html` (served from `BeamDyn/`) runs four cases against
`../../test/run_web_cases.m` (BeamDyn2D.m): frequencies and all dynamic and static peaks agree within
1e-8 relative, φ′ within 1e-7 (HSLM-A, clamped supports, 3-/5-point and triangular load distribution, a
722-dof model on the subspace path). Exception: φ′ of M_min on a single simply supported span, where the
static M_min is round-off zero on both sides.

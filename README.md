# BeamDyn2D — web version

Browser port of the MATLAB program. Dynamic response of a multi-span
Bernoulli–Euler beam to passing trains, swept over a speed range.

Everything runs in the browser. No server, no build step, no dependencies,
no licence.

## Files

| | |
|---|---|
| `index.html` | the page: UI, SVG plots, print-to-PDF report |
| `beamdyn.js` | numerics core — works in the browser *and* in Node |
| `xlsx.js` | minimal .xlsx writer with scatter charts, no dependencies |
| `test.js` | validation suite (Node) |

## Publish on GitHub Pages

Drop `index.html` and `beamdyn.js` in a repo, Settings → Pages → deploy from
branch. That's all — there is nothing to build.

Locally: `python3 -m http.server` then open `localhost:8000`.
(Opening `index.html` directly with `file://` works too.)

## Validate

```
node test.js
```

Ten checks against analytical solutions: simply-supported, clamped–clamped
and two-span continuous frequencies; mass normalisation; Nigam–Jennings
against the exact ramp response; static and velocity gains; the
consistent-load / direct-φ identity; the static deflection limit.

## Excel export

**Save** writes a workbook with three sheets — mode shapes, displacement
envelope, acceleration envelope — each with a scatter chart (straight lines,
no markers) beside its data.

The envelopes carry one column per accumulated run, over the union of all
speeds; a run with no sample at a given speed leaves a blank, which the chart
draws as a gap.

An HSLM run writes **eleven** columns: `HSLM-A1` … `HSLM-A10` with each train
separately, then `HSLM envelope` — the maximum over the ten. Only the envelope
is plotted; the per-train columns are there for inspection.

Mode shapes are those of the *current* model, i.e. the ones on screen, scaled
to max|φ| = 1. Note that this scaling means they are **not** mass-normalised
and cannot be used to reconstruct response amplitudes — use `out.Phi` from
`beamdyn.js` for that. All values are rounded to 3 decimals; frequencies in
the column headers to 2.

**Save as** uses the browser's native file dialog where it exists
(Chrome/Edge/Opera, via the File System Access API — needs `https` or
`localhost`, so GitHub Pages qualifies). Firefox and Safari have no such API:
there it falls back to an ordinary download, which still shows a Save-as
dialog if the browser is set to *always ask where to save files*.

`xlsx.js` builds the file from scratch: a store-only ZIP plus hand-written
SpreadsheetML and DrawingML. No library, and unlike the free build of SheetJS
it can emit charts. Verified by reading the output back with openpyxl.

## Method

Identical to `BeamDyn2D.m`. See `BeamDyn2D_theory.tex` for the derivations.

- 2-node Bernoulli–Euler elements, consistent mass, dense assembly.
- Generalised eigenproblem reduced to standard symmetric form via a Cholesky
  factorisation of **M**, then solved by cyclic Jacobi. Mass normalisation
  falls out of the reduction.
- Modal load built in *position* space, so it is speed-independent; speed
  enters only through dt = ds/v.
- Nigam–Jennings exact recurrence for a piecewise-linear load, applied as a
  two-pole IIR filter. Modal acceleration straight from the ODE.

## Performance

Measured in Node, 2×20 m bridge, 3 modes, 101 speeds:

| | |
|---|---|
| D2, 20 cars | ~0.4 s |
| HSLM-A, 10 trains | ~4 s |

The run yields to the event loop between trains, so the page stays
responsive and the progress label updates.

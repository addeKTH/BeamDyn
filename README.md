# Railway bridge dynamics — web version

A browser port of `Main_gui.m` and the solvers behind it. Multi-span 2D
Bernoulli-Euler beams, modal dynamics, Nigam-Jennings time integration, HSLM-A
and the EN 15528 line categories. No build step, no dependencies, no server-side
anything — static files you can drop straight into a GitHub Pages repository.

## Running it

The page uses Web Workers, and browsers refuse to create those from `file://`.
Serve the folder:

```
python -m http.server 8000      # then open http://localhost:8000
```

or push it to GitHub Pages. `Save` uses the File System Access API for a proper
Save-as dialog in Chrome and Edge, and falls back to a plain download elsewhere.
Opening `index.html` directly still works — it
detects the missing workers and runs everything on the main thread instead, so
it is correct but serial and the interface freezes during a run. The `par`
button reads `no par` in that case.

## The interface

Same layout as the MATLAB figure, 540 × 680 px, same workflow.

| Control | Behaviour |
|---|---|
| `L, m, EI, ζ` | bridge properties; `L` takes a list, e.g. `24 30 24` |
| `Nmod, fmax, dL` | modes kept, frequency cut-off, element size |
| `load` | HSLM-A (all ten Universal Trains) or a line category A1…E5 |
| `cars` | wagons in the train set; disabled for HSLM-A, which defines its own |
| `v` | speeds in km/h; MATLAB range syntax works, `100:1:200` |
| `W` | width of the triangular load distribution; 0 = point loads |
| `par` | one worker or every core — this is what replaces `parpool` |
| `Run` | dynamic + static analysis of the selected model |
| `Save` | xlsx workbook: Input, ModeShapes, Results — opens a Save-as dialog |
| `Clear` / `Reset` | drop the accumulated runs / restore the defaults |
| `Report` | builds an A4 report and opens the print dialog — "Save as PDF" |
| Section-forces | compute M, V, reactions and φ′; off is roughly twice as fast |
| Cumulative | running peak over speed |
| Envelope | on Run: HSLM-A against all twelve categories |
| Separate A1-A10 | draw the ten HSLM-A trains individually instead of their envelope |

Click a support in the beam sketch to switch it between pinned and clamped; the
model, frequencies and mode shapes update immediately. Click a result curve for
a datatip, drag it along the curve, double-click to remove it.

The results object is on `window.out`, and `window.inp` / `window.mdl` hold the
input and the model — the equivalent of the `assignin('base', ...)` calls in the
m-file. `setSupport(i, true)` clamps support `i` from the console.

## Files

| File | MATLAB counterpart |
|---|---|
| `index.html` | `Main_gui.m` — layout, callbacks, plotting, export |
| `beamdyn2d.js` | `BeamDyn2D.m` |
| `beamstat2d.js` | `BeamStat2D.m` |
| `dynamicfactor.js` | `DynamicFactor.m` |
| `trainloads.js` | `TrainLoadHSR.m` (HSLMA) + the EN 15528 branch of `TrainLoadSOU.m` |
| `beammodel.js` | the model block both solvers share — assembly, b.c., recovery |
| `linalg.js` | banded Cholesky, subspace iteration, the `filter` recurrence |
| `plot.js` | `PlotBeam2D.m` plus the axes work MATLAB does for free |
| `solver.worker.js` | the `parfor` body |
| `xlsx.js` | `writecell` / `writematrix` |

## What changed, and why

**Banded everything.** K and M have half bandwidth 3, so the factorisation is
O(n·bw²) and each solve O(n·bw). The eigenvalue problem uses subspace iteration
against that factorisation and returns only the modes asked for, instead of a
dense `eig` that computes the whole spectrum. `updateModel` stays interactive
down to dL = 0.05 m.

**Peaks without the response.** Only extremes leave the solver, so the recovery
`Psi*q` is fused with the max/min reduction — the [output points × time steps]
array is never formed. That is the difference between a few hundred kB and a few
hundred MB in envelope mode, and it made the inner loop about 2.6× faster.

**The shear is evaluated once per element**, not at both ends: `w'''` is constant
inside an element, so the two end values were the same number.

**ds and the record length come from the full speed list**, not from the subset a
worker happens to get. Splitting a run over eight workers therefore gives
bit-identical results to a serial one — there is a test for it.

Rough timings on a modern 8-core laptop, 2 × 20 m, dL = 1, v = 100:1:200, section
forces on: one line category ≈ 0.8 s serial; HSLM-A (ten trains) ≈ 10 s serial and
under 2 s across the workers; a full envelope ≈ 25 s serial, 4-5 s parallel.

## Validation

`node test.js` (no dependencies) checks the engine against the MATLAB reference
case L = 25 m pinned, m = 15e3, EI = 30e9, ζ = 1 %, dL = 1.25, D2 with 6 wagons,
Nmod = 3:

- natural frequencies 3.554308 / 14.217321 / 31.989846 Hz — match to the printed
  precision of the reference, and mode 1 matches the closed form for a simply
  supported beam
- static `dmax` = 1.123341e-2 m and `Mmax` = 5.077227e6 Nm — match
- `fs` = 640 Hz and `t_fr` = 0.06 s — match
- HSLM-A comes out 375.5 to 398.5 m long and 7140 to 9600 kN, exactly the range
  EN 1991-2 states for the ten Universal Trains
- a train at 10 km/h gives φ′ ≈ 0.014 against the static solver, where EN 1991-2
  Annex C gives 0.016 at that K — the two solvers agree independently
- the static correction shifts `dmax` by 0.03 %, as the m-file documents

`node test_gui.js` and `node test_gui2.js` load `index.html` in jsdom
(`npm i jsdom` first) and drive the actual interface: input parsing, the model
update, a full Run, envelope mode, clamping a support, the load distribution, the
report markup, and the xlsx file, which is unzipped and checked.

## Known gaps against the MATLAB version

- **No `.mat` file.** `Save` writes the xlsx; `window.out` covers the rest.
- **No native Excel charts.** The MATLAB version drove Excel through COM to build
  scatter charts inside the workbook, which a browser cannot do. The data and the
  peak summary are all there, ready to chart in two clicks.
- **The report is a print stylesheet**, not a vector PDF writer: the tables are
  real text, the four plots are 3× PNGs. Print to PDF from the browser. It
  carries two peak tables, one over the line categories and one over the ten
  HSLM-A trains, and quotes EI in GNm².
- Mode shape signs may be flipped relative to MATLAB — `eig` fixes no sign
  either, and nothing downstream depends on it.
- The influence-line stage of `beamstat2d.js` is O(nNode²) in memory, like the
  MATLAB `conv2` version. It is fine to dL ≈ 0.1 m on a 40 m bridge and slow
  below that; the dynamic path has no such limit.

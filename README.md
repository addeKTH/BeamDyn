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

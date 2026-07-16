/*--------------------------------------------------------------------------
  BeamDyn2D - browser port of the MATLAB program
  Numerics core: no dependencies, works in Node and in the browser.

  Andreas Andersson, 2026-07-16
  assisted by Claude Opus 4.8

  Dense matrices are row-major Float64Array: A[i*n+j].
  Mode shape arrays Phi are column-major: Phi[n*ndof + i]  (mode n, dof i).
--------------------------------------------------------------------------*/
(function (root) {
'use strict';

/* ===================== LINEAR ALGEBRA ================================== */

/** Cholesky M = L*L', M symmetric positive definite. Returns lower L. */
function cholesky(M, n) {
  const L = new Float64Array(n * n);
  for (let j = 0; j < n; j++) {
    let s = M[j * n + j];
    for (let k = 0; k < j; k++) s -= L[j * n + k] * L[j * n + k];
    if (s <= 0) throw new Error('mass matrix is not positive definite');
    L[j * n + j] = Math.sqrt(s);
    const ljj = L[j * n + j];
    for (let i = j + 1; i < n; i++) {
      let t = M[i * n + j];
      for (let k = 0; k < j; k++) t -= L[i * n + k] * L[j * n + k];
      L[i * n + j] = t / ljj;
    }
  }
  return L;
}

/** Solve L*X = B for X, L lower triangular, B and X n-by-n. */
function forwardSolveMat(L, B, n) {
  const X = new Float64Array(n * n);
  for (let c = 0; c < n; c++) {
    for (let i = 0; i < n; i++) {
      let s = B[i * n + c];
      for (let k = 0; k < i; k++) s -= L[i * n + k] * X[k * n + c];
      X[i * n + c] = s / L[i * n + i];
    }
  }
  return X;
}

/** Solve L'*x = b for x, L lower triangular (so L' is upper). */
function backSolveVec(L, b, n) {
  const x = new Float64Array(n);
  for (let i = n - 1; i >= 0; i--) {
    let s = b[i];
    for (let k = i + 1; k < n; k++) s -= L[k * n + i] * x[k];
    x[i] = s / L[i * n + i];
  }
  return x;
}

function transpose(A, n) {
  const T = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) T[j * n + i] = A[i * n + j];
  return T;
}

/**
 * Cyclic Jacobi eigenvalue decomposition of a symmetric matrix.
 * Returns { lam, V } with V column-major: eigenvector k is V[k*n .. k*n+n-1].
 * Chosen over Householder+QL because it is short, needs no extra storage,
 * and delivers orthonormal eigenvectors to machine precision.
 */
function jacobiEig(Ain, n) {
  const A = Float64Array.from(Ain);
  const V = new Float64Array(n * n);          // row-major during rotation
  for (let i = 0; i < n; i++) V[i * n + i] = 1;

  let scale = 0;
  for (let i = 0; i < n; i++) scale += A[i * n + i] * A[i * n + i];
  const tol = 1e-22 * Math.max(scale, 1e-300);

  for (let sweep = 0; sweep < 60; sweep++) {
    let off = 0;
    for (let i = 0; i < n; i++)
      for (let j = i + 1; j < n; j++) off += A[i * n + j] * A[i * n + j];
    if (off <= tol) break;

    for (let p = 0; p < n - 1; p++) {
      for (let q = p + 1; q < n; q++) {
        const apq = A[p * n + q];
        if (apq === 0) continue;
        const app = A[p * n + p], aqq = A[q * n + q];
        // skip rotations that cannot change anything at this precision
        if (Math.abs(apq) < 1e-18 * Math.sqrt(Math.abs(app * aqq))) continue;
        const theta = (aqq - app) / (2 * apq);
        const sgn = theta >= 0 ? 1 : -1;
        const t = sgn / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;

        for (let k = 0; k < n; k++) {              // columns p,q
          const akp = A[k * n + p], akq = A[k * n + q];
          A[k * n + p] = c * akp - s * akq;
          A[k * n + q] = s * akp + c * akq;
        }
        for (let k = 0; k < n; k++) {              // rows p,q
          const apk = A[p * n + k], aqk = A[q * n + k];
          A[p * n + k] = c * apk - s * aqk;
          A[q * n + k] = s * apk + c * aqk;
        }
        for (let k = 0; k < n; k++) {              // accumulate V
          const vkp = V[k * n + p], vkq = V[k * n + q];
          V[k * n + p] = c * vkp - s * vkq;
          V[k * n + q] = s * vkp + c * vkq;
        }
      }
    }
  }
  const lam = new Float64Array(n);
  for (let i = 0; i < n; i++) lam[i] = A[i * n + i];
  // repack V to column-major (eigenvector k contiguous)
  const Vc = new Float64Array(n * n);
  for (let k = 0; k < n; k++) for (let i = 0; i < n; i++) Vc[k * n + i] = V[i * n + k];
  return { lam, V: Vc };
}

/* ===================== FINITE ELEMENT MODEL ============================ */

/**
 * Assemble K and M on the free dofs and solve the eigenproblem.
 *   K*phi = w^2*M*phi   is reduced to a standard symmetric problem via
 *   M = L*L',  A = inv(L)*K*inv(L'),  phi = inv(L')*phiTilde.
 * Because phiTilde is orthonormal, phi' * M * phi = I falls out of the
 * reduction: the mass normalisation needs no separate step.
 */
function buildModel(inp) {
  const dL = inp.dL;
  const L = inp.L.map(x => dL * Math.round(x / dL));
  const Ltot = L.reduce((a, b) => a + b, 0);
  const EL = Math.round(Ltot / dL);
  const ndof = 2 * (EL + 1);

  // element matrices
  const k1 = 12 * inp.EI / (dL * dL * dL);
  const k2 = 6 * inp.EI / (dL * dL);
  const k3 = 4 * inp.EI / dL;
  const ke = [
    k1, k2, -k1, k2,
    k2, k3, -k2, k3 / 2,
    -k1, -k2, k1, -k2,
    k2, k3 / 2, -k2, k3];
  const c = inp.m * dL / 420, d = dL, d2 = dL * dL;
  const me = [
    156 * c, 22 * d * c, 54 * c, -13 * d * c,
    22 * d * c, 4 * d2 * c, 13 * d * c, -3 * d2 * c,
    54 * c, 13 * d * c, 156 * c, -22 * d * c,
    -13 * d * c, -3 * d2 * c, -22 * d * c, 4 * d2 * c];

  // assembly (dense: ndof is small enough that sparsity buys nothing here)
  const K = new Float64Array(ndof * ndof);
  const M = new Float64Array(ndof * ndof);
  for (let e = 0; e < EL; e++) {
    const g = 2 * e;                                  // first dof of element e
    for (let a = 0; a < 4; a++) {
      for (let b = 0; b < 4; b++) {
        K[(g + a) * ndof + (g + b)] += ke[a * 4 + b];
        M[(g + a) * ndof + (g + b)] += me[a * 4 + b];
      }
    }
  }

  // boundary conditions: dof = 2*(x/dL) + (dir-1)
  const fixed = new Set();
  for (const [x, dir] of inp.bc) {
    if (x < -1e-9 || x > Ltot + 1e-9) throw new Error('support outside the beam');
    if (dir !== 1 && dir !== 2) throw new Error('dir must be 1 (vertical) or 2 (rotation)');
    const nd = x / dL;
    if (Math.abs(nd - Math.round(nd)) > 1e-6)
      throw new Error('support at x = ' + x + ' is not on a node (dL = ' + dL + ')');
    fixed.add(2 * Math.round(nd) + (dir - 1));
  }
  const fdof = [];
  for (let i = 0; i < ndof; i++) if (!fixed.has(i)) fdof.push(i);
  const nf = fdof.length;

  // extract the free blocks
  const Kf = new Float64Array(nf * nf), Mf = new Float64Array(nf * nf);
  for (let i = 0; i < nf; i++) {
    for (let j = 0; j < nf; j++) {
      Kf[i * nf + j] = K[fdof[i] * ndof + fdof[j]];
      Mf[i * nf + j] = M[fdof[i] * ndof + fdof[j]];
    }
  }

  // reduce to standard symmetric form
  const Lc = cholesky(Mf, nf);
  const Y = forwardSolveMat(Lc, Kf, nf);              // inv(L)*K
  const Z = forwardSolveMat(Lc, transpose(Y, nf), nf); // inv(L)*(inv(L)*K)' = A'
  for (let i = 0; i < nf; i++)                        // symmetrise
    for (let j = i + 1; j < nf; j++) {
      const v = 0.5 * (Z[i * nf + j] + Z[j * nf + i]);
      Z[i * nf + j] = v; Z[j * nf + i] = v;
    }

  const { lam, V } = jacobiEig(Z, nf);

  // sort ascending
  const order = Array.from({ length: nf }, (_, i) => i).sort((a, b) => lam[a] - lam[b]);

  // retain the requested modes below fmax, but never fewer than one
  const nreq = Math.min(inp.Nmod, nf);
  const fAll = order.map(i => Math.sqrt(Math.max(lam[i], 0)) / (2 * Math.PI));
  let keep = [];
  for (let i = 0; i < nreq; i++) if (fAll[i] < inp.fmax) keep.push(i);
  if (keep.length === 0) keep = [0];
  const nm = keep.length;

  // back-transform and scatter into the full dof set
  const f = new Float64Array(nm);
  const Phi = new Float64Array(nm * ndof);
  for (let n = 0; n < nm; n++) {
    const src = order[keep[n]];
    f[n] = fAll[keep[n]];
    const pt = V.subarray(src * nf, src * nf + nf);
    const p = backSolveVec(Lc, pt, nf);              // phi = inv(L')*phiTilde
    for (let i = 0; i < nf; i++) Phi[n * ndof + fdof[i]] = p[i];
  }

  // explicit mass normalisation: cheap, and independent of the reduction
  for (let n = 0; n < nm; n++) {
    const b = n * ndof;
    let mm = 0;
    for (let i = 0; i < ndof; i++) {
      let s = 0;
      for (let j = 0; j < ndof; j++) s += M[i * ndof + j] * Phi[b + j];
      mm += Phi[b + i] * s;
    }
    const sc = 1 / Math.sqrt(mm);
    for (let i = 0; i < ndof; i++) Phi[b + i] *= sc;
  }

  return { L, Ltot, EL, ndof, dL, f, Phi, nm, K, M, fdof };
}

/* ===================== MODE SHAPES ===================================== */

/** Cubic Hermite shape functions of the 2-node beam element. */
function shapeFun(s, L) {
  const s2 = s * s, s3 = s * s2;
  return [1 - 3 * s2 + 2 * s3, L * (s - 2 * s2 + s3), 3 * s2 - 2 * s3, L * (s3 - s2)];
}

/** Mode shapes at one point x. Returns Float64Array of length nm. */
function hermiteEval(Phi, nm, ndof, dL, EL, x) {
  let e = Math.floor(x / dL);
  if (e < 0) e = 0; else if (e > EL - 1) e = EL - 1;
  const N = shapeFun(x / dL - e, dL);
  const i1 = 2 * e;
  const out = new Float64Array(nm);
  for (let n = 0; n < nm; n++) {
    const b = n * ndof;
    out[n] = N[0] * Phi[b + i1] + N[1] * Phi[b + i1 + 1] +
             N[2] * Phi[b + i1 + 2] + N[3] * Phi[b + i1 + 3];
  }
  return out;
}

/** Mode shapes at many points. Returns [nm x nx] row-major. */
function hermiteEvalMany(Phi, nm, ndof, dL, EL, xs) {
  const nx = xs.length;
  const P = new Float64Array(nm * nx);
  for (let j = 0; j < nx; j++) {
    const v = hermiteEval(Phi, nm, ndof, dL, EL, xs[j]);
    for (let n = 0; n < nm; n++) P[n * nx + j] = v[n];
  }
  return P;
}

/* ===================== TIME INTEGRATION =============================== */

/**
 * Nigam-Jennings exact recurrence for a unit-mass oscillator
 *     qdd + 2*z*w*qd + w^2*q = p(t),  p linear over each step dt,
 * returned as digital filter coefficients sharing the denominator a.
 */
function njCoeff(w, z, dt) {
  const sq = Math.sqrt(1 - z * z);
  const wd = w * sq;
  const e = Math.exp(-z * w * dt);
  const S = Math.sin(wd * dt);
  const Co = Math.cos(wd * dt);
  const zz = z / sq;
  const k = w * w;                                  // unit modal mass

  const A = e * (zz * S + Co);
  const B = e * (S / wd);
  const C = (2 * z / (w * dt) + e * (((1 - 2 * z * z) / (wd * dt) - zz) * S
            - (1 + 2 * z / (w * dt)) * Co)) / k;
  const D = (1 - 2 * z / (w * dt) + e * (((2 * z * z - 1) / (wd * dt)) * S
            + (2 * z / (w * dt)) * Co)) / k;

  const Ap = -e * (w / sq) * S;
  const Bp = e * (Co - zz * S);
  const Cp = (-1 / dt + e * ((w / sq + zz / dt) * S + Co / dt)) / k;
  const Dp = (1 - e * (zz * S + Co)) / (k * dt);

  return {
    a: [1, -(A + Bp), A * Bp - Ap * B],
    bq: [D, C - Bp * D + B * Dp, B * Cp - Bp * C],
    bv: [Dp, Cp + Ap * D - A * Dp, Ap * C - A * Cp]
  };
}

/** Direct-form-I two-pole IIR, y = filter(b, a, x). */
function filter2(b, a, x, y) {
  const n = x.length;
  if (!y) y = new Float64Array(n);
  let x1 = 0, x2 = 0, y1 = 0, y2 = 0;
  const b0 = b[0], b1 = b[1], b2 = b[2], a1 = a[1], a2 = a[2];
  for (let i = 0; i < n; i++) {
    const xi = x[i];
    const yi = b0 * xi + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
    y[i] = yi;
    x2 = x1; x1 = xi; y2 = y1; y1 = yi;
  }
  return y;
}

/* ===================== TRAIN LOADS ==================================== */

const LINE_CAT = {
  name: ['A', 'B1', 'B2', 'C2', 'C3', 'C4', 'D2', 'D3', 'D4', 'D5', 'E4', 'E5'],
  c: [6.2, 7.8, 4.65, 5.9, 4.5, 3.4, 7.45, 5.9, 4.65, 3.62, 5.9, 4.75],
  P: [160, 180, 180, 200, 200, 200, 225, 225, 225, 225, 250, 250].map(v => v * 1e3),
  a: 1.8, b: 1.5, N_axle: 4
};

/** triang(N) without the Signal Processing Toolbox. */
function triangWin(N) {
  const w = new Float64Array(N);
  for (let i = 1; i <= N; i++) {
    w[i - 1] = (N % 2)
      ? Math.min(2 * i / (N + 1), 2 * (N + 1 - i) / (N + 1))
      : Math.min((2 * i - 1) / N, (2 * (N + 1 - i) - 1) / N);
  }
  return w;
}

/** Axle positions and loads: [[pos, P], ...], first axle at pos = 0. */
function trainLoad(inp, nr) {
  let sp, P;
  if (inp.load === 'HSLM') {
    const N = [18, 17, 16, 15, 14, 13, 13, 12, 11, 11];
    const D = [18, 19, 20, 21, 22, 23, 24, 25, 26, 27];
    const d = [2, 3.5, 2, 3, 2, 2, 2, 2.5, 2, 2];
    const Pa = [170, 200, 180, 190, 170, 180, 190, 190, 210, 210].map(v => v * 1e3);
    const dn = d[nr], Dn = D[nr];
    sp = [0, 3, 11, 3, 3.525, dn, Dn - 1.5 * dn - 1.7625];
    for (let i = 0; i < N[nr]; i++) sp.push(dn, Dn - dn);
    sp.push(dn, Dn - 1.5 * dn - 1.7625, dn, 3.525, 3, 11, 3);
    P = Pa[nr];
  } else {
    const k = LINE_CAT.name.indexOf(inp.load);
    if (k < 0) throw new Error('unknown load model "' + inp.load + '"');
    sp = [];
    for (let i = 0; i < inp.N_car; i++)
      sp.push(2 * LINE_CAT.b, LINE_CAT.a, LINE_CAT.c[k], LINE_CAT.a);
    P = LINE_CAT.P[k];
  }
  const F = [];
  let x = 0;
  for (let i = 0; i < sp.length; i++) { x += sp[i]; F.push([x, P]); }
  const x0 = F[0][0];
  for (const r of F) r[0] -= x0;                    // first axle at 0
  return F;
}

/** Triangular load distribution over width W, then merge coincident loads. */
function distribute(F, W, dL) {
  if (W > 0) {
    const nW = Math.floor(W / dL) + 1;
    const xs = [];
    for (let i = 0; i < nW; i++) xs.push(-W / 2 + i * dL);
    const t = triangWin(nW);
    let ts = 0; for (const v of t) ts += v;
    const G = [];
    for (const [x, P] of F)
      for (let i = 0; i < nW; i++) G.push([x + xs[i], P * t[i] / ts]);
    let mn = Infinity; for (const r of G) mn = Math.min(mn, r[0]);
    for (const r of G) r[0] -= mn;
    F = G;
  }
  const map = new Map();                            // merge, like accumarray
  for (const [x, P] of F) {
    const key = Math.round(x * 1e6) / 1e6;
    map.set(key, (map.get(key) || 0) + P);
  }
  return Array.from(map.entries()).sort((a, b) => a[0] - b[0]);
}

/* ===================== SOLVER ========================================= */

/** Linear interpolation of a [nm x ns] row-major array onto new abscissae. */
function resample(F, nm, s0, ds, ns, sq, nq) {
  const out = new Float64Array(nm * nq);
  for (let j = 0; j < nq; j++) {
    const t = (sq[j] - s0) / ds;
    let i = Math.floor(t);
    if (i < 0 || i >= ns - 1) {
      if (i === ns - 1 && t - i < 1e-9) { i = ns - 2; }
      else continue;                                // outside: leave zero
    }
    const w = t - i;
    for (let n = 0; n < nm; n++)
      out[n * nq + j] = F[n * ns + i] * (1 - w) + F[n * ns + i + 1] * w;
  }
  return out;
}

/** Set up everything that does not depend on the train. */
function prepare(inp) {
  inp = Object.assign({ Nmod: 3, fmax: Infinity, N_car: 20, W: 0 }, inp);
  const mdl = buildModel(inp);
  const { dL, EL, ndof, Ltot, Phi, nm } = mdl;
  const f = mdl.f;

  const fs = inp.fs || Math.ceil(20 * f[nm - 1]);
  const t_fr = inp.t_fr || Math.round(2 / f[nm - 1] * 100) / 100;
  const xs = inp.x || Array.from({ length: EL + 1 }, (_, i) => i * dL);

  const v_arr = inp.v.map(v => v / 3.6);
  return {
    inp, mdl, f, fs, t_fr, xs, nx: xs.length,
    Phi_x: hermiteEvalMany(Phi, nm, ndof, dL, EL, xs),
    zeta: inp.xi * 0.01,
    w: Array.from(f, v => 2 * Math.PI * v),
    v_arr, nv: v_arr.length,
    vmax: Math.max(...v_arr),
    ds: Math.min(...v_arr) / fs,
    nTrain: (inp.load === 'HSLM') ? 10 : 1
  };
}

/** One train: sweep every speed, return the two envelopes. */
function runTrain(ctx, nr) {
  const { inp, mdl, fs, t_fr, nx, Phi_x, zeta, w, v_arr, nv, vmax, ds } = ctx;
  const { dL, EL, ndof, Ltot, Phi, nm } = mdl;

  const F = distribute(trainLoad(inp, nr), inp.W, dL);
  const Ltrain = F[F.length - 1][0];

  // reference modal load in position space: built once, reused at every speed
  const smax = Ltrain + Ltot + vmax * t_fr;
  const ns = Math.floor(smax / ds) + 1;
  const Fmod = new Float64Array(nm * ns);
  for (const [d, P] of F) {
    const j0 = Math.max(0, Math.ceil(d / ds));
    const j1 = Math.min(ns - 1, Math.floor((d + Ltot) / ds));
    for (let j = j0; j <= j1; j++) {
      const ph = hermiteEval(Phi, nm, ndof, dL, EL, j * ds - d);
      for (let n = 0; n < nm; n++) Fmod[n * ns + j] -= P * ph[n];
    }
  }

  const dr = new Float64Array(nv), ar = new Float64Array(nv);
  for (let vi = 0; vi < nv; vi++) {
    const v = v_arr[vi];
    const dt = 1 / fs;
    // resample onto this speed's own grid, so dt = 1/fs everywhere
    const nt = Math.floor(smax * fs / v) + 1;
    const sq = new Float64Array(nt);
    for (let j = 0; j < nt; j++) sq[j] = j * v / fs;
    const p = resample(Fmod, nm, 0, ds, ns, sq, nt);

    const q = new Float64Array(nm * nt), qdd = new Float64Array(nm * nt);
    for (let n = 0; n < nm; n++) {
      const co = njCoeff(w[n], zeta, dt);
      const pn = p.subarray(n * nt, (n + 1) * nt);
      const qn = filter2(co.bq, co.a, pn);
      const vn = filter2(co.bv, co.a, pn);
      const c1 = 2 * zeta * w[n], c2 = w[n] * w[n];
      for (let j = 0; j < nt; j++) {
        q[n * nt + j] = qn[j];
        qdd[n * nt + j] = pn[j] - c1 * vn[j] - c2 * qn[j];   // straight from the ODE
      }
    }

    // superposition + running extrema, without ever forming the nx-by-nt array
    let dm = -Infinity, am = 0;
    for (let j = 0; j < nt; j++) {
      for (let i = 0; i < nx; i++) {
        let wd = 0, wa = 0;
        for (let n = 0; n < nm; n++) {
          wd += Phi_x[n * nx + i] * q[n * nt + j];
          wa += Phi_x[n * nx + i] * qdd[n * nt + j];
        }
        if (wd > dm) dm = wd;
        const aa = wa < 0 ? -wa : wa;
        if (aa > am) am = aa;
      }
    }
    dr[vi] = dm; ar[vi] = am;
  }
  return { dr, ar };
}

function finish(ctx, dmax, amax) {
  const m = ctx.mdl;
  return {
    f: Array.from(m.f), Phi: m.Phi, nm: m.nm, ndof: m.ndof, dL: m.dL, EL: m.EL,
    L: m.L, Ltot: m.Ltot, x: ctx.xs, phi_x: ctx.Phi_x, v: ctx.inp.v,
    dmax, amax, fs: ctx.fs, t_fr: ctx.t_fr, inp: ctx.inp
  };
}

/** Full analysis, synchronous. */
function run(inp, onProgress) {
  const ctx = prepare(inp);
  const dmax = [], amax = [];
  for (let nr = 0; nr < ctx.nTrain; nr++) {
    if (onProgress) onProgress(nr / ctx.nTrain,
      inp.load === 'HSLM' ? 'HSLM-A' + (nr + 1) : inp.load);
    const r = runTrain(ctx, nr);
    dmax.push(r.dr); amax.push(r.ar);
  }
  if (onProgress) onProgress(1, 'done');
  return finish(ctx, dmax, amax);
}

/** Full analysis, yielding to the event loop between trains so the UI paints. */
async function runAsync(inp, onProgress) {
  const ctx = prepare(inp);
  const dmax = [], amax = [];
  for (let nr = 0; nr < ctx.nTrain; nr++) {
    if (onProgress) onProgress(nr / ctx.nTrain,
      inp.load === 'HSLM' ? 'HSLM-A' + (nr + 1) : inp.load);
    await new Promise(r => setTimeout(r, 0));
    const r = runTrain(ctx, nr);
    dmax.push(r.dr); amax.push(r.ar);
  }
  if (onProgress) onProgress(1, 'done');
  return finish(ctx, dmax, amax);
}

/* ===================== EXPORTS ======================================== */

const API = {
  cholesky, forwardSolveMat, backSolveVec, jacobiEig,
  buildModel, shapeFun, hermiteEval, hermiteEvalMany,
  njCoeff, filter2, triangWin, trainLoad, distribute, resample,
  prepare, runTrain, run, runAsync,
  LINE_CAT
};

if (typeof module !== 'undefined' && module.exports) module.exports = API;
else root.BeamDyn = API;

})(typeof self !== 'undefined' ? self : this);

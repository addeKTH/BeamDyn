/* beamdyn2d.js --------------------------------------------------------------
   Dynamic response of a 2D beam to passing trains: JavaScript port of
   BeamDyn2D.m (2026-11-01) as used by the GUI (mode-acceleration correction
   and section forces always on).

   METHOD
     Bernoulli-Euler beam elements, consistent mass.  Moving constant forces:
     no train-bridge interaction, no track irregularity.  The modal load and
     the static response are built once per train as functions of TRAIN
     POSITION, sampled at ds = min(v)/fs, and resampled (linearly) onto each
     speed's own time grid, dt = 1/fs.  Each modal equation is integrated with
     the exact Nigam-Jennings recurrence for a piecewise-linear load, applied
     as a 2-pole IIR filter; the modal acceleration follows from the ODE.

   MODE ACCELERATION
     u = inv(K)*f(t) + sum_n phi_n*(q_n - p_n/w_n^2) for d, rotation, M, V and
     the reactions; the acceleration is the plain modal sum.  The static part
     inv(K)*f also gives the static peaks for phi' = y_dyn/y_stat - 1.

   SECTION FORCES
     M and V at both ends of every element.  In the static part the element
     end forces are exact, f = ke*u - f_eq: the consistent load of the axles
     inside an element is included, so M and V at loads and supports do not
     depend on dL.

   COORDINATES
     x along the beam to the right, w vertical and positive upwards, axle
     loads downwards, so a sagging moment is positive.  dmax, amax, rmax,
     Vmax and RMmax are magnitudes; Mmax, RFmax signed maxima; Mmin, RFmin
     signed minima.

   TRAIN LOADS AND LOAD DISTRIBUTION: see trainLoads and loadDistribution.

   Andreas Andersson, 2026-11-01 (assisted by Claude)
---------------------------------------------------------------------------*/
(function (root) {
'use strict';

const LA = root.LA || (typeof require === 'function' ? require('./linalg.js') : null);

const BD = {};
const BW = 3;                                       // half bandwidth
BD.FIELDS = ['dmax', 'amax', 'rmax', 'Mmax', 'Mmin', 'Vmax', 'RFmax', 'RFmin', 'RMmax'];

// ================================================================ the model
// Geometry, assembly and supports.  inp: {L: [m], m, EI, dL, clamped:
// [bool per support]} (clamped may be shorter; missing = pinned).
BD.buildModel = function (inp) {
    const dL = inp.dL;
    const L = inp.L.map(v => dL * Math.round(v / dL));        // spans on the grid
    if (L.some(v => !(v > 0))) throw new Error('a span is shorter than dL/2');
    const sumL = L.reduce((a, b) => a + b, 0);
    const EL = Math.round(sumL / dL), ndof = 2 * (EL + 1);
    const xs = [0];
    for (const v of L) xs.push(xs[xs.length - 1] + v);
    const clamped = xs.map((_, i) => !!(inp.clamped && inp.clamped[i]));
    const bc = xs.map(x => [x, 1]).concat(xs.filter((_, i) => clamped[i]).map(x => [x, 2]));

    const { K, M } = assemble(EL, dL, inp.EI, inp.m);
    const fixed = new Set(bc.map(([x, dir]) => Math.round(2 * x / dL) + dir - 1));
    const fdof = [];
    for (let i = 0; i < ndof; i++) if (!fixed.has(i)) fdof.push(i);
    return { L, dL, EL, ndof, sumL, xs, clamped, bc, K, M, fdof, EI: inp.EI, m: inp.m };
};

function assemble(EL, dL, EI, m) {
    const ndof = 2 * (EL + 1);
    const k1 = 12 * EI / dL ** 3, k2 = 6 * EI / dL ** 2, k3 = 4 * EI / dL;
    const ke = [[k1, k2, -k1, k2], [k2, k3, -k2, k3 / 2], [-k1, -k2, k1, -k2], [k2, k3 / 2, -k2, k3]];
    const c = m * dL / 420;
    const me = [[156, 22 * dL, 54, -13 * dL], [22 * dL, 4 * dL * dL, 13 * dL, -3 * dL * dL],
                [54, 13 * dL, 156, -22 * dL], [-13 * dL, -3 * dL * dL, -22 * dL, 4 * dL * dL]];
    const K = LA.bandAlloc(ndof, BW), M = LA.bandAlloc(ndof, BW);
    for (let e = 0; e < EL; e++)
        for (let a = 0; a < 4; a++)
            for (let b = a; b < 4; b++) {
                LA.bandAdd(K, BW, 2 * e + a, 2 * e + b, ke[a][b]);
                LA.bandAdd(M, BW, 2 * e + a, 2 * e + b, c * me[a][b]);
            }
    return { K, M };
}

// =================================================== eigenvalue analysis
// Lowest modes of K phi = lam M phi on the free dofs.  opt: {nreq: number of
// modes wanted, fmax: keep at least every mode up to fmax (Hz) when the
// problem is large}.  Small problems (<= 600 free dofs) are solved densely
// (whole spectrum, as eig in MATLAB); large ones by subspace iteration for
// the modes below fmax (Sturm count), at least 3.  Returns f (Hz, ascending)
// and Phi [ndof x n] column major, mass-normalised, the largest vertical
// component of each mode positive.
BD.eigen = function (mdl, opt) {
    const nf = mdl.fdof.length;
    if (nf < 1) throw new Error('every dof is restrained');
    const Kf = LA.bandReduce(mdl.K, BW, mdl.fdof), Mf = LA.bandReduce(mdl.M, BW, mdl.fdof);
    let nreq = Math.min(Math.max(1, opt.nreq), nf), lam, Vf;
    if (nf <= 600) {
        const r = LA.genEig(LA.bandToDense(Kf, nf, BW), LA.bandToDense(Mf, nf, BW), nf);
        lam = r.lam.subarray(0, nreq);
        Vf = new Float64Array(nf * nreq);                       // column major
        for (let c = 0; c < nreq; c++) for (let i = 0; i < nf; i++) Vf[c * nf + i] = r.Y[i * nf + c];
    } else {
        let nb = nreq;
        if (opt.fmaxFn) {                                       // f_max from the lowest modes
            const r3 = subspace(Kf, Mf, nf, Math.min(3, nf));
            const fmax = opt.fmaxFn(Array.from(r3.lam, l => Math.sqrt(Math.max(l, 0)) / (2 * Math.PI)));
            nb = Math.min(nreq, Math.max(3, sturmCount(Kf, Mf, nf, (2 * Math.PI * fmax) ** 2 * (1 + 1e-9))));
        }
        const r = subspace(Kf, Mf, nf, Math.min(nb, nf));
        lam = r.lam; Vf = r.X; nreq = lam.length;
    }
    const ndof = mdl.ndof;
    const f = Float64Array.from(lam, l => Math.sqrt(Math.max(l, 0)) / (2 * Math.PI));
    const Phi = new Float64Array(ndof * nreq);
    const tmp = new Float64Array(ndof);
    for (let c = 0; c < nreq; c++) {
        const col = Phi.subarray(c * ndof, (c + 1) * ndof);
        for (let i = 0; i < nf; i++) col[mdl.fdof[i]] = Vf[c * nf + i];
        // explicit mass normalisation, sign: largest vertical component positive
        let s = 0;
        for (let i = 0; i < ndof; i++) tmp[i] = LA.bandRowDot(mdl.M, ndof, BW, i, col);
        for (let i = 0; i < ndof; i++) s += col[i] * tmp[i];
        let big = 0, sg = 1;
        for (let i = 0; i < ndof; i += 2) if (Math.abs(col[i]) > big) { big = Math.abs(col[i]); sg = Math.sign(col[i]) || 1; }
        const sc = sg / Math.sqrt(s);
        for (let i = 0; i < ndof; i++) col[i] *= sc;
    }
    return { f, Phi };
};

// number of eigenvalues below sigma: negative pivots of LDL' of K - sigma*M
function sturmCount(Kf, Mf, n, sigma) {
    const w = BW + 1, A = new Float64Array(Kf.length);
    for (let i = 0; i < A.length; i++) A[i] = Kf[i] - sigma * Mf[i];
    let neg = 0;
    for (let i = 0; i < n; i++) {
        let d = A[i * w];
        if (d === 0) d = 1e-300;
        if (d < 0) neg++;
        for (let j = i + 1; j <= Math.min(i + BW, n - 1); j++) {
            const l = A[i * w + (j - i)] / d;
            for (let k = j; k <= Math.min(i + BW, n - 1); k++) A[j * w + (k - j)] -= l * A[i * w + (k - i)];
        }
    }
    return neg;
}

// subspace iteration for the p lowest modes (Bathe), dense projected problem
function subspace(Kf, Mf, n, p) {
    const q = Math.min(n, Math.max(2 * p, p + 8));
    const R = LA.bandChol(Kf, n, BW);
    const mul = (A, x, y) => { for (let i = 0; i < n; i++) y[i] = LA.bandRowDot(A, n, BW, i, x); return y; };
    let X = new Float64Array(n * q);
    const ratio = new Float64Array(n);
    for (let i = 0; i < n; i++) {
        const ki = Kf[i * (BW + 1)], mi = Mf[i * (BW + 1)];
        ratio[i] = ki ? mi / ki : 0; X[i] = mi;
    }
    const order = Array.from({ length: n }, (_, i) => i).sort((a, b) => ratio[b] - ratio[a]);
    for (let c = 1; c < q - 1; c++) X[c * n + order[c - 1]] = 1;
    for (let i = 0; i < n; i++) X[(q - 1) * n + i] = (i % 2 ? -1 : 1) * (1 + i / n);
    const Xb = new Float64Array(n * q), MX = new Float64Array(n * q), tmp = new Float64Array(n);
    const Kr = new Float64Array(q * q), Mr = new Float64Array(q * q);
    let lamOld = new Float64Array(q).fill(Infinity), lam = null;
    for (let it = 0; it < 500; it++) {
        for (let c = 0; c < q; c++) {
            mul(Mf, X.subarray(c * n, (c + 1) * n), tmp);
            LA.bandCholSolve(R, n, BW, tmp, Xb.subarray(c * n, (c + 1) * n));
        }
        for (let c = 0; c < q; c++) { mul(Kf, Xb.subarray(c * n, (c + 1) * n), tmp); MX.set(tmp, c * n); }
        for (let a = 0; a < q; a++) for (let b = a; b < q; b++) {
            let s = 0; for (let i = 0; i < n; i++) s += Xb[a * n + i] * MX[b * n + i];
            Kr[a * q + b] = s; Kr[b * q + a] = s;
        }
        for (let c = 0; c < q; c++) { mul(Mf, Xb.subarray(c * n, (c + 1) * n), tmp); MX.set(tmp, c * n); }
        for (let a = 0; a < q; a++) for (let b = a; b < q; b++) {
            let s = 0; for (let i = 0; i < n; i++) s += Xb[a * n + i] * MX[b * n + i];
            Mr[a * q + b] = s; Mr[b * q + a] = s;
        }
        const e = LA.genEig(Kr, Mr, q);
        lam = e.lam;
        const Xn = new Float64Array(n * q);
        for (let c = 0; c < q; c++)
            for (let k = 0; k < q; k++) {
                const qk = e.Y[k * q + c];
                if (qk) for (let i = 0; i < n; i++) Xn[c * n + i] += Xb[k * n + i] * qk;
            }
        X = Xn;
        let done = it > 1;
        for (let c = 0; c < p && done; c++)
            if (Math.abs(lam[c] - lamOld[c]) > 1e-13 * Math.abs(lam[c])) done = false;
        lamOld = Float64Array.from(lam);
        if (done) break;
    }
    return { lam: lam.subarray(0, p), X: X.subarray(0, n * p) };
}

// Mode shape at arbitrary x (element Hermite polynomials: exact between
// nodes).  Phi column major [ndof x n]; returns Float64Array(x.length).
BD.modeShape = function (mdl, Phi, n, x) {
    const { ndof, dL, EL } = mdl, c = n * ndof, out = new Float64Array(x.length);
    for (let p = 0; p < x.length; p++) {
        const e = Math.min(Math.max(Math.floor(x[p] / dL), 0), EL - 1);
        const s = x[p] / dL - e, s2 = s * s, s3 = s * s2, i1 = c + 2 * e;
        out[p] = (1 - 3 * s2 + 2 * s3) * Phi[i1] + dL * (s - 2 * s2 + s3) * Phi[i1 + 1] +
                 (3 * s2 - 2 * s3) * Phi[i1 + 2] + dL * (s3 - s2) * Phi[i1 + 3];
    }
    return out;
};

// ===================================================== recovery operators
// Per point the first dof and four Hermite coefficients for der = 0 (w),
// 1 (w'), 2 (w''), 3 (w''').  e/sl give the element (0 based) and the local
// coordinate explicitly; otherwise the element to the right (the last one
// for the end node).
function recov(x, dL, EL, der, scale, eIn, slIn) {
    const np = x.length, dof = new Int32Array(np), N = new Float64Array(4 * np);
    for (let p = 0; p < np; p++) {
        let e, sl;
        if (eIn) { e = eIn[p]; sl = slIn[p]; }
        else { e = Math.min(Math.max(Math.floor(x[p] / dL), 0), EL - 1); sl = x[p] / dL - e; }
        const s2 = sl * sl, s3 = sl * s2;
        let n;
        if (der === 0) n = [1 - 3 * s2 + 2 * s3, dL * (sl - 2 * s2 + s3), 3 * s2 - 2 * s3, dL * (s3 - s2)];
        else if (der === 1) n = [(-6 * sl + 6 * s2) / dL, 1 - 4 * sl + 3 * s2, (6 * sl - 6 * s2) / dL, 3 * s2 - 2 * sl];
        else if (der === 2) n = [(-6 + 12 * sl) / dL ** 2, (-4 + 6 * sl) / dL, (6 - 12 * sl) / dL ** 2, (-2 + 6 * sl) / dL];
        else n = [12 / dL ** 3, 6 / dL ** 2, -12 / dL ** 3, 6 / dL ** 2];
        dof[p] = 2 * e;
        for (let k = 0; k < 4; k++) N[4 * p + k] = scale * n[k];
    }
    return { np, dof, N };
}

function applyRecov(R, u, out, o) {
    for (let p = 0; p < R.np; p++) {
        const d = R.dof[p], q = 4 * p;
        out[o + p] = R.N[q] * u[d] + R.N[q + 1] * u[d + 1] + R.N[q + 2] * u[d + 2] + R.N[q + 3] * u[d + 3];
    }
}

// [np x Nmod] row major: the recovery applied to every mode
function recovModal(R, Phi, ndof, nm) {
    const out = new Float64Array(R.np * nm);
    for (let c = 0; c < nm; c++) {
        const o = c * ndof;
        for (let p = 0; p < R.np; p++) {
            const d = o + R.dof[p], q = 4 * p;
            out[p * nm + c] = R.N[q] * Phi[d] + R.N[q + 1] * Phi[d + 1] + R.N[q + 2] * Phi[d + 2] + R.N[q + 3] * Phi[d + 3];
        }
    }
    return out;
}

// ================================================ model for the dynamics
// Everything a train needs (also in a worker).  spec: {L (on the grid), dL,
// EI, m, clamped, f, Phi, xi (%), fs, t_fr, v (km/h), loadDist, b_Q}
BD.dynModel = function (spec) {
    const mdl = BD.buildModel({ L: spec.L, dL: spec.dL, EI: spec.EI, m: spec.m || 1, clamped: spec.clamped });
    const { dL, EL, ndof, sumL, K, bc, fdof } = mdl;
    const f = Float64Array.from(spec.f), Phi = Float64Array.from(spec.Phi), Nmod = f.length;
    const nfr = fdof.length;
    const Rchol = LA.bandChol(LA.bandReduce(K, BW, fdof), nfr, BW);

    const xd = Float64Array.from({ length: EL + 1 }, (_, i) => i * dL);
    const xr = Float64Array.from([...new Set(bc.map(b => b[0]))].sort((a, b) => a - b));
    const xF = new Float64Array(2 * EL), eF = new Int32Array(2 * EL), sF = new Float64Array(2 * EL);
    for (let e = 0; e < EL; e++) {
        xF[2 * e] = e * dL; eF[2 * e] = e; sF[2 * e] = 0;
        xF[2 * e + 1] = (e + 1) * dL; eF[2 * e + 1] = e; sF[2 * e + 1] = 1;
    }
    const Rd = recov(xd, dL, EL, 0, 1), Rr = recov(xr, dL, EL, 1, 1);
    const RM = recov(xF, dL, EL, 2, spec.EI, eF, sF), RV = recov(xF, dL, EL, 3, spec.EI, eF, sF);
    const pdofF = Int32Array.from([...new Set(bc.filter(b => b[1] === 1).map(b => Math.round(2 * b[0] / dL)))].sort((a, b) => a - b));
    const pdofM = Int32Array.from([...new Set(bc.filter(b => b[1] === 2).map(b => Math.round(2 * b[0] / dL) + 1))].sort((a, b) => a - b));
    const rowsTimesPhi = (rows) => {
        const out = new Float64Array(rows.length * Nmod);
        for (let c = 0; c < Nmod; c++) {
            const col = Phi.subarray(c * ndof, (c + 1) * ndof);
            for (let k = 0; k < rows.length; k++) out[k * Nmod + c] = LA.bandRowDot(K, ndof, BW, rows[k], col);
        }
        return out;
    };
    const v = Float64Array.from(spec.v, x => x / 3.6);
    let vmin = Infinity, vmax = -Infinity;
    for (const x of v) { vmin = Math.min(vmin, x); vmax = Math.max(vmax, x); }
    return {
        dL, EL, ndof, sumL, K, fdof, nfr, Rchol, f, Phi, Nmod,
        w2: Float64Array.from(f, x => (2 * Math.PI * x) ** 2),
        xi: spec.xi * 0.01, fs: spec.fs, t_fr: spec.t_fr, v, vmax, ds: vmin / spec.fs,
        xd, xr, xF, Rd, Rr, RM, RV, pdofF, pdofM,
        Psi_d: recovModal(Rd, Phi, ndof, Nmod), Psi_r: recovModal(Rr, Phi, ndof, Nmod),
        Psi_M: recovModal(RM, Phi, ndof, Nmod), Psi_V: recovModal(RV, Phi, ndof, Nmod),
        Psi_rF: rowsTimesPhi(pdofF), Psi_rM: rowsTimesPhi(pdofM),
        Dq: BD.loadDistribution(spec.loadDist || '1-point', spec.b_Q || 0, dL)
    };
};

// ============================================================ one train
// Speed-independent part of one train ({x, P}: axle positions from the
// first axle (m), axle loads (N)): the modal load and the static response
// as functions of train position, and the static peaks.
BD.trainPrep = function (train, P) {
    const { dL, EL, ndof, sumL, Nmod, Phi, ds } = P;
    // load distribution, then merge loads landing on the same position
    let xs = Array.from(train.x), Ps = Array.from(train.P);
    if (P.Dq.x.length > 1) {
        const xj = [], Pj = [];
        for (let k = 0; k < P.Dq.x.length; k++)
            for (let a = 0; a < xs.length; a++) { xj.push(xs[a] + P.Dq.x[k]); Pj.push(Ps[a] * P.Dq.w[k]); }
        let xmin = Infinity; for (const x of xj) xmin = Math.min(xmin, x);
        xs = xj.map(x => x - xmin); Ps = Pj;
    }
    const map = new Map();
    for (let i = 0; i < xs.length; i++) {
        const key = Math.round(xs[i] * 1e6) / 1e6;
        map.set(key, (map.get(key) || 0) + Ps[i]);
    }
    const ax = Float64Array.from([...map.keys()].sort((a, b) => a - b));
    const aP = Float64Array.from(ax, k => map.get(k));
    const na = ax.length;

    // s = position of the first axle; axle i is on the bridge while
    // 0 <= s - x_i <= sumL; the tail lets the last axle leave and the beam
    // vibrate freely for t_fr at the highest speed
    const smax = ax[na - 1] + sumL + P.vmax * P.t_fr;
    const Ns = Math.floor(smax / ds + 1e-9) + 1;

    const nx = P.Rd.np, nr = P.Rr.np, nF = P.RM.np, nsF = P.pdofF.length, nsM = P.pdofM.length;
    const F_mod = new Float64Array(Nmod * Ns);
    const d_st = new Float64Array(nx * Ns), r_st = new Float64Array(nr * Ns);
    const M_st = new Float64Array(nF * Ns), V_st = new Float64Array(nF * Ns);
    const rF_st = new Float64Array(nsF * Ns), rM_st = new Float64Array(nsM * Ns);
    const st = { dmax: 0, rmax: 0, Mmax: -Infinity, Mmin: Infinity, Vmax: 0, RFmax: 0, RFmin: 0, RMmax: 0 };
    if (nsF) { st.RFmax = -Infinity; st.RFmin = Infinity; }

    const rhs = new Float64Array(ndof), u = new Float64Array(ndof), uf = new Float64Array(P.nfr);
    const CM = new Float64Array(nF), CV = new Float64Array(nF);
    const touched = [], touchedF = [];
    for (let i = 0; i < Ns; i++) {
        const s = i * ds;
        for (let a = 0; a < na; a++) {
            const si = s - ax[a];
            if (si < 0 || si > sumL) continue;
            const e = Math.min(Math.max(Math.floor(si / dL), 0), EL - 1);
            const sl = si / dL - e, s2 = sl * sl, s3 = sl * s2;
            const N0 = 1 - 3 * s2 + 2 * s3, N1 = dL * (sl - 2 * s2 + s3), N2 = 3 * s2 - 2 * s3, N3 = dL * (s3 - s2);
            const d0 = 2 * e, Pa = aP[a];
            for (let n = 0; n < Nmod; n++) {                    // modal load, -P*phi_n(si)
                const c = n * ndof + d0;
                F_mod[n * Ns + i] -= Pa * (N0 * Phi[c] + N1 * Phi[c + 1] + N2 * Phi[c + 2] + N3 * Phi[c + 3]);
            }
            rhs[d0] -= Pa * N0; rhs[d0 + 1] -= Pa * N1; rhs[d0 + 2] -= Pa * N2; rhs[d0 + 3] -= Pa * N3;
            touched.push(d0);
            // exact end forces of the loaded element, f = ke*u - f_eq
            CM[2 * e] -= Pa * N1; CM[2 * e + 1] += Pa * N3;
            CV[2 * e] += Pa * N0; CV[2 * e + 1] -= Pa * N2;
            touchedF.push(2 * e);
        }
        // static solution and its recovery
        for (let j = 0; j < P.nfr; j++) uf[j] = rhs[P.fdof[j]];
        LA.bandCholSolve(P.Rchol, P.nfr, BW, uf, uf);
        u.fill(0);
        for (let j = 0; j < P.nfr; j++) u[P.fdof[j]] = uf[j];
        applyRecov(P.Rd, u, d_st, i * nx);
        applyRecov(P.Rr, u, r_st, i * nr);
        applyRecov(P.RM, u, M_st, i * nF);
        applyRecov(P.RV, u, V_st, i * nF);
        for (let p = 0; p < nF; p++) { M_st[i * nF + p] += CM[p]; V_st[i * nF + p] += CV[p]; }
        for (let k = 0; k < nsF; k++) rF_st[i * nsF + k] = LA.bandRowDot(P.K, ndof, BW, P.pdofF[k], u) - rhs[P.pdofF[k]];
        for (let k = 0; k < nsM; k++) rM_st[i * nsM + k] = LA.bandRowDot(P.K, ndof, BW, P.pdofM[k], u) - rhs[P.pdofM[k]];
        // static peaks
        for (let p = 0; p < nx; p++) st.dmax = Math.max(st.dmax, Math.abs(d_st[i * nx + p]));
        for (let p = 0; p < nr; p++) st.rmax = Math.max(st.rmax, Math.abs(r_st[i * nr + p]));
        for (let p = 0; p < nF; p++) {
            const m = M_st[i * nF + p];
            if (m > st.Mmax) st.Mmax = m;
            if (m < st.Mmin) st.Mmin = m;
            st.Vmax = Math.max(st.Vmax, Math.abs(V_st[i * nF + p]));
        }
        for (let k = 0; k < nsF; k++) { const r = rF_st[i * nsF + k]; if (r > st.RFmax) st.RFmax = r; if (r < st.RFmin) st.RFmin = r; }
        for (let k = 0; k < nsM; k++) st.RMmax = Math.max(st.RMmax, Math.abs(rM_st[i * nsM + k]));
        // reset the touched entries
        for (const d of touched) { rhs[d] = 0; rhs[d + 1] = 0; rhs[d + 2] = 0; rhs[d + 3] = 0; }
        for (const p of touchedF) { CM[p] = 0; CM[p + 1] = 0; CV[p] = 0; CV[p + 1] = 0; }
        touched.length = 0; touchedF.length = 0;
    }
    if (!nsM) st.RMmax = null;
    return { F_mod, Ns, smax, st, d_st, r_st, M_st, V_st, rF_st, rM_st };
};

// Nigam-Jennings recurrence for qdd + 2 z w qd + w^2 q = p(t), p linear over
// each step dt (unit modal mass), as digital filter coefficients
function njCoeff(w, z, dt) {
    const sq = Math.sqrt(1 - z * z), wd = w * sq, e = Math.exp(-z * w * dt);
    const s = Math.sin(wd * dt), c = Math.cos(wd * dt), zz = z / sq, k = w * w;
    const A = e * (zz * s + c), B = e * (s / wd);
    const Cc = (2 * z / (w * dt) + e * (((1 - 2 * z * z) / (wd * dt) - zz) * s - (1 + 2 * z / (w * dt)) * c)) / k;
    const Dd = (1 - 2 * z / (w * dt) + e * (((2 * z * z - 1) / (wd * dt)) * s + (2 * z / (w * dt)) * c)) / k;
    const Ap = -e * (w / sq) * s, Bp = e * (c - zz * s);
    const Cp = (-1 / dt + e * ((w / sq + zz / dt) * s + c / dt)) / k;
    const Dp = (1 - e * (zz * s + c)) / (k * dt);
    return {
        a: [1, -(A + Bp), A * Bp - Ap * B],
        bq: [Dd, Cc - Bp * Dd + B * Dp, B * Cp - Bp * Cc],
        bv: [Dp, Cp + Ap * Dd - A * Dp, Ap * Cc - A * Cp]
    };
}

// max and min over all points and times of lin(stat) + Psi*Q, without
// forming the response; stat null = modal sum only
function peakOf(Psi, np, nm, Q, Nt, stat, i0, wg) {
    let mx = -Infinity, mn = Infinity;
    if (!np) return [0, 0];
    const qn = new Float64Array(nm);
    for (let t = 0; t < Nt; t++) {
        for (let n = 0; n < nm; n++) qn[n] = Q[n * Nt + t];
        const b0 = stat ? i0[t] * np : 0, b1 = b0 + np, w1 = stat ? wg[t] : 0, w0 = 1 - w1;
        for (let p = 0, r = 0; p < np; p++, r += nm) {
            let s = stat ? stat[b0 + p] * w0 + stat[b1 + p] * w1 : 0;
            for (let n = 0; n < nm; n++) s += Psi[r + n] * qn[n];
            if (s > mx) mx = s;
            if (s < mn) mn = s;
        }
    }
    return [mx, mn];
}
const absMax = r => Math.max(Math.abs(r[0]), Math.abs(r[1]));

// Peaks at the speeds P.v[iv] of the train prepared in Q (the other speeds
// are left 0).  check(): optional, true stops with an Error 'stopped'.
BD.speedPeaks = function (Q, P, iv, check) {
    const { Nmod, f, xi, fs, ds } = P, nv = P.v.length, Ns = Q.Ns;
    const T = {};
    for (const k of BD.FIELDS) T[k] = new Float64Array(nv);
    const nx = P.Rd.np, nr = P.Rr.np, nF = P.RM.np, nsF = P.pdofF.length, nsM = P.pdofM.length;
    const dt = 1 / fs;
    const co = Array.from(f, x => njCoeff(2 * Math.PI * x, xi, dt));
    for (const vi of iv) {
        const v = P.v[vi];
        const Nt = Math.floor(Q.smax * fs / v) + 1;
        const i0 = new Int32Array(Nt), wg = new Float64Array(Nt);
        for (let t = 0; t < Nt; t++) {
            const tt = (t * (v / fs)) / ds;
            const j = Math.min(Math.floor(tt), Ns - 2);
            i0[t] = j; wg[t] = tt - j;
        }
        const p = new Float64Array(Nt), q = new Float64Array(Nt), vn = new Float64Array(Nt);
        const qc = new Float64Array(Nmod * Nt), qdd = new Float64Array(Nmod * Nt);
        for (let n = 0; n < Nmod; n++) {
            const off = n * Ns, w = 2 * Math.PI * f[n];
            for (let t = 0; t < Nt; t++) p[t] = Q.F_mod[off + i0[t]] * (1 - wg[t]) + Q.F_mod[off + i0[t] + 1] * wg[t];
            LA.filter3(co[n].bq, co[n].a, p, Nt, q);
            LA.filter3(co[n].bv, co[n].a, p, Nt, vn);
            const o = n * Nt;
            for (let t = 0; t < Nt; t++) {
                qdd[o + t] = p[t] - 2 * xi * w * vn[t] - w * w * q[t];
                qc[o + t] = q[t] - p[t] / P.w2[n];
            }
        }
        T.dmax[vi] = absMax(peakOf(P.Psi_d, nx, Nmod, qc, Nt, Q.d_st, i0, wg));
        T.amax[vi] = absMax(peakOf(P.Psi_d, nx, Nmod, qdd, Nt, null));
        T.rmax[vi] = absMax(peakOf(P.Psi_r, nr, Nmod, qc, Nt, Q.r_st, i0, wg));
        const m = peakOf(P.Psi_M, nF, Nmod, qc, Nt, Q.M_st, i0, wg);
        T.Mmax[vi] = m[0]; T.Mmin[vi] = m[1];
        T.Vmax[vi] = absMax(peakOf(P.Psi_V, nF, Nmod, qc, Nt, Q.V_st, i0, wg));
        if (nsF) { const r = peakOf(P.Psi_rF, nsF, Nmod, qc, Nt, Q.rF_st, i0, wg); T.RFmax[vi] = r[0]; T.RFmin[vi] = r[1]; }
        if (nsM) T.RMmax[vi] = absMax(peakOf(P.Psi_rM, nsM, Nmod, qc, Nt, Q.rM_st, i0, wg));
        if (check && check()) throw new Error('stopped');
    }
    if (!nsM) T.RMmax = null;
    T.stat = Q.st;
    return T;
};

// ===================================================== dynamic factor
// phi' = y_dyn/y_stat - 1 per train and speed; dyn fields [ntr][nv], stat
// fields [ntr]; non-finite ratios (zero static peak) give NaN
BD.dynamicFactor = function (stat, dyn) {
    const phi = {};
    for (const k of ['dmax', 'rmax', 'Mmax', 'Mmin', 'Vmax', 'RFmax', 'RFmin', 'RMmax']) {
        if (!dyn[k] || !stat[k] || !dyn[k].length || !dyn[k][0] || stat[k][0] === null) { phi[k] = null; continue; }
        phi[k] = dyn[k].map((row, t) => Float64Array.from(row, y => {
            const r = y / stat[k][t];
            return isFinite(r) ? r - 1 : NaN;
        }));
    }
    return phi;
};

// ======================================================== train loads
// Array of trains {name, x: axle positions from the leading axle (m), P:
// axle loads (N)}.  'HSLM-A' (cars = train numbers, default 1-10): EN
// 1991-2:2024 Table 8.4, Fig. 8.12 - power car with 4 axles at 3, 11, 3 m;
// 3.525 m to the outer bogie of the end coach; N coaches of length D on
// articulation bogies of axle spacing d; 2N+14 axles of P.  A line category
// A1-E5 (cars = wagons, default 10): EN 15528:2021 Table A.1, identical
// loaded reference wagons, 1 t = 9.81 kN.
BD.LOADS = ['HSLM-A', 'A1', 'B1', 'B2', 'C2', 'C3', 'C4', 'D2', 'D3', 'D4', 'D5', 'E4', 'E5'];
BD.trainLoads = function (name, cars) {
    const g = 9.81;
    if (name === 'HSLM-A') {
        const A = [[18, 18, 2.0, 170], [17, 19, 3.5, 200], [16, 20, 2.0, 180], [15, 21, 3.0, 190],
                   [14, 22, 2.0, 170], [13, 23, 2.0, 180], [13, 24, 2.0, 190], [12, 25, 2.5, 190],
                   [11, 26, 2.0, 210], [11, 27, 2.0, 210]];
        const list = (cars && cars.length) ? cars : [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        return list.map(k => {
            const [N, D, d, Pk] = A[k - 1];
            const pc = [0, 3, 14, 17], c0 = 17 + 3.525 / 2;
            const Lt = c0 + (N + 2) * D + 3.525 / 2 + 17;
            const x = pc.concat([17 + 3.525, 17 + 3.525 + d]);
            for (let b = 1; b <= N + 1; b++) { const sb = c0 + b * D; x.push(sb - d / 2, sb + d / 2); }
            x.push(Lt - 17 - 3.525 - d, Lt - 17 - 3.525);
            for (const p of pc.slice().reverse()) x.push(Lt - p);
            return { name: 'A' + k, x: Float64Array.from(x), P: new Float64Array(x.length).fill(Pk * 1e3) };
        });
    }
    //        P (t)  front  axle spacings        rear (m)
    const W = { A1: [16.0, 1.5, [1.8, 6.20, 1.8], 1.5], B1: [18.0, 1.5, [1.8, 7.80, 1.8], 1.5],
                B2: [18.0, 1.5, [1.8, 4.65, 1.8], 1.5], C2: [20.0, 1.5, [1.8, 5.90, 1.8], 1.5],
                C3: [20.0, 1.5, [1.8, 4.50, 1.8], 1.5], C4: [20.0, 1.5, [1.8, 3.40, 1.8], 1.5],
                D2: [22.5, 1.5, [1.8, 7.45, 1.8], 1.5], D3: [22.5, 1.5, [1.8, 5.90, 1.8], 1.5],
                D4: [22.5, 1.5, [1.8, 4.65, 1.8], 1.5], D5: [22.5, 1.5, [1.8, 3.62, 1.8], 1.5],
                E4: [25.0, 1.5, [1.8, 5.90, 1.8], 1.5], E5: [25.0, 1.5, [1.8, 4.75, 1.8], 1.5] };
    const w = W[name];
    if (!w) throw new Error('unknown load model ' + name);
    const xw = [0]; for (const s of w[2]) xw.push(xw[xw.length - 1] + s);
    const Lw = w[1] + xw[xw.length - 1] + w[3];
    const nW = Math.max(Math.round(cars || 10), 1), x = [];
    for (let k = 0; k < nW; k++) for (const a of xw) x.push(a + k * Lw);
    return [{ name, x: Float64Array.from(x), P: new Float64Array(x.length).fill(w[0] * 1e3 * g) }];
};

// =================================================== load distribution
// Offsets x (m) from the axle and weights w (sum 1) of one axle load:
//   '1-point'     the axle load
//   '3-point'     25/50/25 % at -a, 0, a, sleeper spacing a = 0.6 m
//                 (EN 1991-2:2024 8.3.6.2)
//   '5-point'     sleeper loads at 0, +-a, +-2a from the rail deflection of a
//                 Winkler beam, w ~ exp(-beta|x|)(cos + sin)(beta|x|), beta =
//                 (K/(8 EI_r))^(1/3), track stiffness K = 100 MN/m per rail,
//                 60E1 rail EI_r = 210E9*30.38E-6 Nm^2; weights scaled to 1
//   'Triangular'  triangle of base b_Q, step <= dL and b_Q/10
BD.loadDistribution = function (type, b_Q, dL) {
    const a = 0.6;
    switch (String(type).toLowerCase()) {
        case '1-point': return { x: [0], w: [1] };
        case '3-point': return { x: [-a, 0, a], w: [0.25, 0.5, 0.25] };
        case '5-point': {
            const K = 100e6, EIr = 210e9 * 30.38e-6, be = Math.cbrt(K / (8 * EIr));
            const x = [-2 * a, -a, 0, a, 2 * a];
            const w = x.map(v => Math.exp(-be * Math.abs(v)) * (Math.cos(be * Math.abs(v)) + Math.sin(be * Math.abs(v))));
            const s = w.reduce((p, c) => p + c, 0);
            return { x, w: w.map(v => v / s) };
        }
        case 'triangular': {
            if (!(b_Q > 0)) return { x: [0], w: [1] };
            const n = 2 * Math.ceil(b_Q / (2 * Math.min(dL, b_Q / 10))) + 1;
            const x = [], w = [];
            for (let i = 0; i < n; i++) { const v = -b_Q / 2 + i * b_Q / (n - 1); x.push(v); w.push(1 - Math.abs(v) / (b_Q / 2)); }
            const s = w.reduce((p, c) => p + c, 0);
            return { x, w: w.map(v => v / s) };
        }
        default: throw new Error('unknown load distribution ' + type);
    }
};

// ========================================================== Eurocode
BD.EN = {
    // determinant length L_Phi (m), EN 1991-2:2024 Table 8.2, 5.2: k*L_m for n
    // continuous spans (k = 1.2/1.3/1.4/1.5 for n = 2/3/4/>=5), not less than
    // the longest span; a single span gives L
    LPhi(L) {
        if (L.length === 1) return L[0];
        const k = [1.2, 1.3, 1.4, 1.5][Math.min(L.length, 5) - 2];
        return Math.max(k * L.reduce((a, b) => a + b, 0) / L.length, Math.max(...L));
    },
    // lower bound damping (%), EN 1991-2:2024 Table 8.6
    xi(type, L) {
        const t = { 'Steel': [0.5, 0.125], 'Concrete': [1.5, 0.07], 'Pre-stressed': [1.0, 0.07] }[type];
        return t[0] + t[1] * Math.max(20 - L, 0);
    },
    // phi'' for track defects, Annex C (C.6)-(C.7); v in km/h, alpha_v =
    // v/22 up to 22 m/s, else 1
    phiBis(v, n0, Lphi) {
        const av = Math.min(v / 3.6 / 22, 1);
        const p = av / 100 * (56 * Math.exp(-((Lphi / 10) ** 2)) + 50 * (Lphi * n0 / 80 - 1) * Math.exp(-((Lphi / 20) ** 2)));
        return Math.max(p, 0);
    }
};

if (typeof module !== 'undefined' && module.exports) module.exports = BD;
root.BD = BD;
})(typeof self !== 'undefined' ? self : this);

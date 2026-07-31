/* beamdyn2d.js --------------------------------------------------------------
   Dynamic response of a 2D beam to passing trains.  Port of BeamDyn2D.m.

   METHOD
     Bernoulli-Euler beam elements, consistent mass.  Moving constant forces:
     no train-bridge interaction, no track irregularity.  The modal load is
     built as a function of TRAIN POSITION, sampled at ds = min(v)/fs, so speed
     enters only through the time step dt = ds/v and the load vector is
     speed-independent.  Each modal equation is integrated with the exact
     Nigam-Jennings recurrence for a piecewise-linear load, applied as a 2-pole
     IIR filter.  Modal acceleration follows algebraically from the ODE.

   MODE ACCELERATION (statCorr)
     Modal truncation converges as sum(n^-4) for displacement but only
     sum(n^-2) for the moment and sum(n^-1) for the shear.  The static
     correction adds back the exact static response of the truncated modes,
         u = inv(K)*f(t) + sum_n phi_n*(q_n - p_n/w_n^2)
     which is a function of train position like the modal load, so it is built
     once per train and resampled per speed.

   PEAKS
     Only extremes leave the solver, so the recovery is fused with the
     reduction: the [output points x time steps] response is never formed.
     dmax, amax, rmax, Vmax and RMmax are magnitudes; Mmax and RFmax are signed
     maxima, Mmin and RFmin signed minima (a negative RFmin is net uplift).
---------------------------------------------------------------------------*/
(function (root) {
'use strict';

const req = (typeof require === 'function' && typeof module !== 'undefined');
const LA = req ? require('./linalg.js') : root.LA;
const BM = req ? require('./beammodel.js') : root.BM;

const BD = {};

// ---------------------------------------------------------------- the model
// Everything that depends on the bridge but not on the train or the speed:
// assembly, boundary conditions, the eigenvalue problem and the recovery
// operators.  This is what the GUI recomputes on every input change.
BD.buildModel = function (inp) {
    const g = BM.geometry(inp);
    const { dL, EL, ndof, sumL } = g;
    const bw = BM.BW;
    const { K, M } = BM.assemble(g, inp.EI, inp.m);
    const bcd = BM.bcDofs(g, inp.bc);
    if (bcd.nf < 1) throw new Error('every dof is restrained');

    const Kf = BM.reduce(K, ndof, bw, bcd.f2g);
    const Mf = BM.reduce(M, ndof, bw, bcd.f2g);

    const nreq = Math.min(Math.max(1, Math.round(inp.Nmod)), bcd.nf);
    const e = LA.eigBand(Kf, Mf, bcd.nf, bw, nreq);

    let f = Array.from(e.lam, l => Math.sqrt(Math.max(l, 0)) / (2 * Math.PI));
    let keep = f.map(v => v < (inp.fmax === undefined ? Infinity : inp.fmax));
    if (!keep.some(Boolean)) keep[0] = true;
    const idx = keep.map((k, i) => k ? i : -1).filter(i => i >= 0);
    const Nmod = idx.length;
    f = Float64Array.from(idx, i => f[i]);

    const Phi = new Float64Array(ndof * Nmod);          // column major
    for (let c = 0; c < Nmod; c++) {
        const src = idx[c] * bcd.nf;
        for (let j = 0; j < bcd.nf; j++) Phi[c * ndof + bcd.f2g[j]] = e.Phi[src + j];
    }

    // output points
    const x = inp.x ? Float64Array.from(inp.x)
                    : Float64Array.from({ length: EL + 1 }, (_, i) => i * dL);
    // M and V at both ends of every element: M is linear and V constant inside
    // an element, so the ends bracket every extremum and keep the jump at a node
    const xF = new Float64Array(2 * EL), eF = new Int32Array(2 * EL), sF = new Float64Array(2 * EL);
    for (let el = 0; el < EL; el++) {
        xF[2 * el] = el * dL; eF[2 * el] = el; sF[2 * el] = 0;
        xF[2 * el + 1] = (el + 1) * dL; eF[2 * el + 1] = el; sF[2 * el + 1] = 1;
    }
    const Rd = BM.recov(x, dL, EL, 0);
    const Rr = BM.recov(x, dL, EL, 1);
    const RM = BM.recov(xF, dL, EL, 2, eF, sF); RM.scale = inp.EI;
    // The shear is constant within an element - the value at both ends is the
    // same number - so it is evaluated once per element instead of twice.
    const xV = new Float64Array(EL), eV = new Int32Array(EL), sV = new Float64Array(EL);
    for (let el = 0; el < EL; el++) { xV[el] = el * dL; eV[el] = el; }
    const RV = BM.recov(xV, dL, EL, 3, eV, sV); RV.scale = inp.EI;

    // reaction dofs, split into vertical (force) and rotational (moment)
    const setF = new Set(), setM = new Set();
    for (const [pos, dir] of inp.bc) {
        const d = Math.round(2 * pos / dL) + dir - 1;
        (dir === 1 ? setF : setM).add(d);
    }
    const pdofF = Int32Array.from(Array.from(setF).sort((a, b) => a - b));
    const pdofM = Int32Array.from(Array.from(setM).sort((a, b) => a - b));

    const fs = inp.fs !== undefined && inp.fs ? inp.fs : Math.ceil(20 * f[Nmod - 1]);
    const t_fr = inp.t_fr !== undefined && inp.t_fr !== null
        ? inp.t_fr : Math.round(2 / f[Nmod - 1] * 100) / 100;

    return {
        g, bw, K, M, Kf, Mf, bcd, ndof, dL, EL, sumL,
        f, Phi, Nmod, x, xF, Rd, Rr, RM, RV, pdofF, pdofM, fs, t_fr,
        EI: inp.EI, Rchol: null
    };
};

// Mode shape n (0 based) evaluated at arbitrary x, for plotting.
BD.modeShape = function (model, n, x) {
    const col = model.Phi.subarray(n * model.ndof, (n + 1) * model.ndof);
    return BM.hermiteEval(col, model.ndof, 1, model.dL, model.EL, x);
};

// -------------------------------------- Nigam-Jennings recurrence, as an IIR
// qdd + 2*z*w*qd + w^2*q = p(t), p linear over each step dt, returned as
// digital filter coefficients sharing one denominator.
function njCoeff(w, z, dt) {
    const sq = Math.sqrt(1 - z * z);
    const wd = w * sq;
    const e = Math.exp(-z * w * dt);
    const s = Math.sin(wd * dt), c = Math.cos(wd * dt);
    const zz = z / sq;
    const k = w * w;

    const A = e * (zz * s + c);
    const B = e * (s / wd);
    const Cc = (2 * z / (w * dt) + e * (((1 - 2 * z * z) / (wd * dt) - zz) * s - (1 + 2 * z / (w * dt)) * c)) / k;
    const Dd = (1 - 2 * z / (w * dt) + e * (((2 * z * z - 1) / (wd * dt)) * s + (2 * z / (w * dt)) * c)) / k;

    const Ap = -e * (w / sq) * s;
    const Bp = e * (c - zz * s);
    const Cp = (-1 / dt + e * ((w / sq + zz / dt) * s + c / dt)) / k;
    const Dp = (1 - e * (zz * s + c)) / (k * dt);

    return {
        a: [1, -(A + Bp), (A * Bp - Ap * B)],
        bq: [Dd, Cc - Bp * Dd + B * Dp, B * Cp - Bp * Cc],
        bv: [Dp, Cp + Ap * Dd - A * Dp, Ap * Cc - A * Cp]
    };
}

// ----------------------------------------------------- fused recovery + peak
// val(p,t) = sum_n Psi[p][n]*Q[n][t] + sign*interp(stat[:, i0(t)], i0+1)
// reduced to (max, min) over every output point and the whole record, without
// ever forming the response array.
function peakOf(Psi, np, nm, Q, Nt, stat, i0, wg, sign) {
    let mx = -Infinity, mn = Infinity;
    if (np === 0) return [0, 0];

    if (nm === 3 && !stat) {                         // the common case, unrolled
        for (let t = 0; t < Nt; t++) {
            const q0 = Q[t], q1 = Q[Nt + t], q2 = Q[2 * Nt + t];
            for (let p = 0, r = 0; p < np; p++, r += 3) {
                const s = Psi[r] * q0 + Psi[r + 1] * q1 + Psi[r + 2] * q2;
                if (s > mx) mx = s;
                if (s < mn) mn = s;
            }
        }
    } else if (nm === 3) {
        for (let t = 0; t < Nt; t++) {
            const q0 = Q[t], q1 = Q[Nt + t], q2 = Q[2 * Nt + t];
            const b0 = i0[t] * np, b1 = b0 + np, w1 = wg[t], w0 = 1 - w1;
            for (let p = 0, r = 0; p < np; p++, r += 3) {
                const s = Psi[r] * q0 + Psi[r + 1] * q1 + Psi[r + 2] * q2
                        + sign * (stat[b0 + p] * w0 + stat[b1 + p] * w1);
                if (s > mx) mx = s;
                if (s < mn) mn = s;
            }
        }
    } else if (!stat) {
        const qn = new Float64Array(nm);
        for (let t = 0; t < Nt; t++) {
            for (let n = 0; n < nm; n++) qn[n] = Q[n * Nt + t];
            for (let p = 0, r = 0; p < np; p++, r += nm) {
                let s = 0;
                for (let n = 0; n < nm; n++) s += Psi[r + n] * qn[n];
                if (s > mx) mx = s;
                if (s < mn) mn = s;
            }
        }
    } else {
        const qn = new Float64Array(nm);
        for (let t = 0; t < Nt; t++) {
            for (let n = 0; n < nm; n++) qn[n] = Q[n * Nt + t];
            const b0 = i0[t] * np, b1 = b0 + np, w1 = wg[t], w0 = 1 - w1;
            for (let p = 0, r = 0; p < np; p++, r += nm) {
                let s = sign * (stat[b0 + p] * w0 + stat[b1 + p] * w1);
                for (let n = 0; n < nm; n++) s += Psi[r + n] * qn[n];
                if (s > mx) mx = s;
                if (s < mn) mn = s;
            }
        }
    }
    if (mx === -Infinity) { mx = 0; mn = 0; }
    if (mn === Infinity) mn = mx;
    return [mx, mn];
}

const absMax = (r) => Math.max(Math.abs(r[0]), Math.abs(r[1]));

// ------------------------------------------------------------------- solver
// One train over a list of speeds.  vAll (km/h) is the FULL speed list of the
// run: ds and the record length are taken from it, so splitting the speeds
// over several workers gives bit-identical results to one call.
BD.runTrain = function (model, inp, train, vList, vAll, onProgress) {
    const { dL, EL, ndof, sumL, bw, Nmod, Phi, f } = model;
    const secForce = inp.secForces !== false;
    const statCorr = inp.statCorr !== false;
    const xi = inp.xi * 0.01;
    const fs = model.fs, t_fr = model.t_fr;

    const vArr = Float64Array.from(vList, v => v / 3.6);
    const vAllA = Float64Array.from(vAll || vList, v => v / 3.6);
    let vmin = Infinity, vmax = -Infinity;
    for (const v of vAllA) { if (v < vmin) vmin = v; if (v > vmax) vmax = v; }
    const ds = vmin / fs;

    const F = BM.prepareTrain(train, inp.W || 0, dL);
    const na = F.x.length;
    const smax = F.x[na - 1] + sumL + vmax * t_fr;
    const Ns = Math.floor(smax / ds) + 1;

    const nx = model.Rd.np, nF = model.RM.np, nV = model.RV.np;
    const nsF = model.pdofF.length, nsM = model.pdofM.length;
    const hasF = nsF > 0, hasM = nsM > 0;

    // modal shapes of each response quantity, [np x Nmod] row major
    const Psi_d = BM.recovModal(model.Rd, Phi, ndof, Nmod);
    const Psi_r = BM.recovModal(model.Rr, Phi, ndof, Nmod);
    const Psi_M = secForce ? BM.recovModal(model.RM, Phi, ndof, Nmod) : null;
    const Psi_V = secForce ? BM.recovModal(model.RV, Phi, ndof, Nmod) : null;
    const Psi_rF = new Float64Array(nsF * Nmod), Psi_rM = new Float64Array(nsM * Nmod);
    for (let c = 0; c < Nmod; c++) {
        const col = Phi.subarray(c * ndof, (c + 1) * ndof);
        for (let k = 0; k < nsF; k++) Psi_rF[k * Nmod + c] = LA.bandRowDot(model.K, ndof, bw, model.pdofF[k], col);
        for (let k = 0; k < nsM; k++) Psi_rM[k * Nmod + c] = LA.bandRowDot(model.K, ndof, bw, model.pdofM[k], col);
    }

    const w2 = new Float64Array(Nmod);
    for (let n = 0; n < Nmod; n++) { const w = 2 * Math.PI * f[n]; w2[n] = w * w; }

    // ---- load as a function of train position -----------------------------
    const F_mod = new Float64Array(Nmod * Ns);
    let d_st = null, r_st = null, M_st = null, V_st = null, rF_st = null, rM_st = null;
    let apF = null, apM = null;
    const needStatic = secForce && statCorr;
    const needApplied = secForce && !statCorr;

    if (needStatic) {
        d_st = new Float64Array(nx * Ns); r_st = new Float64Array(nx * Ns);
        M_st = new Float64Array(nF * Ns); V_st = new Float64Array(nV * Ns);
        rF_st = new Float64Array(nsF * Ns); rM_st = new Float64Array(nsM * Ns);
    }
    if (needApplied || needStatic) {
        apF = new Float64Array(nsF * Ns); apM = new Float64Array(nsM * Ns);
    }

    const rhs = new Float64Array(ndof);
    const touched = new Int32Array(4 * na);
    const uf = new Float64Array(model.bcd.nf), bf = new Float64Array(model.bcd.nf);
    const u = new Float64Array(ndof);
    const colBuf = new Float64Array(Math.max(nx, nF, nV));
    let Rchol = null;
    if (needStatic) {
        Rchol = model.Rchol || (model.Rchol = LA.bandChol(model.Kf, model.bcd.nf, bw));
    }

    for (let i = 0; i < Ns; i++) {
        const s = i * ds;
        let nt = 0;
        for (let a = 0; a < na; a++) {
            const si = s - F.x[a];
            if (si < 0 || si > sumL) continue;
            const el = Math.min(Math.max(Math.floor(si / dL), 0), EL - 1);
            const sl = si / dL - el, s2 = sl * sl, s3 = sl * s2;
            const N0 = 1 - 3 * s2 + 2 * s3, N1 = dL * (sl - 2 * s2 + s3);
            const N2 = 3 * s2 - 2 * s3, N3 = dL * (s3 - s2);
            const d0 = 2 * el, P = F.P[a];
            // modal load: -P * phi_n(si)
            for (let n = 0; n < Nmod; n++) {
                const c = n * ndof;
                F_mod[n * Ns + i] -= P * (N0 * Phi[c + d0] + N1 * Phi[c + d0 + 1] +
                                          N2 * Phi[c + d0 + 2] + N3 * Phi[c + d0 + 3]);
            }
            if (secForce) {                          // consistent nodal load
                rhs[d0] -= P * N0; rhs[d0 + 1] -= P * N1;
                rhs[d0 + 2] -= P * N2; rhs[d0 + 3] -= P * N3;
                touched[nt++] = d0; touched[nt++] = d0 + 1;
                touched[nt++] = d0 + 2; touched[nt++] = d0 + 3;
            }
        }
        if (!secForce) continue;

        for (let k = 0; k < nsF; k++) apF[i * nsF + k] = rhs[model.pdofF[k]];
        for (let k = 0; k < nsM; k++) apM[i * nsM + k] = rhs[model.pdofM[k]];

        if (needStatic) {
            for (let j = 0; j < model.bcd.nf; j++) bf[j] = rhs[model.bcd.f2g[j]];
            LA.bandCholSolve(Rchol, model.bcd.nf, bw, bf, uf);
            u.fill(0);
            for (let j = 0; j < model.bcd.nf; j++) u[model.bcd.f2g[j]] = uf[j];
            BM.applyRecov(model.Rd, u, colBuf); d_st.set(colBuf.subarray(0, nx), i * nx);
            BM.applyRecov(model.Rr, u, colBuf); r_st.set(colBuf.subarray(0, nx), i * nx);
            BM.applyRecov(model.RM, u, colBuf); M_st.set(colBuf.subarray(0, nF), i * nF);
            BM.applyRecov(model.RV, u, colBuf); V_st.set(colBuf.subarray(0, nV), i * nV);
            for (let k = 0; k < nsF; k++)
                rF_st[i * nsF + k] = LA.bandRowDot(model.K, ndof, bw, model.pdofF[k], u) - apF[i * nsF + k];
            for (let k = 0; k < nsM; k++)
                rM_st[i * nsM + k] = LA.bandRowDot(model.K, ndof, bw, model.pdofM[k], u) - apM[i * nsM + k];
        }
        for (let j = 0; j < nt; j++) rhs[touched[j]] = 0;
    }

    // ---- speed loop -------------------------------------------------------
    const nv = vList.length;
    const res = {
        dmax: new Float64Array(nv), amax: new Float64Array(nv), rmax: new Float64Array(nv),
        Mmax: new Float64Array(nv), Mmin: new Float64Array(nv), Vmax: new Float64Array(nv),
        RFmax: new Float64Array(nv), RFmin: new Float64Array(nv), RMmax: new Float64Array(nv)
    };

    // The slowest speed of this job gives the longest record; allocate for it
    // once and reuse the buffers, so the speed loop allocates nothing.
    let vLo = Infinity; for (const v of vArr) if (v < vLo) vLo = v;
    const NtMax = Math.floor(smax * fs / vLo) + 1;
    const i0 = new Int32Array(NtMax), wg = new Float64Array(NtMax);
    const p_v = new Float64Array(Nmod * NtMax);
    const q = new Float64Array(Nmod * NtMax), qdd = new Float64Array(Nmod * NtMax);
    const qc = (statCorr && secForce) ? new Float64Array(Nmod * NtMax) : null;
    const vn = new Float64Array(NtMax);

    for (let vi = 0; vi < nv; vi++) {
        const v = vArr[vi];
        const Nt = Math.floor(smax * fs / v) + 1;
        const dtds = (v / fs) / ds;
        for (let t = 0; t < Nt; t++) {
            const tt = t * dtds;
            const j = Math.min(Math.floor(tt), Ns - 2);
            i0[t] = j; wg[t] = tt - j;
        }
        for (let n = 0; n < Nmod; n++) {
            const off = n * Ns, po = n * Nt;
            for (let t = 0; t < Nt; t++) {
                const j = i0[t], w = wg[t];
                p_v[po + t] = F_mod[off + j] * (1 - w) + F_mod[off + j + 1] * w;
            }
        }

        const dt = 1 / fs;
        const vnT = vn.subarray(0, Nt);
        for (let n = 0; n < Nmod; n++) {
            const w = 2 * Math.PI * f[n];
            const co = njCoeff(w, xi, dt);
            const p = p_v.subarray(n * Nt, (n + 1) * Nt);
            const qn = q.subarray(n * Nt, (n + 1) * Nt);
            LA.filter3(co.bq, co.a, p, qn);
            LA.filter3(co.bv, co.a, p, vnT);
            const qddn = qdd.subarray(n * Nt, (n + 1) * Nt);
            for (let t = 0; t < Nt; t++) qddn[t] = p[t] - 2 * xi * w * vn[t] - w * w * qn[t];
            if (qc) {
                const qcn = qc.subarray(n * Nt, (n + 1) * Nt);
                for (let t = 0; t < Nt; t++) qcn[t] = qn[t] - p[t] / w2[n];
            }
        }

        const Q = qc || q;
        res.amax[vi] = absMax(peakOf(Psi_d, nx, Nmod, qdd, Nt, null, null, null, 1));
        if (!secForce) {
            res.dmax[vi] = absMax(peakOf(Psi_d, nx, Nmod, q, Nt, null, null, null, 1));
            res.rmax[vi] = absMax(peakOf(Psi_r, nx, Nmod, q, Nt, null, null, null, 1));
        } else if (statCorr) {
            res.dmax[vi] = absMax(peakOf(Psi_d, nx, Nmod, Q, Nt, d_st, i0, wg, 1));
            res.rmax[vi] = absMax(peakOf(Psi_r, nx, Nmod, Q, Nt, r_st, i0, wg, 1));
            const rM_ = peakOf(Psi_M, nF, Nmod, Q, Nt, M_st, i0, wg, 1);
            res.Mmax[vi] = rM_[0]; res.Mmin[vi] = rM_[1];
            res.Vmax[vi] = absMax(peakOf(Psi_V, nV, Nmod, Q, Nt, V_st, i0, wg, 1));
            if (hasF) {
                const r = peakOf(Psi_rF, nsF, Nmod, Q, Nt, rF_st, i0, wg, 1);
                res.RFmax[vi] = r[0]; res.RFmin[vi] = r[1];
            }
            if (hasM) res.RMmax[vi] = absMax(peakOf(Psi_rM, nsM, Nmod, Q, Nt, rM_st, i0, wg, 1));
        } else {
            res.dmax[vi] = absMax(peakOf(Psi_d, nx, Nmod, q, Nt, null, null, null, 1));
            res.rmax[vi] = absMax(peakOf(Psi_r, nx, Nmod, q, Nt, null, null, null, 1));
            const rM_ = peakOf(Psi_M, nF, Nmod, q, Nt, null, null, null, 1);
            res.Mmax[vi] = rM_[0]; res.Mmin[vi] = rM_[1];
            res.Vmax[vi] = absMax(peakOf(Psi_V, nV, Nmod, q, Nt, null, null, null, 1));
            if (hasF) {
                const r = peakOf(Psi_rF, nsF, Nmod, q, Nt, apF, i0, wg, -1);
                res.RFmax[vi] = r[0]; res.RFmin[vi] = r[1];
            }
            if (hasM) res.RMmax[vi] = absMax(peakOf(Psi_rM, nsM, Nmod, q, Nt, apM, i0, wg, -1));
        }
        if (onProgress) onProgress(vi + 1, nv);
    }
    if (!hasM) res.RMmax = null;
    return res;
};

// Full run: every train of the set over every speed, serial.  The GUI splits
// this over workers instead, but this is the reference path.
BD.run = function (inp) {
    const model = BD.buildModel(inp);
    const trains = inp.F_train;
    const nv = inp.v.length;
    const out = { f: model.f, x: model.x, xF: model.xF, v: inp.v, model };
    const flds = ['dmax', 'amax', 'rmax', 'Mmax', 'Mmin', 'Vmax', 'RFmax', 'RFmin', 'RMmax'];
    for (const fl of flds) out[fl] = [];
    for (let k = 0; k < trains.length; k++) {
        const r = BD.runTrain(model, inp, trains[k], inp.v, inp.v);
        for (const fl of flds) out[fl].push(r[fl]);
    }
    if (out.RMmax[0] === null) out.RMmax = [];
    if (inp.secForces === false) { out.Mmax = []; out.Mmin = []; out.Vmax = []; out.RFmax = []; out.RFmin = []; out.RMmax = []; }
    return out;
};

if (typeof module !== 'undefined' && module.exports) module.exports = BD;
root.BD = BD;
})(typeof self !== 'undefined' ? self : this);

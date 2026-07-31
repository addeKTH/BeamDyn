/* beamstat2d.js -------------------------------------------------------------
   Static response of a 2D beam to passing trains, via influence lines.
   Port of BeamStat2D.m - the static companion to beamdyn2d.js: same model,
   same sign conventions, same peak definitions, no speed axis.

   Two stages, as in the MATLAB version:
     1. INFLUENCE LINES.  A unit downward load at every node in turn, all
        solved together against one Cholesky factorisation of K.  The recovery
        operators then give an influence surface IL(load position, output point)
        for each response quantity.
     2. TRAIN RESPONSE.  Moving a set of point loads over an influence line is
        a convolution of the axle load pattern with the line, so the whole
        envelope is one convolution per quantity - no loop over positions.
        Each axle is split linearly between its two neighbouring nodes, so no
        position is snapped to the grid.

   The peak is the moving-train quasi-static value, i.e. exactly what the
   dynamic run is compared against by dynamicfactor.js.
---------------------------------------------------------------------------*/
(function (root) {
'use strict';

const req = (typeof require === 'function' && typeof module !== 'undefined');
const LA = req ? require('./linalg.js') : root.LA;
const BM = req ? require('./beammodel.js') : root.BM;

const BS = {};

// Influence surfaces, stored column major: IL[c*nL + r], r = load position
// (node), c = output point.
BS.influence = function (inp) {
    const g = BM.geometry(inp);
    const { dL, EL, ndof, sumL } = g;
    const bw = BM.BW;
    const { K } = BM.assemble(g, inp.EI, 0);
    const bcd = BM.bcDofs(g, inp.bc);
    const Kf = BM.reduce(K, ndof, bw, bcd.f2g);
    const R = LA.bandChol(Kf, bcd.nf, bw);

    const xn = Float64Array.from({ length: EL + 1 }, (_, i) => i * dL);
    const Rd = BM.recov(xn, dL, EL, 0);
    const Rr = BM.recov(xn, dL, EL, 1);
    const RM = BM.recov(xn, dL, EL, 2); RM.scale = inp.EI;
    const RV = BM.recov(xn, dL, EL, 3); RV.scale = inp.EI;

    const setF = new Set(), setM = new Set();
    for (const [pos, dir] of inp.bc) {
        const d = Math.round(2 * pos / dL) + dir - 1;
        (dir === 1 ? setF : setM).add(d);
    }
    const pdofF = Int32Array.from(Array.from(setF).sort((a, b) => a - b));
    const pdofM = Int32Array.from(Array.from(setM).sort((a, b) => a - b));

    const nL = EL + 1, nC = EL + 1;
    const nsF = pdofF.length, nsM = pdofM.length;
    const IL = {
        s: xn, x: xn, nL, nC, nsF, nsM,
        d: new Float64Array(nL * nC), r: new Float64Array(nL * nC),
        M: new Float64Array(nL * nC), V: new Float64Array(nL * nC),
        RF: new Float64Array(nL * nsF), RM: nsM ? new Float64Array(nL * nsM) : null
    };

    const bf = new Float64Array(bcd.nf), uf = new Float64Array(bcd.nf);
    const u = new Float64Array(ndof);
    const buf = new Float64Array(nC);

    for (let r = 0; r < nL; r++) {                   // unit DOWNWARD load at node r
        const dofw = 2 * r;
        bf.fill(0);
        if (bcd.g2f[dofw] >= 0) bf[bcd.g2f[dofw]] = -1;
        LA.bandCholSolve(R, bcd.nf, bw, bf, uf);
        u.fill(0);
        for (let j = 0; j < bcd.nf; j++) u[bcd.f2g[j]] = uf[j];

        BM.applyRecov(Rd, u, buf); for (let c = 0; c < nC; c++) IL.d[c * nL + r] = buf[c];
        BM.applyRecov(Rr, u, buf); for (let c = 0; c < nC; c++) IL.r[c * nL + r] = buf[c];
        BM.applyRecov(RM, u, buf); for (let c = 0; c < nC; c++) IL.M[c * nL + r] = buf[c];
        BM.applyRecov(RV, u, buf); for (let c = 0; c < nC; c++) IL.V[c * nL + r] = buf[c];
        for (let k = 0; k < nsF; k++) {
            const applied = (pdofF[k] === dofw) ? -1 : 0;
            IL.RF[k * nL + r] = LA.bandRowDot(K, ndof, bw, pdofF[k], u) - applied;
        }
        for (let k = 0; k < nsM; k++) {
            const applied = (pdofM[k] === dofw) ? -1 : 0;
            IL.RM[k * nL + r] = LA.bandRowDot(K, ndof, bw, pdofM[k], u) - applied;
        }
    }
    // tidy round-off in the section forces and reactions, as in the m-file
    const clip = (A, tol) => { if (A) for (let i = 0; i < A.length; i++) if (Math.abs(A[i]) < tol) A[i] = 0; };
    clip(IL.M, 1e-9); clip(IL.V, 1e-9); clip(IL.RF, 1e-8); clip(IL.RM, 1e-8);
    return IL;
};

// Axle load pattern on the node grid: each axle contributes to the two
// bracketing nodes with linear weights - exact placement, no snapping.
function axlePattern(F, dL) {
    const na = F.x.length;
    let nmax = 0;
    for (let a = 0; a < na; a++) nmax = Math.max(nmax, Math.floor(F.x[a] / dL) + 2);
    const ell = new Float64Array(nmax);
    for (let a = 0; a < na; a++) {
        const gp = F.x[a] / dL;
        const g0 = Math.floor(gp), fr = gp - g0;
        ell[g0] += F.P[a] * (1 - fr);
        ell[g0 + 1] += F.P[a] * fr;
    }
    return ell;
}

// Convolve every column of A [nL x nC] (column major) with ell and reduce the
// result to (max, min).  The convolution output is never stored in full.
function convPeak(A, nL, nC, ell) {
    if (!A || nC === 0) return null;
    const ne = ell.length, nR = nL + ne - 1;
    let mx = -Infinity, mn = Infinity;
    for (let c = 0; c < nC; c++) {
        const off = c * nL;
        for (let k = 0; k < nR; k++) {
            let s = 0;
            const m0 = Math.max(0, k - nL + 1), m1 = Math.min(ne - 1, k);
            for (let m = m0; m <= m1; m++) s += A[off + (k - m)] * ell[m];
            if (s > mx) mx = s;
            if (s < mn) mn = s;
        }
    }
    return [mx, mn];
}

const absMax = (r) => (r ? Math.max(Math.abs(r[0]), Math.abs(r[1])) : 0);

// Peak static response of one train, all quantities.
BS.runTrain = function (IL, inp, train) {
    const F = BM.prepareTrain(train, inp.W || 0, inp.dL);
    const ell = axlePattern(F, inp.dL);
    const { nL, nC, nsF, nsM } = IL;
    const rd = convPeak(IL.d, nL, nC, ell);
    const rr = convPeak(IL.r, nL, nC, ell);
    const out = {
        dmax: absMax(rd), rmax: absMax(rr),
        Mmax: 0, Mmin: 0, Vmax: 0, RFmax: 0, RFmin: 0, RMmax: null
    };
    if (inp.secForces !== false) {
        const rM = convPeak(IL.M, nL, nC, ell);
        const rV = convPeak(IL.V, nL, nC, ell);
        const rF = convPeak(IL.RF, nL, nsF, ell);
        out.Mmax = rM[0]; out.Mmin = rM[1]; out.Vmax = absMax(rV);
        if (rF) { out.RFmax = rF[0]; out.RFmin = rF[1]; }
        if (nsM) out.RMmax = absMax(convPeak(IL.RM, nL, nsM, ell));
    }
    return out;
};

// Every train of a set.
BS.run = function (inp) {
    const IL = BS.influence(inp);
    const flds = ['dmax', 'rmax', 'Mmax', 'Mmin', 'Vmax', 'RFmax', 'RFmin', 'RMmax'];
    const out = {}; for (const f of flds) out[f] = [];
    for (const tr of inp.F_train) {
        const r = BS.runTrain(IL, inp, tr);
        for (const f of flds) out[f].push(r[f]);
    }
    out.IL = IL;
    return out;
};

if (typeof module !== 'undefined' && module.exports) module.exports = BS;
root.BS = BS;
})(typeof self !== 'undefined' ? self : this);

/* beammodel.js --------------------------------------------------------------
   The model block that BeamDyn2D.m and BeamStat2D.m share, factored out once:
   element matrices, assembly into a symmetric band, boundary conditions, the
   Hermite recovery operators and the axle-pattern preparation.

   Bernoulli-Euler beam elements, consistent mass, uniform section.
   Dof order (w1,fi1,w2,fi2,...), so the dof at position x in direction dir
   (1 = vertical, 2 = rotation) is 2*x/dL + dir, one based.

   COORDINATE SYSTEM: x runs along the beam to the right, w is vertical and
   positive upwards, and the axle loads act downwards.  A sagging moment is
   therefore positive.
---------------------------------------------------------------------------*/
(function (root) {
'use strict';

const LA = (typeof require === 'function' && typeof module !== 'undefined')
    ? require('./linalg.js') : root.LA;

const BM = {};
BM.BW = 3;                                          // half bandwidth

// ------------------------------------------------------------------ geometry
BM.geometry = function (inp) {
    const dL = inp.dL;
    const L = inp.L.map(v => dL * Math.round(v / dL));   // spans on the grid
    const sumL = L.reduce((a, b) => a + b, 0);
    const EL = Math.round(sumL / dL);
    if (!(EL >= 1)) throw new Error('the beam has no elements - check L and dL');
    return { dL, L, sumL, EL, ndof: 2 * (EL + 1) };
};

// -------------------------------------------------------------- assembly
// K and M as symmetric bands over all dofs.  Every element is identical, so
// the local matrices are formed once.
BM.assemble = function (g, EI, m) {
    const { dL, EL, ndof } = g, bw = BM.BW;
    const k1 = 12 * EI / (dL * dL * dL), k2 = 6 * EI / (dL * dL), k3 = 4 * EI / dL;
    const ke = [
        [k1, k2, -k1, k2],
        [k2, k3, -k2, k3 / 2],
        [-k1, -k2, k1, -k2],
        [k2, k3 / 2, -k2, k3]
    ];
    const c = m * dL / 420;
    const me = [
        [156 * c, 22 * dL * c, 54 * c, -13 * dL * c],
        [22 * dL * c, 4 * dL * dL * c, 13 * dL * c, -3 * dL * dL * c],
        [54 * c, 13 * dL * c, 156 * c, -22 * dL * c],
        [-13 * dL * c, -3 * dL * dL * c, -22 * dL * c, 4 * dL * dL * c]
    ];
    const K = LA.bandAlloc(ndof, bw), M = LA.bandAlloc(ndof, bw);
    for (let e = 0; e < EL; e++) {
        const o = 2 * e;                            // first dof of the element
        for (let a = 0; a < 4; a++)
            for (let b = a; b < 4; b++) {
                LA.bandAdd(K, bw, o + a, o + b, ke[a][b]);
                if (m) LA.bandAdd(M, bw, o + a, o + b, me[a][b]);
            }
    }
    return { K, M, bw };
};

// ------------------------------------------------------ boundary conditions
// bc = [[pos, dir], ...].  Returns the fixed and free dof lists (0 based) and
// the mapping used to compact the free system.
BM.bcDofs = function (g, bc) {
    const { dL, sumL, ndof } = g;
    const fixed = new Set();
    for (const [pos, dir] of bc) {
        if (pos < -1e-6 || pos > sumL + 1e-6) throw new Error('support outside the beam');
        if (dir !== 1 && dir !== 2) throw new Error('support direction must be 1 or 2');
        const nd = pos / dL;
        if (Math.abs(nd - Math.round(nd)) > 1e-6)
            throw new Error('support at x = ' + pos + ' is not on a node (dL = ' + dL + ')');
        fixed.add(Math.round(2 * pos / dL) + dir - 1);   // 0 based
    }
    const pdof = Array.from(fixed).sort((a, b) => a - b);
    const g2f = new Int32Array(ndof).fill(-1);
    const f2g = [];
    for (let i = 0; i < ndof; i++) if (!fixed.has(i)) { g2f[i] = f2g.length; f2g.push(i); }
    return { pdof, f2g: Int32Array.from(f2g), g2f, nf: f2g.length };
};

// Compact a global band matrix onto the free dofs.  Deleting rows and columns
// cannot increase the bandwidth, so the band layout is kept.
BM.reduce = function (A, ndof, bw, f2g) {
    const nf = f2g.length;
    const Af = LA.bandAlloc(nf, bw);
    for (let i = 0; i < nf; i++) {
        const gi = f2g[i];
        for (let j = i; j < Math.min(i + bw + 1, nf); j++) {
            const gj = f2g[j];
            if (gj - gi > bw) break;
            Af[i * (bw + 1) + (j - i)] = LA.bandGet(A, bw, gi, gj);
        }
    }
    return Af;
};

// ------------------------------------------------------- recovery operators
// A compact [np x ndof] operator mapping the dof vector to
//   der = 0  displacement w      der = 1  rotation w'
//   der = 2  curvature w''       der = 3  w'''
// at the points x, from the analytical derivatives of the Hermite shape
// functions.  Optional e/sl give the owning element and the local coordinate
// explicitly, which is how a point on a node is assigned to one side.
// Stored as, per point, the first dof index and four coefficients.
BM.recov = function (x, dL, EL, der, eIn, slIn) {
    const np = x.length;
    const dof = new Int32Array(np);
    const N = new Float64Array(4 * np);
    for (let p = 0; p < np; p++) {
        let e, sl;
        if (eIn) { e = eIn[p]; sl = slIn[p]; }
        else {
            e = Math.min(Math.max(Math.floor(x[p] / dL), 0), EL - 1);
            sl = x[p] / dL - e;
        }
        const s2 = sl * sl, s3 = sl * s2;
        let n0, n1, n2, n3;
        switch (der) {
            case 0:
                n0 = 1 - 3 * s2 + 2 * s3; n1 = dL * (sl - 2 * s2 + s3);
                n2 = 3 * s2 - 2 * s3; n3 = dL * (s3 - s2); break;
            case 1:
                n0 = (-6 * sl + 6 * s2) / dL; n1 = 1 - 4 * sl + 3 * s2;
                n2 = (6 * sl - 6 * s2) / dL; n3 = 3 * s2 - 2 * sl; break;
            case 2:
                n0 = (-6 + 12 * sl) / (dL * dL); n1 = (-4 + 6 * sl) / dL;
                n2 = (6 - 12 * sl) / (dL * dL); n3 = (-2 + 6 * sl) / dL; break;
            case 3:
                n0 = 12 / (dL * dL * dL); n1 = 6 / (dL * dL);
                n2 = -12 / (dL * dL * dL); n3 = 6 / (dL * dL); break;
        }
        dof[p] = 2 * e;
        N[4 * p] = n0; N[4 * p + 1] = n1; N[4 * p + 2] = n2; N[4 * p + 3] = n3;
    }
    return { np, dof, N, scale: 1 };
};

// Apply a recovery operator to one dof vector.
BM.applyRecov = function (R, u, out) {
    const np = R.np, dof = R.dof, N = R.N, sc = R.scale;
    out = out || new Float64Array(np);
    for (let p = 0; p < np; p++) {
        const d = dof[p], q = 4 * p;
        out[p] = sc * (N[q] * u[d] + N[q + 1] * u[d + 1] + N[q + 2] * u[d + 2] + N[q + 3] * u[d + 3]);
    }
    return out;
};

// Apply a recovery operator to every column of Phi (ndof x nc, column major).
// Returns [np x nc] row major, which is the order the response kernels read.
BM.recovModal = function (R, Phi, ndof, nc) {
    const np = R.np, out = new Float64Array(np * nc);
    for (let c = 0; c < nc; c++) {
        const col = Phi.subarray(c * ndof, (c + 1) * ndof);
        for (let p = 0; p < np; p++) {
            const d = R.dof[p], q = 4 * p;
            out[p * nc + c] = R.scale * (R.N[q] * col[d] + R.N[q + 1] * col[d + 1] +
                                         R.N[q + 2] * col[d + 2] + R.N[q + 3] * col[d + 3]);
        }
    }
    return out;
};

// ------------------------------------------------------ mode shape at any x
// Within an element the FE mode shape IS the element cubic, so this is exact
// - no spline fit through the nodal values.
BM.hermiteEval = function (Phi, ndof, nm, dL, EL, x, out) {
    const nx = x.length;
    out = out || new Float64Array(nm * nx);         // [nm x nx] row major
    for (let p = 0; p < nx; p++) {
        const e = Math.min(Math.max(Math.floor(x[p] / dL), 0), EL - 1);
        const s = x[p] / dL - e, s2 = s * s, s3 = s * s2;
        const N1 = 1 - 3 * s2 + 2 * s3, N2 = dL * (s - 2 * s2 + s3);
        const N3 = 3 * s2 - 2 * s3, N4 = dL * (s3 - s2);
        const i1 = 2 * e;
        for (let n = 0; n < nm; n++) {
            const c = n * ndof;
            out[n * nx + p] = N1 * Phi[c + i1] + N2 * Phi[c + i1 + 1] +
                              N3 * Phi[c + i1 + 2] + N4 * Phi[c + i1 + 3];
        }
    }
    return out;
};

// --------------------------------------------------------- axle preparation
// Optional triangular load distribution over the width W, then merge loads
// that land on the same position, as in both MATLAB solvers.
BM.triangWin = function (N) {
    const w = new Float64Array(N);
    for (let k = 1; k <= N; k++) {
        w[k - 1] = (N % 2)
            ? Math.min(2 * k / (N + 1), 2 * (N + 1 - k) / (N + 1))
            : Math.min((2 * k - 1) / N, (2 * (N + 1 - k) - 1) / N);
    }
    return w;
};

BM.prepareTrain = function (train, W, dL) {
    let xs = Array.from(train.x), Ps = Array.from(train.P);
    if (W > 0) {
        const nW = Math.floor(W / dL) + 1;          // -W/2 : dL : W/2
        const xi = [];
        for (let i = 0; i < nW; i++) xi.push(-W / 2 + i * dL);
        const wi = BM.triangWin(nW);
        let sw = 0; for (const v of wi) sw += v;
        const xj = [], Pj = [];
        for (let a = 0; a < xs.length; a++)
            for (let i = 0; i < nW; i++) { xj.push(xs[a] + xi[i]); Pj.push(Ps[a] * wi[i] / sw); }
        const xmin = Math.min(...xj);
        xs = xj.map(v => v - xmin); Ps = Pj;
    }
    // round to 1e-6 and merge coincident positions, ascending
    const map = new Map();
    for (let i = 0; i < xs.length; i++) {
        const key = Math.round(xs[i] * 1e6) / 1e6;
        map.set(key, (map.get(key) || 0) + Ps[i]);
    }
    const keys = Array.from(map.keys()).sort((a, b) => a - b);
    return { x: Float64Array.from(keys), P: Float64Array.from(keys, k => map.get(k)) };
};

if (typeof module !== 'undefined' && module.exports) module.exports = BM;
root.BM = BM;
})(typeof self !== 'undefined' ? self : this);

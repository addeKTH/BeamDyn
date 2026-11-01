/* linalg.js -----------------------------------------------------------------
   Linear algebra for the beam models: symmetric band storage, banded
   Cholesky, a dense symmetric eigensolver and the IIR filter of the
   Nigam-Jennings recurrence.

   The beam dofs are ordered (w1,fi1,w2,fi2,...), so an element couples four
   consecutive dofs and K and M have half bandwidth 3.  The static solves use
   the banded Cholesky factor, O(n*bw) per solve.

   Band storage:  A[i*(bw+1) + (j-i)] = a(i,j)  for  i <= j <= i+bw
   (upper half only; the matrix is symmetric).

   The generalized eigenvalue problem K*phi = lam*M*phi is reduced to a
   standard one with the Cholesky factor of M and solved densely
   (Householder tridiagonalisation + implicit QL, after JAMA / EISPACK tred2
   and tql2), which returns the whole spectrum as MATLAB's eig() does.

   Andreas Andersson, 2026-11-01 (assisted by Claude)
---------------------------------------------------------------------------*/
(function (root) {
'use strict';

const LA = {};

// ---------------------------------------------------------------- band basics
LA.bandAlloc = (n, bw) => new Float64Array(n * (bw + 1));

LA.bandGet = function (A, bw, i, j) {
    if (j < i) { const t = i; i = j; j = t; }
    if (j - i > bw) return 0;
    return A[i * (bw + 1) + (j - i)];
};

LA.bandAdd = function (A, bw, i, j, v) {
    if (j < i) { const t = i; i = j; j = t; }
    A[i * (bw + 1) + (j - i)] += v;
};

// One row of a symmetric band matrix times a full vector
LA.bandRowDot = function (A, n, bw, i, x) {
    let s = 0;
    const j0 = Math.max(0, i - bw), j1 = Math.min(i + bw, n - 1);
    for (let j = j0; j <= j1; j++) s += LA.bandGet(A, bw, i, j) * x[j];
    return s;
};

// Band matrix restricted to the dofs in idx (ascending); deleting rows and
// columns cannot increase the bandwidth
LA.bandReduce = function (A, bw, idx) {
    const nf = idx.length, Af = LA.bandAlloc(nf, bw);
    for (let i = 0; i < nf; i++)
        for (let j = i; j < Math.min(i + bw + 1, nf); j++) {
            if (idx[j] - idx[i] > bw) break;
            Af[i * (bw + 1) + (j - i)] = LA.bandGet(A, bw, idx[i], idx[j]);
        }
    return Af;
};

// Dense row-major copy of a symmetric band matrix
LA.bandToDense = function (A, n, bw) {
    const D = new Float64Array(n * n);
    for (let i = 0; i < n; i++)
        for (let j = i; j <= Math.min(i + bw, n - 1); j++) {
            const a = A[i * (bw + 1) + (j - i)];
            D[i * n + j] = a; D[j * n + i] = a;
        }
    return D;
};

// --------------------------------------------------------- banded Cholesky
// A = R'R with R upper triangular, returned in the same band layout
LA.bandChol = function (A, n, bw) {
    const R = new Float64Array(A);
    const w = bw + 1;
    for (let i = 0; i < n; i++) {
        const jm = Math.min(i + bw, n - 1);
        for (let j = i; j <= jm; j++) {
            let s = R[i * w + (j - i)];
            for (let k = Math.max(0, j - bw); k < i; k++) s -= R[k * w + (i - k)] * R[k * w + (j - k)];
            if (j === i) {
                if (!(s > 0)) throw new Error('the stiffness matrix is singular (check the supports)');
                R[i * w] = Math.sqrt(s);
            } else {
                R[i * w + (j - i)] = s / R[i * w];
            }
        }
    }
    return R;
};

// Solve A x = b with the band factor R (x may be b)
LA.bandCholSolve = function (R, n, bw, b, x) {
    x = x || new Float64Array(n);
    if (x !== b) x.set(b);
    const w = bw + 1;
    for (let i = 0; i < n; i++) {                       // R' y = b
        let s = x[i];
        for (let k = Math.max(0, i - bw); k < i; k++) s -= R[k * w + (i - k)] * x[k];
        x[i] = s / R[i * w];
    }
    for (let i = n - 1; i >= 0; i--) {                  // R x = y
        let s = x[i];
        const jm = Math.min(i + bw, n - 1);
        for (let j = i + 1; j <= jm; j++) s -= R[i * w + (j - i)] * x[j];
        x[i] = s / R[i * w];
    }
    return x;
};

// ------------------------------------------- dense symmetric eigenproblem
// Eigenvalues (ascending) and eigenvectors (columns of V, row major) of the
// dense symmetric [n x n] A, by Householder tridiagonalisation (tred2) and
// the implicit QL algorithm (tql2), as in JAMA.
LA.symEig = function (A, n) {
    const V = Float64Array.from(A), d = new Float64Array(n), e = new Float64Array(n);
    const at = (i, j) => i * n + j;

    // tred2
    for (let j = 0; j < n; j++) d[j] = V[at(n - 1, j)];
    for (let i = n - 1; i > 0; i--) {
        let scale = 0, h = 0;
        for (let k = 0; k < i; k++) scale += Math.abs(d[k]);
        if (scale === 0) {
            e[i] = d[i - 1];
            for (let j = 0; j < i; j++) { d[j] = V[at(i - 1, j)]; V[at(i, j)] = 0; V[at(j, i)] = 0; }
        } else {
            for (let k = 0; k < i; k++) { d[k] /= scale; h += d[k] * d[k]; }
            let f = d[i - 1], g = Math.sqrt(h);
            if (f > 0) g = -g;
            e[i] = scale * g; h -= f * g; d[i - 1] = f - g;
            for (let j = 0; j < i; j++) e[j] = 0;
            for (let j = 0; j < i; j++) {
                f = d[j]; V[at(j, i)] = f; g = e[j] + V[at(j, j)] * f;
                for (let k = j + 1; k <= i - 1; k++) { g += V[at(k, j)] * d[k]; e[k] += V[at(k, j)] * f; }
                e[j] = g;
            }
            f = 0;
            for (let j = 0; j < i; j++) { e[j] /= h; f += e[j] * d[j]; }
            const hh = f / (h + h);
            for (let j = 0; j < i; j++) e[j] -= hh * d[j];
            for (let j = 0; j < i; j++) {
                f = d[j]; g = e[j];
                for (let k = j; k <= i - 1; k++) V[at(k, j)] -= (f * e[k] + g * d[k]);
                d[j] = V[at(i - 1, j)]; V[at(i, j)] = 0;
            }
        }
        d[i] = h;
    }
    for (let i = 0; i < n - 1; i++) {
        V[at(n - 1, i)] = V[at(i, i)]; V[at(i, i)] = 1;
        const h = d[i + 1];
        if (h !== 0) {
            for (let k = 0; k <= i; k++) d[k] = V[at(k, i + 1)] / h;
            for (let j = 0; j <= i; j++) {
                let g = 0;
                for (let k = 0; k <= i; k++) g += V[at(k, i + 1)] * V[at(k, j)];
                for (let k = 0; k <= i; k++) V[at(k, j)] -= g * d[k];
            }
        }
        for (let k = 0; k <= i; k++) V[at(k, i + 1)] = 0;
    }
    for (let j = 0; j < n; j++) { d[j] = V[at(n - 1, j)]; V[at(n - 1, j)] = 0; }
    V[at(n - 1, n - 1)] = 1; e[0] = 0;

    // tql2
    for (let i = 1; i < n; i++) e[i - 1] = e[i];
    e[n - 1] = 0;
    let f = 0, tst1 = 0;
    const eps = Math.pow(2, -52);
    for (let l = 0; l < n; l++) {
        tst1 = Math.max(tst1, Math.abs(d[l]) + Math.abs(e[l]));
        let m = l;
        while (m < n) { if (Math.abs(e[m]) <= eps * tst1) break; m++; }
        if (m > l) {
            let iter = 0;
            do {
                if (++iter > 100) throw new Error('symEig: no convergence');
                let g = d[l];
                let p = (d[l + 1] - g) / (2 * e[l]);
                let r = Math.hypot(p, 1);
                if (p < 0) r = -r;
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                const dl1 = d[l + 1];
                let h = g - d[l];
                for (let i = l + 2; i < n; i++) d[i] -= h;
                f += h;
                p = d[m];
                let c = 1, c2 = c, c3 = c, s = 0, s2 = 0;
                const el1 = e[l + 1];
                for (let i = m - 1; i >= l; i--) {
                    c3 = c2; c2 = c; s2 = s;
                    g = c * e[i]; h = c * p;
                    r = Math.hypot(p, e[i]);
                    e[i + 1] = s * r; s = e[i] / r; c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    for (let k = 0; k < n; k++) {
                        h = V[at(k, i + 1)];
                        V[at(k, i + 1)] = s * V[at(k, i)] + c * h;
                        V[at(k, i)] = c * V[at(k, i)] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p; d[l] = c * p;
            } while (Math.abs(e[l]) > eps * tst1);
        }
        d[l] += f; e[l] = 0;
    }
    // ascending order
    for (let i = 0; i < n - 1; i++) {
        let k = i, p = d[i];
        for (let j = i + 1; j < n; j++) if (d[j] < p) { k = j; p = d[j]; }
        if (k !== i) {
            d[k] = d[i]; d[i] = p;
            for (let j = 0; j < n; j++) { const t = V[at(j, i)]; V[at(j, i)] = V[at(j, k)]; V[at(j, k)] = t; }
        }
    }
    return { lam: d, V };
};

// Generalized problem K y = lam M y (dense row major, M positive definite):
// M = L L', C = L^-1 K L^-T, C z = lam z, y = L^-T z.  Returns ascending lam
// and Y (row major, columns M-orthonormal).
LA.genEig = function (K, M, n) {
    const L = new Float64Array(n * n);
    for (let i = 0; i < n; i++)
        for (let j = 0; j <= i; j++) {
            let s = M[i * n + j];
            for (let k = 0; k < j; k++) s -= L[i * n + k] * L[j * n + k];
            if (i === j) {
                if (!(s > 0)) throw new Error('the mass matrix is not positive definite');
                L[i * n + i] = Math.sqrt(s);
            } else L[i * n + j] = s / L[j * n + j];
        }
    const T = new Float64Array(n * n);                  // T = L^-1 K
    for (let c = 0; c < n; c++)
        for (let i = 0; i < n; i++) {
            let s = K[i * n + c];
            for (let k = 0; k < i; k++) s -= L[i * n + k] * T[k * n + c];
            T[i * n + c] = s / L[i * n + i];
        }
    const C = new Float64Array(n * n);                  // C = T L^-T
    for (let r = 0; r < n; r++)
        for (let i = 0; i < n; i++) {
            let s = T[r * n + i];
            for (let k = 0; k < i; k++) s -= L[i * n + k] * C[r * n + k];
            C[r * n + i] = s / L[i * n + i];
        }
    for (let i = 0; i < n; i++)                          // symmetrise round-off
        for (let j = i + 1; j < n; j++) {
            const m = 0.5 * (C[i * n + j] + C[j * n + i]);
            C[i * n + j] = m; C[j * n + i] = m;
        }
    const { lam, V } = LA.symEig(C, n);
    const Y = new Float64Array(n * n);                  // Y = L^-T Z
    for (let c = 0; c < n; c++)
        for (let i = n - 1; i >= 0; i--) {
            let s = V[i * n + c];
            for (let k = i + 1; k < n; k++) s -= L[k * n + i] * Y[k * n + c];
            Y[i * n + c] = s / L[i * n + i];
        }
    return { lam, Y };
};

// ------------------------------------------------------------------ IIR filter
// MATLAB filter(b,a,x) with zero initial conditions, a(1) = 1, three taps
// (direct form II transposed); y may be a preallocated buffer
LA.filter3 = function (b, a, x, n, y) {
    let z1 = 0, z2 = 0;
    const b0 = b[0], b1 = b[1], b2 = b[2], a1 = a[1], a2 = a[2];
    for (let k = 0; k < n; k++) {
        const xk = x[k];
        const yk = b0 * xk + z1;
        z1 = b1 * xk - a1 * yk + z2;
        z2 = b2 * xk - a2 * yk;
        y[k] = yk;
    }
    return y;
};

if (typeof module !== 'undefined' && module.exports) module.exports = LA;
root.LA = LA;
})(typeof self !== 'undefined' ? self : this);

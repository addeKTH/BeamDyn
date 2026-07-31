/* linalg.js -----------------------------------------------------------------
   Banded symmetric linear algebra for the beam models.

   The beam dofs are ordered (w1,fi1,w2,fi2,...), so an element couples four
   consecutive dofs and both K and M have half bandwidth 3.  Everything here
   works on that packed band, which is what makes a fine mesh affordable:
   the Cholesky factorisation is O(n*bw^2) and each solve O(n*bw), against
   O(n^3) and O(n^2) for the dense equivalents.

   Band storage:  A[i*(bw+1) + (j-i)] = a(i,j)  for  i <= j <= i+bw.
   Only the upper half is stored; the matrix is symmetric.

   Eigenvalues use subspace iteration with the banded factorisation as the
   inner solve, so only the lowest few modes are computed - not the whole
   spectrum, as a dense eig() would.
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

// Symmetric band times full vector.
LA.bandMulVec = function (A, n, bw, x, y) {
    y = y || new Float64Array(n);
    y.fill(0);
    for (let i = 0; i < n; i++) {
        const r = i * (bw + 1);
        y[i] += A[r] * x[i];
        const jm = Math.min(i + bw, n - 1);
        for (let j = i + 1; j <= jm; j++) {
            const a = A[r + (j - i)];
            if (a === 0) continue;
            y[i] += a * x[j];
            y[j] += a * x[i];
        }
    }
    return y;
};

// One row of a symmetric band matrix, as a dot product with a full vector.
LA.bandRowDot = function (A, n, bw, i, x) {
    let s = 0;
    const j0 = Math.max(0, i - bw), j1 = Math.min(i + bw, n - 1);
    for (let j = j0; j <= j1; j++) s += LA.bandGet(A, bw, i, j) * x[j];
    return s;
};

// --------------------------------------------------------- banded Cholesky
// A = R'R with R upper triangular, in place in the same band layout.
LA.bandChol = function (A, n, bw) {
    const R = new Float64Array(A);
    const w = bw + 1;
    for (let i = 0; i < n; i++) {
        const jm = Math.min(i + bw, n - 1);
        for (let j = i; j <= jm; j++) {
            let s = R[i * w + (j - i)];
            const k0 = Math.max(0, i - bw);
            for (let k = k0; k < i; k++) {
                const dki = i - k, dkj = j - k;
                if (dkj > bw) continue;
                s -= R[k * w + dki] * R[k * w + dkj];
            }
            if (j === i) {
                if (s <= 0) throw new Error('bandChol: matrix is not positive definite');
                R[i * w + 0] = Math.sqrt(s);
            } else {
                R[i * w + (j - i)] = s / R[i * w];
            }
        }
    }
    return R;
};

// Solve A x = b given the Cholesky band factor R.  b is overwritten if x === b.
LA.bandCholSolve = function (R, n, bw, b, x) {
    x = x || new Float64Array(n);
    if (x !== b) x.set(b);
    const w = bw + 1;
    for (let i = 0; i < n; i++) {                       // R' y = b
        let s = x[i];
        const k0 = Math.max(0, i - bw);
        for (let k = k0; k < i; k++) s -= R[k * w + (i - k)] * x[k];
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

// ------------------------------------------------ small dense symmetric eig
// Cyclic Jacobi on a dense symmetric [n x n] (row major).  Returns ascending
// eigenvalues and the matching eigenvectors as columns of V.
LA.jacobiEig = function (Ain, n) {
    const A = Float64Array.from(Ain);
    const V = new Float64Array(n * n);
    for (let i = 0; i < n; i++) V[i * n + i] = 1;
    for (let sweep = 0; sweep < 100; sweep++) {
        let off = 0;
        for (let p = 0; p < n; p++)
            for (let q = p + 1; q < n; q++) off += A[p * n + q] * A[p * n + q];
        if (off < 1e-32) break;
        for (let p = 0; p < n - 1; p++) {
            for (let q = p + 1; q < n; q++) {
                const apq = A[p * n + q];
                if (Math.abs(apq) < 1e-300) continue;
                const theta = (A[q * n + q] - A[p * n + p]) / (2 * apq);
                const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
                const c = 1 / Math.sqrt(t * t + 1), s = t * c;
                for (let k = 0; k < n; k++) {
                    const akp = A[k * n + p], akq = A[k * n + q];
                    A[k * n + p] = c * akp - s * akq;
                    A[k * n + q] = s * akp + c * akq;
                }
                for (let k = 0; k < n; k++) {
                    const apk = A[p * n + k], aqk = A[q * n + k];
                    A[p * n + k] = c * apk - s * aqk;
                    A[q * n + k] = s * apk + c * aqk;
                }
                for (let k = 0; k < n; k++) {
                    const vkp = V[k * n + p], vkq = V[k * n + q];
                    V[k * n + p] = c * vkp - s * vkq;
                    V[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    const lam = new Float64Array(n);
    for (let i = 0; i < n; i++) lam[i] = A[i * n + i];
    const ix = Array.from({ length: n }, (_, i) => i).sort((a, b) => lam[a] - lam[b]);
    const lamS = new Float64Array(n), VS = new Float64Array(n * n);
    for (let c = 0; c < n; c++) {
        lamS[c] = lam[ix[c]];
        for (let r = 0; r < n; r++) VS[r * n + c] = V[r * n + ix[c]];
    }
    return { lam: lamS, V: VS };
};

// Dense generalized symmetric problem K y = lam M y, both [n x n] row major,
// M positive definite.  Used only for the small projected problem.
LA.denseGenEig = function (K, M, n) {
    // Cholesky M = L L'
    const L = new Float64Array(n * n);
    for (let i = 0; i < n; i++) {
        for (let j = 0; j <= i; j++) {
            let s = M[i * n + j];
            for (let k = 0; k < j; k++) s -= L[i * n + k] * L[j * n + k];
            if (i === j) {
                if (s <= 0) throw new Error('denseGenEig: M not positive definite');
                L[i * n + i] = Math.sqrt(s);
            } else {
                L[i * n + j] = s / L[j * n + j];
            }
        }
    }
    // C = L^-1 K L^-T
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
    const { lam, V } = LA.jacobiEig(C, n);
    // y = L^-T z
    const Y = new Float64Array(n * n);
    for (let c = 0; c < n; c++)
        for (let i = n - 1; i >= 0; i--) {
            let s = V[i * n + c];
            for (let k = i + 1; k < n; k++) s -= L[k * n + i] * Y[k * n + c];
            Y[i * n + c] = s / L[i * n + i];
        }
    return { lam, V: Y };
};

// ------------------------------------------------------- subspace iteration
// The nreq lowest modes of K phi = lam M phi, both symmetric band [n x n],
// K positive definite.  Returns ascending lam and M-orthonormal columns Phi
// (n x nreq, column major), i.e. Phi' M Phi = I.
LA.eigBand = function (K, M, n, bw, nreq, opt) {
    opt = opt || {};
    const tol = opt.tol || 1e-14, maxIt = opt.maxIt || 400;
    nreq = Math.min(nreq, n);
    const q = Math.min(n, Math.max(2 * nreq, nreq + 8));
    const R = LA.bandChol(K, n, bw);

    // Starting vectors, after Bathe: diag(M)/diag(K), then unit vectors at the
    // dofs with the largest mass-to-stiffness ratio, then one alternating sign
    // vector to break symmetry.
    let X = new Float64Array(n * q);                    // column major
    const ratio = new Float64Array(n);
    for (let i = 0; i < n; i++) {
        const ki = K[i * (bw + 1)], mi = M[i * (bw + 1)];
        ratio[i] = ki !== 0 ? mi / ki : 0;
        X[i] = mi;
    }
    const order = Array.from({ length: n }, (_, i) => i).sort((a, b) => ratio[b] - ratio[a]);
    for (let c = 1; c < q - 1 && c - 1 < n; c++) X[c * n + order[c - 1]] = 1;
    for (let i = 0; i < n; i++) X[(q - 1) * n + i] = (i % 2 ? -1 : 1);

    const Xb = new Float64Array(n * q);
    const MX = new Float64Array(n * q);
    const Kr = new Float64Array(q * q), Mr = new Float64Array(q * q);
    const tmp = new Float64Array(n);
    let lamOld = new Float64Array(q).fill(1e300), lam = null, Q = null;

    for (let it = 0; it < maxIt; it++) {
        for (let c = 0; c < q; c++) {                   // MX = M X
            LA.bandMulVec(M, n, bw, X.subarray(c * n, c * n + n), tmp);
            MX.set(tmp, c * n);
        }
        for (let c = 0; c < q; c++)                     // K Xb = M X
            LA.bandCholSolve(R, n, bw, MX.subarray(c * n, c * n + n), Xb.subarray(c * n, c * n + n));
        for (let a = 0; a < q; a++) {                   // Kr = Xb' M X = Xb' K Xb
            for (let b = a; b < q; b++) {
                let s = 0;
                for (let i = 0; i < n; i++) s += Xb[a * n + i] * MX[b * n + i];
                Kr[a * q + b] = s; Kr[b * q + a] = s;
            }
        }
        for (let c = 0; c < q; c++) {                   // MX = M Xb  (reused)
            LA.bandMulVec(M, n, bw, Xb.subarray(c * n, c * n + n), tmp);
            MX.set(tmp, c * n);
        }
        for (let a = 0; a < q; a++) {                   // Mr = Xb' M Xb
            for (let b = a; b < q; b++) {
                let s = 0;
                for (let i = 0; i < n; i++) s += Xb[a * n + i] * MX[b * n + i];
                Mr[a * q + b] = s; Mr[b * q + a] = s;
            }
        }
        const e = LA.denseGenEig(Kr, Mr, q);
        lam = e.lam; Q = e.V;
        const Xn = new Float64Array(n * q);             // X = Xb Q
        for (let c = 0; c < q; c++)
            for (let k = 0; k < q; k++) {
                const qk = Q[k * q + c];
                if (qk === 0) continue;
                for (let i = 0; i < n; i++) Xn[c * n + i] += Xb[k * n + i] * qk;
            }
        X = Xn;
        let done = true;
        for (let c = 0; c < nreq; c++) {
            const d = Math.abs(lam[c] - lamOld[c]) / Math.max(Math.abs(lam[c]), 1e-300);
            if (d > tol) { done = false; break; }
        }
        lamOld = Float64Array.from(lam);
        if (done && it > 1) break;
    }

    // M-normalise and fix the sign (largest component positive), so repeated
    // runs give the same picture.
    const Phi = new Float64Array(n * nreq);
    for (let c = 0; c < nreq; c++) {
        const x = X.subarray(c * n, c * n + n);
        LA.bandMulVec(M, n, bw, x, tmp);
        let s = 0;
        for (let i = 0; i < n; i++) s += x[i] * tmp[i];
        let sc = 1 / Math.sqrt(s);
        let big = 0, bi = 0;
        for (let i = 0; i < n; i++) if (Math.abs(x[i]) > big) { big = Math.abs(x[i]); bi = i; }
        if (x[bi] < 0) sc = -sc;
        for (let i = 0; i < n; i++) Phi[c * n + i] = x[i] * sc;
    }
    return { lam: lam.subarray(0, nreq), Phi };
};

// ------------------------------------------------------------------ IIR filter
// MATLAB filter(b,a,x) with zero initial conditions, a(1) = 1, three taps.
// Direct form II transposed.
LA.filter3 = function (b, a, x, y) {
    const n = x.length;
    y = y || new Float64Array(n);
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

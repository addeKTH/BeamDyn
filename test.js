/* Validation of the BeamDyn2D browser port against analytical solutions. */
const B = require('./beamdyn.js');

let fails = 0;
function check(name, got, want, tol, unit) {
  const rel = Math.abs(got - want) / (Math.abs(want) || 1);
  const ok = rel < tol;
  if (!ok) fails++;
  console.log(`${ok ? '  ok  ' : ' FAIL '} ${name.padEnd(42)} ` +
    `got ${got.toPrecision(10)}  want ${want.toPrecision(10)}  rel ${rel.toExponential(2)}`);
}

console.log('\n=== 1. Simply supported span: natural frequencies ===');
console.log('    f_n = (n^2*pi/(2L^2))*sqrt(EI/m)\n');
{
  const L = 20, m = 15e3, EI = 30e9;
  const inp = { L: [L], m, EI, xi: 1.5, dL: 0.25, Nmod: 4, fmax: Infinity,
                bc: [[0, 1], [L, 1]] };
  const mdl = B.buildModel(inp);
  for (let n = 1; n <= 4; n++) {
    const fa = (n * n * Math.PI / (2 * L * L)) * Math.sqrt(EI / m);
    check(`f_${n}`, mdl.f[n - 1], fa, 1e-4, 'Hz');
  }
}

console.log('\n=== 2. Clamped-clamped span: natural frequencies ===');
console.log('    f_n = (b_n^2/(2*pi*L^2))*sqrt(EI/m), b_n from cosh*cos = 1\n');
{
  const L = 20, m = 15e3, EI = 30e9;
  const inp = { L: [L], m, EI, xi: 1.5, dL: 0.25, Nmod: 3, fmax: Infinity,
                bc: [[0, 1], [0, 2], [L, 1], [L, 2]] };
  const mdl = B.buildModel(inp);
  const beta = [4.730040745, 7.853204624, 10.99560784];   // roots of cosh*cos=1
  for (let n = 1; n <= 3; n++) {
    const fa = (beta[n - 1] ** 2 / (2 * Math.PI * L * L)) * Math.sqrt(EI / m);
    check(`f_${n}`, mdl.f[n - 1], fa, 1e-4, 'Hz');
  }
}

console.log('\n=== 3. Two-span continuous beam: natural frequencies ===');
console.log('    equal spans: modes are the SS and CS modes of one span\n');
{
  const L = 20, m = 15e3, EI = 30e9;
  const inp = { L: [L, L], m, EI, xi: 1.5, dL: 0.25, Nmod: 2, fmax: Infinity,
                bc: [[0, 1], [L, 1], [2 * L, 1]] };
  const mdl = B.buildModel(inp);
  // mode 1 = simply supported half-wave; mode 2 = clamped-pinned (beta=3.9266)
  const f1 = (Math.PI / (2 * L * L)) * Math.sqrt(EI / m);
  const f2 = (3.926602312 ** 2 / (2 * Math.PI * L * L)) * Math.sqrt(EI / m);
  check('f_1 (symmetric)', mdl.f[0], f1, 1e-4);
  check('f_2 (antisymmetric)', mdl.f[1], f2, 1e-4);
}

console.log('\n=== 4. Mass normalisation: Phi\'*M*Phi = I ===\n');
{
  const inp = { L: [20, 20], m: 15e3, EI: 30e9, xi: 1.5, dL: 1, Nmod: 3,
                fmax: Infinity, bc: [[0, 1], [20, 1], [40, 1]] };
  const mdl = B.buildModel(inp);
  const { M, Phi, nm, ndof } = mdl;
  let maxOff = 0, maxDiagErr = 0;
  for (let a = 0; a < nm; a++) {
    for (let b = 0; b < nm; b++) {
      let s = 0;
      for (let i = 0; i < ndof; i++) {
        let t = 0;
        for (let j = 0; j < ndof; j++) t += M[i * ndof + j] * Phi[b * ndof + j];
        s += Phi[a * ndof + i] * t;
      }
      if (a === b) maxDiagErr = Math.max(maxDiagErr, Math.abs(s - 1));
      else maxOff = Math.max(maxOff, Math.abs(s));
    }
  }
  check('max |diag - 1|', maxDiagErr, 0, 1e-9);
  check('max |off-diagonal|', maxOff, 0, 1e-9);
}

console.log('\n=== 5. Nigam-Jennings vs analytical SDOF, ramp load p = c*t ===');
console.log('    A ramp IS piecewise linear, so the recurrence must be EXACT');
console.log('    at any dt.  This is the central accuracy claim.\n');
{
  const w = 2 * Math.PI * 5.0, z = 0.015, c = 1000;
  const wd = w * Math.sqrt(1 - z * z);
  const A = 2 * z * c / (w ** 3);
  const Bc = (z * w * A - c / (w * w)) / wd;
  for (const spc of [20, 5, 2.5]) {                  // samples per cycle
    const dt = 1 / (spc * 5.0);
    const nt = Math.ceil(3 / dt);
    const p = new Float64Array(nt);
    for (let i = 0; i < nt; i++) p[i] = c * i * dt;   // p(0) = 0
    const co = B.njCoeff(w, z, dt);
    const q = B.filter2(co.bq, co.a, p);
    const qd = B.filter2(co.bv, co.a, p);
    let eq = 0, ev = 0, rq = 0, rv = 0;
    for (let i = 0; i < nt; i++) {
      const t = i * dt, e = Math.exp(-z * w * t);
      const qa = (c / (w * w)) * (t - 2 * z / w) + e * (A * Math.cos(wd * t) + Bc * Math.sin(wd * t));
      const va = c / (w * w) + e * (-z * w * (A * Math.cos(wd * t) + Bc * Math.sin(wd * t))
                 + wd * (-A * Math.sin(wd * t) + Bc * Math.cos(wd * t)));
      eq = Math.max(eq, Math.abs(q[i] - qa)); rq = Math.max(rq, Math.abs(qa));
      ev = Math.max(ev, Math.abs(qd[i] - va)); rv = Math.max(rv, Math.abs(va));
    }
    check(`q  , ${spc} samples/cycle (rel err)`, eq / rq, 0, 1e-11);
    check(`qd , ${spc} samples/cycle (rel err)`, ev / rv, 0, 1e-11);
  }
}

console.log('\n=== 5b. Sine load: error is O(dt^2), from interpolating p only ===\n');
{
  const w = 2 * Math.PI * 5.0, z = 0.015, p0 = 1000, Om = 2 * Math.PI * 4.0;
  const wd = w * Math.sqrt(1 - z * z), r = Om / w;
  const D = (1 - r * r) ** 2 + (2 * z * r) ** 2;
  const Cs = (p0 / (w * w)) * (1 - r * r) / D, Cc = -(p0 / (w * w)) * (2 * z * r) / D;
  const A = -Cc;
  const Bc = (z * w * A - Om * Cs) / wd;
  let prev = 0;
  for (const spc of [10, 20, 40]) {
    const dt = 1 / (spc * 5.0);
    const nt = Math.ceil(2 / dt);
    const p = new Float64Array(nt);
    for (let i = 0; i < nt; i++) p[i] = p0 * Math.sin(Om * i * dt);
    const co = B.njCoeff(w, z, dt);
    const q = B.filter2(co.bq, co.a, p);
    let e2 = 0, rf = 0;
    for (let i = 0; i < nt; i++) {
      const t = i * dt, ex = Math.exp(-z * w * t);
      const qa = Cs * Math.sin(Om * t) + Cc * Math.cos(Om * t)
               + ex * (A * Math.cos(wd * t) + Bc * Math.sin(wd * t));
      e2 = Math.max(e2, Math.abs(q[i] - qa)); rf = Math.max(rf, Math.abs(qa));
    }
    const rel = e2 / rf;
    const ratio = prev ? prev / rel : NaN;
    console.log(`       ${String(spc).padStart(2)} samples/cycle: rel err ${rel.toExponential(2)}` +
                (prev ? `   (${ratio.toFixed(1)}x better, expect ~4x)` : ''));
    prev = rel;
  }
}

console.log('\n=== 6. Nigam-Jennings: static gain and velocity gain ===\n');
{
  const w = 2 * Math.PI * 7.3, z = 0.02, dt = 1 / 300;
  const { a, bq, bv } = B.njCoeff(w, z, dt);
  const sa = a[0] + a[1] + a[2];
  check('sum(bq)/sum(a) = 1/w^2', (bq[0] + bq[1] + bq[2]) / sa, 1 / (w * w), 1e-10);
  check('sum(bv)/sum(a) = 0', (bv[0] + bv[1] + bv[2]) / sa, 0, 1e-9);
  // undamped denominator must be [1, -2cos(w*dt), 1]
  const u = B.njCoeff(w, 0, dt);
  check('undamped a(2) = -2*cos(w*dt)', u.a[1], -2 * Math.cos(w * dt), 1e-12);
  check('undamped a(3) = 1', u.a[2], 1, 1e-12);
}

console.log('\n=== 7. Consistent load projection = direct phi evaluation ===');
console.log('    phi\'*(P*N(xi)\') must equal P*phi(x)\n');
{
  const inp = { L: [20, 20], m: 15e3, EI: 30e9, xi: 1.5, dL: 1, Nmod: 3,
                fmax: Infinity, bc: [[0, 1], [20, 1], [40, 1]] };
  const mdl = B.buildModel(inp);
  const { Phi, nm, ndof, dL, EL } = mdl;
  const P = 225e3;
  let maxErr = 0, ref = 0;
  for (const x of [0.37, 5.5, 13.91, 20.0, 27.3, 39.6]) {
    const e = Math.min(Math.max(Math.floor(x / dL), 0), EL - 1);
    const N = B.shapeFun(x / dL - e, dL);
    const direct = B.hermiteEval(Phi, nm, ndof, dL, EL, x);
    for (let n = 0; n < nm; n++) {
      // classical: build the nodal load vector, then project
      let proj = 0;
      for (let k = 0; k < 4; k++) proj += Phi[n * ndof + 2 * e + k] * (P * N[k]);
      maxErr = Math.max(maxErr, Math.abs(proj - P * direct[n]));
      ref = Math.max(ref, Math.abs(P * direct[n]));
    }
  }
  check('max |proj - P*phi(x)| (rel)', maxErr / ref, 0, 1e-12);
}

console.log('\n=== 8. Static deflection: quasi-static limit of the solver ===');
console.log('    one axle crossing very slowly, SS span: max w = -P*L^3/(48*EI)\n');
{
  const L = 20, m = 15e3, EI = 30e9, P = 1e5;
  const inp = { L: [L], m, EI, xi: 5, dL: 0.5, Nmod: 12, fmax: Infinity,
                bc: [[0, 1], [L, 1]], load: 'TEST', v: [1], W: 0 };
  // a one-axle "train" via the internal path
  const mdl = B.buildModel(inp);
  const { Phi, nm, ndof, dL, EL } = mdl;
  // static modal solution: q_n = p_n / w_n^2 with the load at midspan
  const ph = B.hermiteEval(Phi, nm, ndof, dL, EL, L / 2);
  let wmid = 0;
  for (let n = 0; n < nm; n++) {
    const pn = -P * ph[n];
    const wn = 2 * Math.PI * mdl.f[n];
    wmid += ph[n] * (pn / (wn * wn));
  }
  check('midspan deflection (12 modes)', wmid, -P * L ** 3 / (48 * EI), 2e-3);
}

console.log('\n=== 9. End-to-end run: D2 over a 2x20 m bridge ===\n');
{
  const inp = {
    L: [20, 20], m: 15e3, EI: 30e9, xi: 1.5, dL: 1, Nmod: 3, fmax: 30,
    bc: [[0, 1], [20, 1], [40, 1]],
    load: 'D2', N_car: 20, v: [], W: 0
  };
  for (let v = 100; v <= 200; v += 1) inp.v.push(v);
  const t0 = Date.now();
  const out = B.run(inp);
  const ms = Date.now() - t0;

  const env = out.amax[0];
  let am = 0, ai = 0;
  for (let i = 0; i < env.length; i++) if (env[i] > am) { am = env[i]; ai = i; }
  const de = out.dmax[0];
  let dm = -Infinity, di = 0;
  for (let i = 0; i < de.length; i++) if (de[i] > dm) { dm = de[i]; di = i; }

  console.log(`    f       = ${out.f.map(v => v.toFixed(3)).join(', ')} Hz`);
  console.log(`    fs      = ${out.fs} Hz`);
  console.log(`    a_max   = ${am.toFixed(4)} m/s^2 at ${out.v[ai]} km/h`);
  console.log(`    d_max   = ${(dm * 1e3).toFixed(4)} mm at ${out.v[di]} km/h`);
  console.log(`    runtime = ${ms} ms  (${inp.v.length} speeds)`);
  if (!isFinite(am) || am <= 0) { console.log(' FAIL  a_max not positive finite'); fails++; }
}

console.log('\n=== 10. Nmod = 1 (the shape bug that bit the MATLAB version) ===\n');
{
  const inp = {
    L: [20, 20], m: 15e3, EI: 30e9, xi: 1.5, dL: 1, Nmod: 1, fmax: 30,
    bc: [[0, 1], [20, 1], [40, 1]], load: 'D2', N_car: 20,
    v: [140, 150, 160], W: 0
  };
  const out = B.run(inp);
  const ok = out.amax[0].every(v => isFinite(v) && v > 0);
  console.log(`  ${ok ? '  ok  ' : ' FAIL '} single-mode run returns non-zero: ` +
              `[${Array.from(out.amax[0], v => v.toFixed(4)).join(', ')}]`);
  if (!ok) fails++;
}

console.log('\n' + (fails === 0 ? 'ALL CHECKS PASSED' : fails + ' CHECK(S) FAILED') + '\n');
process.exit(fails ? 1 : 0);

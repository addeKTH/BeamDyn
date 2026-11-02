/* solver.worker.js ----------------------------------------------------------
   One train (or a group of speeds of one train) per job; replaces the
   parfeval tasks of BeamDyn2D.m.  The model for the dynamics is built once
   per worker and kept while the input is unchanged.

   job: {id, key, spec (BD.dynModel input), train {x, P}, iv (speed indices)}
   reply: {id, ok, T (peaks at iv, static peaks)} or {id, ok: false, message}

   Andreas Andersson, 2026-11-01 (assisted by Claude)
---------------------------------------------------------------------------*/
importScripts('linalg.js?v=2026-11-02', 'beamdyn2d.js?v=2026-11-02');

let key = null, P = null;

self.onmessage = function (ev) {
    const job = ev.data;
    try {
        if (job.key !== key || !P) { P = BD.dynModel(job.spec); key = job.key; }
        const T = BD.speedPeaks(BD.trainPrep(job.train, P), P, job.iv);
        const buf = [];
        for (const k of BD.FIELDS) if (T[k]) buf.push(T[k].buffer);
        self.postMessage({ id: job.id, ok: true, T }, buf);
    } catch (e) {
        self.postMessage({ id: job.id, ok: false, message: e && e.message ? e.message : String(e) });
    }
};

/* solver.worker.js ----------------------------------------------------------
   One solver job at a time.  This is what replaces parfor: the GUI splits the
   run into jobs - one per train, or a block of speeds when there is only one
   train - and hands them to a pool of these workers.

   The bridge model does not change within a run, so it is built once per
   worker and cached; only the train and the speed list vary.  ds and the
   record length are always taken from the FULL speed list (vAll), so a run
   split over eight workers gives exactly the same numbers as a serial one.
---------------------------------------------------------------------------*/
importScripts('linalg.js', 'beammodel.js', 'trainloads.js',
              'beamdyn2d.js', 'beamstat2d.js', 'dynamicfactor.js');

let cacheKey = null, model = null, IL = null;

function modelKey(inp) {
    return JSON.stringify([inp.L, inp.m, inp.EI, inp.dL, inp.bc, inp.Nmod,
                           inp.fmax, inp.fs, inp.t_fr, inp.xi, inp.W, inp.secForces]);
}

function ensureModel(inp) {
    const key = modelKey(inp);
    if (key !== cacheKey) { cacheKey = key; model = null; IL = null; }
    if (!model) model = BD.buildModel(inp);
    return model;
}

self.onmessage = function (ev) {
    const job = ev.data;
    try {
        const inp = job.inp;
        inp.L = Array.from(inp.L);
        if (job.kind === 'dyn') {
            const m = ensureModel(inp);
            const train = { x: Float64Array.from(job.train.x), P: Float64Array.from(job.train.P) };
            const r = BD.runTrain(m, inp, train, job.vList, job.vAll);
            const out = { id: job.id, kind: 'dyn', ok: true, res: r, f: m.f };
            const buf = [];
            for (const k in r) if (r[k] && r[k].buffer) buf.push(r[k].buffer);
            self.postMessage(out, buf);
        } else if (job.kind === 'stat') {
            ensureModel(inp);
            if (!IL) IL = BS.influence(inp);
            const res = [];
            for (const t of job.trains) {
                res.push(BS.runTrain(IL, inp, { x: Float64Array.from(t.x), P: Float64Array.from(t.P) }));
            }
            self.postMessage({ id: job.id, kind: 'stat', ok: true, res: res });
        }
    } catch (e) {
        self.postMessage({ id: job.id, ok: false, message: e && e.message ? e.message : String(e) });
    }
};

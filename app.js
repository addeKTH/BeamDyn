/* app.js --------------------------------------------------------------------
   Railway bridge dynamics, web version of Main_gui.m (2026-11-01): input,
   model and eigenvalue analysis on every change, Run (Web Workers, Stop after
   5 s), result plots with tick boxes, HSLM-A vs A1-E5 envelope ('All'),
   phi'' (EN 1991-2:2024 Annex C), Open/Save of the input (JSON), Export
   (xlsx), Print (report) and Print of the result plots (browser print to PDF).

   Andreas Andersson, 2026-11-01 (assisted by Claude)
---------------------------------------------------------------------------*/
(function () {
'use strict';
const $ = id => document.getElementById(id);

// ---------------------------------------------------------- default input
const def = {
    L: '20 20', m: '16000', EI: '35E9', xi: '1.5', xiType: 'Concrete',
    Nmod: '3', fmax: '30', dL: '1', load: 'HSLM-A', N_car: '20', v: '100:2:300',
    bQType: '1-point', bQLast: '1-point', b_Q: '3.0'
};
const catList = ['A1', 'B1', 'B2', 'C2', 'C3', 'C4', 'D2', 'D3', 'D4', 'D5', 'E4', 'E5'];
const loadList = ['HSLM-A', ...catList, 'All'];
const xiTypes = ['Steel', 'Concrete', 'Pre-stressed', 'Other'];
const bQTypes = ['1-point', '3-point', '5-point', 'Triangular'];

// shared state (as Main_gui)
let inp = {}, mdl = null, runs = [], envData = null, sup = [], isEnv = false;
const ui = { autoN: true, autoF: true, xiType: def.xiType, xiLast: def.xiType,
             bQType: def.bQType, bQLast: def.bQLast, bQVal: def.b_Q, cumOn: false,
             busy: false, stopReq: false, expOn: false };
const prog = { done: 0, total: 1 };

// ============================================================ helpers
// numbers as MATLAB str2num would read them: '20 20', '[10 15 10]', '15E3',
// '100:1:200', simple expressions (+ - * / ^ and parentheses)
function evalExpr(s) {
    let i = 0;
    const peek = () => s[i], ws = () => { while (s[i] === ' ') i++; };
    function num() {
        ws();
        const m = /^(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?/.exec(s.slice(i));
        if (m) { i += m[0].length; return parseFloat(m[0]); }
        if (peek() === '(') { i++; const v = add(); ws(); if (s[i++] !== ')') throw 0; return v; }
        if (peek() === '-') { i++; return -pow(); }
        if (peek() === '+') { i++; return pow(); }
        if (s.slice(i, i + 2) === 'pi') { i += 2; return Math.PI; }
        throw 0;
    }
    function pow() { let v = num(); ws(); while (peek() === '^') { i++; v = Math.pow(v, num()); ws(); } return v; }
    function mul() { let v = pow(); ws(); while (peek() === '*' || peek() === '/') { const o = s[i++]; const w = pow(); v = o === '*' ? v * w : v / w; ws(); } return v; }
    function add() { let v = mul(); ws(); while (peek() === '+' || peek() === '-') { const o = s[i++]; const w = mul(); v = o === '+' ? v + w : v - w; ws(); } return v; }
    const v = add(); ws();
    if (i !== s.length || !isFinite(v)) throw 0;
    return v;
}
function parseVec(str) {
    try {
        let s = String(str).trim().replace(/^\[|\]$/g, '').replace(/,/g, ' ').trim();
        if (!s) return null;
        s = s.replace(/\s*:\s*/g, ':').replace(/([*/^+(-])\s+/g, '$1').replace(/\s+([*/^)])/g, '$1');
        const out = [];
        for (const tok of s.split(/\s+/)) {
            const p = tok.split(':').map(evalExpr);
            if (p.length === 1) out.push(p[0]);
            else {
                const [a, st, b] = p.length === 2 ? [p[0], 1, p[1]] : p;
                if (!st) return null;
                const n = Math.floor((b - a) / st + 1e-10) + 1;
                for (let k = 0; k < n; k++) out.push(Number((a + k * st).toPrecision(12)));
            }
        }
        return out;
    } catch (e) { return null; }
}
const scalar = s => { const v = parseVec(s); return v && v.length === 1 ? v[0] : null; };
const fmtFreq = f => (Math.round(f * 10) / 10 >= 100 ? f.toFixed(0) : f.toFixed(1));
const gfmt = x => String(Number(x.toPrecision(6)));           // %g
const maxrow = X => (X && X.length ? Float64Array.from(X[0], (_, j) => Math.max(...X.map(r => r[j]))) : null);
const minrow = X => (X && X.length ? Float64Array.from(X[0], (_, j) => Math.min(...X.map(r => r[j]))) : null);
function applyCum(y, v, on) {
    if (!on) return Array.from(y);
    const is = Array.from(v, (_, i) => i).sort((a, b) => v[a] - v[b]), out = Array.from(y);
    let m = -Infinity;
    for (const i of is) { if (y[i] > m || isNaN(m)) m = y[i]; out[i] = m; }
    return out;
}

// ============================================================ symbols
const C = { green: '#21963f', red: '#cc2626', blue: '#2659a6', yellow: '#f5c242', grey: '#737373' };
const ICON = {
    run: `<circle cx="12" cy="12" r="11" fill="${C.green}"/><polygon points="9.5,7 9.5,17 17.5,12" fill="#fff"/>`,
    stop: `<circle cx="12" cy="12" r="11" fill="${C.red}"/><rect x="8" y="8" width="8" height="8" rx="1" fill="#fff"/>`,
    reset: `<path d="M18.06 8.5 A7 7 0 1 1 12 5" fill="none" stroke="${C.red}" stroke-width="3"/><polygon points="12,0.8 12,9.2 17.5,5" fill="${C.red}"/>`,
    save: `<polygon points="3,3 18,3 21,6 21,21 3,21" fill="${C.blue}"/><rect x="6" y="3" width="10" height="6" fill="#d9e0ed"/>` +
          `<rect x="12.5" y="4.2" width="2" height="3.6" fill="${C.blue}"/><rect x="6" y="13" width="12" height="8" fill="#fff"/>`,
    load: `<polygon points="2,5 9,5 11,7 21,7 21,19 2,19" fill="#c49b35"/><polygon points="5,10 23,10 20,20 2,20" fill="${C.yellow}"/>`,
    print: `<rect x="7.5" y="3" width="9" height="6.5" fill="#fff" stroke="${C.grey}"/><rect x="2" y="8" width="20" height="10" rx="2" fill="#666b73"/>` +
           `<rect x="6" y="14" width="12" height="1.5" fill="#333338"/><rect x="7.5" y="15" width="9" height="6" fill="#fff" stroke="${C.grey}"/>` +
           `<path d="M9 17.5H15M9 19.5H15" stroke="${C.grey}"/><circle cx="18.5" cy="11" r="1" fill="#4dbf59"/>`,
    export: `<rect x="7" y="2" width="15" height="20" rx="1.5" fill="#219e59"/><path d="M14 7H22M14 12H22M14 17H22M18 2V22" stroke="#d9f2e0" stroke-width="1"/>` +
            `<rect x="2" y="6" width="13" height="13" rx="1.5" fill="#0f6636"/><path d="M5.56 8.3L11.44 16.7M11.44 8.3L5.56 16.7" stroke="#fff" stroke-width="2.1"/>`,
    clear: `<path d="M3.5 2V18M3 17.5H20" stroke="#595959" stroke-width="1.4"/>` +
           `<path d="M4.5 11.85 C8 11 9 6 11 6 C13 6 14 11 19.5 11.95" fill="none" stroke="#0072bd" stroke-width="2"/>` +
           `<circle cx="17" cy="17" r="6.5" fill="${C.red}"/><path d="M14.3 14.3L19.7 19.7M19.7 14.3L14.3 19.7" stroke="#fff" stroke-width="1.6"/>`,
    cum: `<path d="M3.5 2V21M3 20.5H22" stroke="#595959" stroke-width="1.4"/><polyline points="4.5,17 9,17 9,12 13,12 13,9 17,9 17,5 21.5,5" fill="none" stroke="#0072bd" stroke-width="2"/>`
};
ICON.cumOn = `<rect x="0" y="0" width="24" height="24" rx="4" fill="#cce0f7"/>` + ICON.cum;
const svgIcon = name => `<svg viewBox="0 0 24 24">${ICON[name]}</svg>`;

const SYMS = {
    sRun: ['run', 'Run: train analysis', () => (ui.busy ? onStop() : onRun())],
    sRst: ['reset', 'Reset: default input, clears the results', () => onReset()],
    sOpen: ['load', 'Open: input from a JSON file', () => $('fileIn').click()],
    sSave: ['save', 'Save: input to a JSON file', () => onSaveInput()],
    sExp: ['export', 'Export: results to *.xlsx', () => onExport()],
    sRep: ['print', 'Print: pdf report of the results', () => onReport()],
    sCum: ['cum', 'Cumulative: running peak over speed (on/off)', () => { ui.cumOn = !ui.cumOn; setIcon('sCum', ui.cumOn ? 'cumOn' : 'cum'); redrawResults(); }],
    sClr: ['clear', 'Clear: remove the results from the plots', () => onClear()],
    sPrp: ['print', 'Print: the result plots on one A4 page (pdf)', () => onPrintPlots()]
};
function setIcon(id, name) { $(id).innerHTML = svgIcon(name); }
function setEnable(ids, on) { for (const id of ids) $(id).classList.toggle('off', !on); }
const enabled = id => !$(id).classList.contains('off');
for (const id in SYMS) {
    const [ic, tip, fcn] = SYMS[id], el = $(id);
    setIcon(id, ic);
    el.dataset.tip = tip;
    el.addEventListener('click', () => {
        if (!enabled(id) || (ui.busy && id !== 'sRun')) return;
        hideTip(); fcn();
    });
    el.addEventListener('mouseenter', () => showTip(el));
    el.addEventListener('mouseleave', hideTip);
}
function showTip(el) {
    const t = $('tip'), r = el.offsetTop;
    t.textContent = el.dataset.tip; t.style.display = 'block';
    const w = t.offsetWidth;
    t.style.left = Math.min(el.offsetLeft, 550 - w) + 'px';
    t.style.top = (r + el.offsetHeight + 4) + 'px';
    el.style.cursor = enabled(el.id) ? 'pointer' : 'default';
}
function hideTip() { $('tip').style.display = 'none'; }

function setMsg(s, col) { const m = $('msg'); m.textContent = s || ''; m.style.color = col || '#4d4d4d'; }
function setProg(fr) {
    if (fr === null) { $('pfill').style.width = '0'; $('pct').textContent = ''; return; }
    $('pfill').style.width = (100 * fr) + '%'; $('pct').textContent = Math.round(100 * fr) + ' %';
}

// ============================================================ input fields
for (const s of loadList) { const o = document.createElement('option'); o.textContent = s; $('iLoad').appendChild(o); }
const fields = { iL: 'L', im: 'm', iEI: 'EI', iNmod: 'Nmod', ifmax: 'fmax', idL: 'dL', iNcar: 'N_car', iv: 'v' };
function setDefaults() {
    for (const id in fields) $(id).value = def[fields[id]];
    $('ixi').value = def.xi; ui.xiType = def.xiType; ui.xiLast = def.xiType;
    ui.bQType = def.bQType; ui.bQLast = def.bQLast; ui.bQVal = def.b_Q; $('ibQ').value = def.b_Q;
    $('iLoad').value = def.load;
    ui.autoN = true; ui.autoF = true;
}
for (const id of ['iL', 'im', 'iEI', 'idL', 'iNcar', 'iv']) $(id).addEventListener('change', () => updateModel());
$('iNmod').addEventListener('change', e => { ui.autoN = !e.target.value.trim(); updateModel(); });
$('ifmax').addEventListener('change', e => { ui.autoF = !e.target.value.trim(); updateModel(); });
$('iLoad').addEventListener('change', () => {
    updateModel();
    if (envMode() && !envData) setMsg('Press Run to compute the envelope');
    redrawResults();
});
$('iModes').addEventListener('change', () => plotBeam());

// damping: one field with a menu button; Table 8.6 (grey) or a typed value
// ('Other', black); an empty field returns to the last Table 8.6 type
$('ixi').addEventListener('change', e => {
    ui.xiType = e.target.value.trim() ? 'Other' : ui.xiLast;
    updateModel();
});
$('bxi').addEventListener('click', ev => {
    ev.stopPropagation();
    const Lv = parseVec($('iL').value), items = [];
    xiTypes.slice(0, 3).forEach(t => items.push({
        text: t + (Lv && Lv.every(x => x > 0) ? '   ' + BD.EN.xi(t, BD.EN.LPhi(Lv)).toFixed(2) + ' %' : ''),
        chk: ui.xiType === t, fn: () => { ui.xiType = t; ui.xiLast = t; updateModel(); } }));
    items.push('-', { text: 'Other (type a value)', chk: ui.xiType === 'Other', fn: () => { ui.xiType = 'Other'; updateModel(); } });
    openMenu('mXi', $('bxi'), items);
});
// load distribution: 1-, 3-, 5-point by name (grey) or a typed b_Q
// (Triangular, black); an empty field returns to the last point type
$('ibQ').addEventListener('change', e => {
    const s = e.target.value.trim();
    if (!s || bQTypes.slice(0, 3).includes(s)) ui.bQType = ui.bQLast;
    else { ui.bQType = 'Triangular'; ui.bQVal = s; }
    updateModel();
});
$('bbQ').addEventListener('click', ev => {
    ev.stopPropagation();
    const names = ['1-point   (axle load)', '3-point   25/50/25 % at 0.6 m', '5-point   Winkler beam, 100 MN/m'];
    const items = names.map((n, k) => ({ text: n, chk: ui.bQType === bQTypes[k],
        fn: () => { ui.bQType = bQTypes[k]; ui.bQLast = bQTypes[k]; updateModel(); } }));
    items.push('-', { text: 'Triangular (type b_Q)', chk: ui.bQType === 'Triangular',
        fn: () => { ui.bQType = 'Triangular'; $('ibQ').value = ui.bQVal; updateModel(); } });
    openMenu('mBQ', $('bbQ'), items);
});
function openMenu(id, btn, items) {
    const m = $(id);
    m.innerHTML = '';
    for (const it of items) {
        if (it === '-') { m.appendChild(document.createElement('hr')); continue; }
        const d = document.createElement('div');
        d.textContent = it.text; if (it.chk) d.className = 'chk';
        d.addEventListener('click', () => { m.style.display = 'none'; it.fn(); });
        m.appendChild(d);
    }
    m.style.left = btn.offsetLeft + 'px'; m.style.top = (btn.offsetTop + btn.offsetHeight) + 'px';
    m.style.display = 'block';
}
document.addEventListener('click', () => { $('mXi').style.display = 'none'; $('mBQ').style.display = 'none'; removeSupSel(); });

// tick boxes
const tick = id => $(id).checked;
function excl(a, b) { $(a).addEventListener('change', () => { $(a).checked = true; $(b).checked = false; redrawResults(); }); }
excl('tD', 'tR'); excl('tR', 'tD'); excl('tMx', 'tMn'); excl('tMn', 'tMx');
$('tPhi').addEventListener('change', () => redrawResults());
$('tPbis').addEventListener('change', () => { pbisLabel(); redrawResults(); });
$('tSep').addEventListener('change', () => { if (tick('tSep')) setMsg('separate results for HSLM A1-A10 trains'); redrawResults(); });

// ============================================================ input / model
function readInput() {
    const L = parseVec($('iL').value), m = scalar($('im').value), EI = scalar($('iEI').value);
    let xi;
    if (ui.xiType === 'Other') { xi = scalar($('ixi').value); $('ixi').classList.remove('auto'); }
    else { xi = L && L.every(x => x > 0) ? BD.EN.xi(ui.xiType, BD.EN.LPhi(L)) : null; $('ixi').classList.add('auto'); }
    $('ixi').title = `Modal damping (%), ${ui.xiType.toLowerCase()}: type a value, or pick EN 1991-2 Table 8.6 ` +
                     '(L = L_Phi) with the button; empty = back to Table 8.6';
    if (xi !== null && isFinite(xi)) {             // two decimals, one from 10 %
        const s = Math.abs(Math.round(xi * 100) / 100) >= 10 ? xi.toFixed(1) : xi.toFixed(2);
        xi = parseFloat(s); $('ixi').value = s;
    }
    // load distribution: point types by name (grey), triangular by b_Q (black)
    let b_Q = 0;
    if (ui.bQType === 'Triangular') {
        let s = $('ibQ').value.trim();
        if (bQTypes.includes(s)) s = ui.bQVal;
        $('ibQ').value = s; $('ibQ').classList.remove('auto');
        b_Q = scalar(s); ui.bQVal = s;
    } else { $('ibQ').value = ui.bQType; $('ibQ').classList.add('auto'); }
    $('ibQ').title = `Load distribution of each axle, ${ui.bQType.toLowerCase()}: pick 1-, 3- or 5-point with the ` +
                     'button, or type the base b_Q (m) of a triangle; empty = back to the point distribution';
    const r = {
        L, m, EI, xi, xiType: ui.xiType,
        Nmod: ui.autoN ? 60 : scalar($('iNmod').value), fmax: ui.autoF ? Infinity : scalar($('ifmax').value),
        dL: scalar($('idL').value), load: $('iLoad').value, N_car: scalar($('iNcar').value),
        v: parseVec($('iv').value), loadDist: ui.bQType, b_Q
    };
    const ok = L && L.length && L.every(x => x > 0) && m > 0 && EI > 0 && xi !== null && xi >= 0 && xi < 100 &&
        r.Nmod >= 1 && r.fmax > 0 && r.dL > 0 && r.v && r.v.length && r.v.every(x => x > 0) &&
        b_Q !== null && b_Q >= 0 && r.N_car >= 1 && r.N_car % 1 === 0;
    if (!ok) { setMsg('input error', '#b30000'); return null; }
    $('iNcar').disabled = r.load === 'HSLM-A';
    const hs = r.load === 'HSLM-A' ? '' : 'none';
    $('tSep').style.display = hs; $('lSep').style.display = hs;
    return r;
}

function updateModel() {
    const r = readInput();
    if (!r) return false;
    inp = r;
    setMsg('');
    let model, eg;
    try {
        model = BD.buildModel({ L: inp.L, m: inp.m, EI: inp.EI, dL: inp.dL,
                                clamped: sup.length === inp.L.length + 1 ? sup : [] });
        if (sup.length !== model.xs.length) sup = model.xs.map(() => false);
        const nf = model.fdof.length;
        const fmaxOf = f => ui.autoF ? Math.ceil(Math.max(30, 1.5 * f[0], f[Math.min(2, f.length - 1)]) - 1e-9) : inp.fmax;
        eg = BD.eigen(model, { nreq: Math.min(Math.max(inp.Nmod, 3), nf), fmaxFn: fmaxOf });
        // f_max: EN 1990:2023 A.2.9.4.2.1(3) unless typed in; keep the modes up
        // to f_max and, if typed in, not more than N_mod; never fewer than one
        if (ui.autoF) { inp.fmax = fmaxOf(eg.f); $('ifmax').value = gfmt(inp.fmax); }
        let keep = Array.from(eg.f, (f, i) => f <= inp.fmax && (ui.autoN || i < inp.Nmod));
        if (!keep.some(Boolean)) keep[0] = true;
        const idx = keep.map((k, i) => (k ? i : -1)).filter(i => i >= 0);
        const nd = model.ndof, Phi = new Float64Array(nd * idx.length);
        idx.forEach((i, c) => Phi.set(eg.Phi.subarray(i * nd, (i + 1) * nd), c * nd));
        inp.Nmod = idx.length;
        if (ui.autoN) $('iNmod').value = String(inp.Nmod);
        $('iNmod').classList.toggle('auto', ui.autoN); $('ifmax').classList.toggle('auto', ui.autoF);
        const f = Float64Array.from(idx, i => eg.f[i]), fe = f[f.length - 1];
        inp.fs = Math.ceil(20 * fe);
        inp.t_fr = Math.round(Math.max(2 / fe, 1 / f[0]) * 100) / 100;
        mdl = Object.assign(model, { f, Phi });
    } catch (e) {
        setMsg('input error: ' + e.message, '#b30000'); return false;
    }
    setEnable(['sExp', 'sRep', 'sPrp'], false);
    plotBeam(); applyVlim(); pbisLabel();
    window.inp = inp; window.mdl = mdl;
    return true;
}

// ============================================================ model plot
// beam, supports, mode shapes (amplitude 0.08 L, four at a time), dimension
// lines, 1:1 aspect; legend to the right (GUI) or on top (report)
function drawModel(svg, box, opt) {
    opt = opt || {};
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    const Lt = mdl.sumL, a = 0.08 * Lt, yd = -(a + 0.045 * Lt);
    const xl = [-0.02 * Lt, 1.02 * Lt];
    const yr = [yd - 0.035 * Lt, a + (opt.top ? 0.075 : 0.03) * Lt];
    const ysp = Math.max((xl[1] - xl[0]) * box.h / box.w, yr[1] - yr[0]);
    const yl = [(yr[0] + yr[1]) / 2 - ysp / 2, (yr[0] + yr[1]) / 2 + ysp / 2];
    const X = x => box.x + (x - xl[0]) / (xl[1] - xl[0]) * box.w;
    const Y = y => box.y + box.h - (y - yl[0]) / (yl[1] - yl[0]) * box.h;
    const sz = opt.size || 14.7;
    PL.el('line', { x1: X(0), y1: Y(0), x2: X(Lt), y2: Y(0), stroke: '#8c8c8c', 'stroke-width': 4 }, svg);
    // supports
    mdl.xs.forEach((x, i) => {
        if (mdl.clamped[i]) PL.el('line', { x1: X(x), y1: Y(0.03 * Lt), x2: X(x), y2: Y(-0.03 * Lt), stroke: '#000', 'stroke-width': 1 }, svg);
        else {
            PL.el('polygon', { points: `${X(x - 0.02 * Lt)},${Y(-0.03 * Lt)} ${X(x + 0.02 * Lt)},${Y(-0.03 * Lt)} ${X(x)},${Y(0)}`,
                               fill: '#fff', stroke: '#000', 'stroke-width': 0.7 }, svg);
            PL.el('circle', { cx: X(x), cy: Y(0), r: 1.8, fill: '#000' }, svg);
        }
    });
    // mode shapes of the chosen group, on top of the supports
    const nf = mdl.f.length, g = opt.group !== undefined ? opt.group : (+$('iModes').value || 0);
    const im = []; for (let n = 4 * g; n < Math.min(4 * g + 4, nf); n++) im.push(n);
    const xp = Float64Array.from({ length: 1001 }, (_, i) => i * Lt / 1000);
    const lg = [];
    im.forEach((n, j) => {
        const ph = BD.modeShape(mdl, mdl.Phi, n, xp);
        let mx = 0; for (const v of ph) mx = Math.max(mx, Math.abs(v));
        let d = '';
        for (let i = 0; i < xp.length; i++) d += (i ? 'L' : 'M') + X(xp[i]).toFixed(2) + ' ' + Y(a * ph[i] / mx).toFixed(2);
        PL.el('path', { d, fill: 'none', stroke: PL.color(j), 'stroke-width': 1 }, svg);
        lg.push({ c: PL.color(j), s: `$f_{${n + 1}}$ = ${opt.top ? mdl.f[n].toFixed(1) : fmtFreq(mdl.f[n])} Hz` });
    });
    // dimension lines below the lowest mode-shape point
    const pos = [...new Set([0, ...mdl.xs, Lt])].sort((p, q) => p - q);
    const ln = (x1, y1, x2, y2) => PL.el('line', { x1: X(x1), y1: Y(y1), x2: X(x2), y2: Y(y2), stroke: '#000', 'stroke-width': 0.5 }, svg);
    ln(0, yd, Lt, yd);
    for (const p of pos) { ln(p, yd - 0.016 * Lt, p, yd + 0.024 * Lt); ln(p - 0.012 * Lt, yd - 0.012 * Lt, p + 0.012 * Lt, yd + 0.012 * Lt); }
    for (let i = 0; i < pos.length - 1; i++)
        PL.text(svg, (pos[i + 1] - pos[i]).toFixed(1), X((pos[i] + pos[i + 1]) / 2), Y(yd + 0.004 * Lt), { size: sz, anchor: 'middle', baseline: 'text-after-edge' });
    // legend
    const ls = opt.legendSize || 13;
    if (opt.top) {
        let w = 0; const ws = lg.map(e => 15 + 5 + PL.textWidth(e.s, ls) + 14);
        w = ws.reduce((p, q) => p + q, 0);
        let x0 = box.x + (box.w - w) / 2;
        lg.forEach((e, k) => {
            PL.el('line', { x1: x0, y1: box.y + 8, x2: x0 + 15, y2: box.y + 8, stroke: e.c, 'stroke-width': 1 }, svg);
            PL.text(svg, e.s, x0 + 20, box.y + 8, { size: ls }); x0 += ws[k];
        });
    } else {
        const rowH = ls * 1.1, w = Math.max(...lg.map(e => 20 + PL.textWidth(e.s, ls)));
        const x0 = Math.min(box.x + box.w + 4, 548 - w - 2), y0 = box.y + (box.h - lg.length * rowH) / 2;
        lg.forEach((e, k) => {
            const y = y0 + (k + 0.5) * rowH;
            PL.el('line', { x1: x0, y1: y, x2: x0 + 15, y2: y, stroke: e.c, 'stroke-width': 1 }, svg);
            PL.text(svg, e.s, x0 + 20, y, { size: ls });
        });
    }
    return { X, Y };
}
let modelMap = null;
function plotBeam() {
    // Modes pull-down: groups of four, shown with more than four modes
    const nf = mdl.f.length, ng = Math.ceil(nf / 4), sel = $('iModes'), cur = Math.min(+sel.value || 0, ng - 1);
    sel.innerHTML = '';
    for (let g = 0; g < ng; g++) {
        const o = document.createElement('option'); o.value = g;
        o.textContent = 4 * g + 1 < nf ? `${4 * g + 1}-${Math.min(4 * g + 4, nf)}` : String(4 * g + 1);
        sel.appendChild(o);
    }
    sel.value = cur;
    const vis = nf > 4 ? '' : 'none';
    sel.style.display = vis; $('lModes').style.display = vis;
    modelMap = drawModel($('model'), { x: 48, y: 5, w: 410, h: 120 });
}
// click near a support: pull-down pinned / clamped
$('model').addEventListener('click', ev => {
    ev.stopPropagation();
    if (ui.busy || !mdl) return;
    const r = $('model').getBoundingClientRect(), px = ev.clientX - r.left;
    let bi = -1, bd = Infinity;
    mdl.xs.forEach((x, i) => { const d = Math.abs(modelMap.X(x) - px); if (d < bd) { bd = d; bi = i; } });
    if (bd > 0.03 * 410 + 4) return;
    removeSupSel();
    const s = document.createElement('select');
    s.id = 'supsel';
    for (const t of ['pinned', 'clamped']) { const o = document.createElement('option'); o.textContent = t; s.appendChild(o); }
    s.value = sup[bi] ? 'clamped' : 'pinned';
    s.style.left = (10 + modelMap.X(mdl.xs[bi]) - 32) + 'px';
    s.style.top = (220 + modelMap.Y(0) + 12) + 'px';
    s.addEventListener('click', e => e.stopPropagation());
    s.addEventListener('change', () => { sup[bi] = s.value === 'clamped'; removeSupSel(); updateModel(); });
    $('fig').appendChild(s);
    s.focus();
});
function removeSupSel() { const s = $('supsel'); if (s) s.remove(); }

// ============================================================ result plots
const resSvg = $('res');
const AX = [10, 126, 242, 358].map((y, k) => new PL.Axes(resSvg, { x: 48, y, w: 410, h: 110 }, { xTickLabels: k === 3 }));
AX[0].ylabel = '$a_\\mathrm{max}$ (m/s$^2$)';
AX[3].xlabel = '$v$ (km/h)';
AX.forEach((ax, k) => {
    ax.onTipText = (x, y, s) => {
        let vs;
        if (k === 0) vs = `${y.toFixed(1)} m/s$^2$`;
        else if (tick('tPhi')) vs = `$\\varphi'$ = ${y.toFixed(2)}`;
        else {                                      // unit from the y-label
            const u = (/\(([^)]*)\)\s*$/.exec(ax.ylabel) || ['', ''])[1];
            vs = u === 'MNm' ? `${(y * 1e3).toFixed(0)} kNm` : u === 'MN' ? `${(y * 1e3).toFixed(0)} kN` : `${y.toFixed(1)} ${u}`;   // mm, mrad
        }
        // one of A1-A10: the train on the first row
        return [(s && s.train ? s.train + '@' : '') + `${x.toFixed(0)} km/h`, vs];
    };
});
const envMode = () => inp.load === 'All';
const sepShown = () => tick('tSep') && inp.load === 'HSLM-A';
function inputVlim() {
    const v = inp.v;
    return v && v.length > 1 && Math.max(...v) > Math.min(...v) ? [Math.min(...v), Math.max(...v)] : null;
}
function applyVlim() {
    const vl = inputVlim();
    if (!vl) return;
    AX.forEach(ax => { ax.xlim = vl; });
    redrawResults();
}

// y-limits from the data range: from zero (padded on the data side), or
// fitted to the data (fit, padded on both sides)
function setLimFrom(ax, rg, fit) {
    if (!rg.every(isFinite)) { ax.ylim = [0, 1]; return; }
    if (fit) {
        let pad = Math.max(0.1 * (rg[1] - rg[0]), 0.02 * Math.max(Math.abs(rg[0]), Math.abs(rg[1])));
        if (pad === 0) pad = 1;
        ax.ylim = [rg[0] - pad, rg[1] + pad]; return;
    }
    const lo = Math.min(rg[0], 0); let hi = Math.max(rg[1], 0);
    if (hi <= lo) hi = lo + 1e-12;
    const pad = 0.1 * (hi - lo);
    ax.ylim = [lo - pad * (lo < 0), hi + pad * (hi > 0)];
}

// draw a_max, d_max / phi_max, M_max / M_min, V_max into the axes A (the GUI
// axes or the report's); with the phi'' tick on scaled results (withPbis)
function redrawResults(A) {
    A = A || AX;
    withPbis(() => drawResults(A));
    if (A === AX) AX.forEach(ax => ax.draw());
}
function drawResults(A) {
    const cumOn = ui.cumOn, sepOn = sepShown();
    A.forEach(ax => ax.clear());
    const ax3 = A[0];
    let vlim = null;
    if (envMode() && envData) {
        const v = envData.v; vlim = [Math.min(...v), Math.max(...v)];
        ax3.plot(v, applyCum(envData.ha, v, cumOn), { color: PL.color(0) });
        ax3.plot(v, applyCum(envData.oa, v, cumOn), { color: PL.color(1) });
        ax3.legend = { names: ['HSLM-A: ' + peakText('a', peakOf(envData.ha)), 'A1-E5: ' + peakText('a', peakOf(envData.oa))],
                       colors: [PL.color(0), PL.color(1)] };
    } else {
        let ci = 0;
        const names = [], cols = [];
        for (const rr of runs) {
            const v = rr.v;
            vlim = vlim ? [Math.min(vlim[0], ...v), Math.max(vlim[1], ...v)] : [Math.min(...v), Math.max(...v)];
            names.push(runText(rr, 'a', peakOf(rr.a)));
            if (sepOn && rr.a.length > 1) {         // the ten HSLM-A trains; the peak in the legend,
                rr.a.forEach((row, t) => ax3.plot(v, applyCum(row, v, cumOn), { color: PL.color(ci++), train: 'A' + (t + 1) }));   // the train in the tip
                cols.push(null);
            } else { const c = PL.color(ci++); ax3.plot(v, applyCum(maxrow(rr.a), v, cumOn), { color: c }); cols.push(c); }
        }
        if (runs.length) ax3.legend = { names, colors: cols, size: names.length > 4 ? 10.7 : undefined };
    }
    let ymax = -Infinity;
    for (const s of ax3.series) for (const y of s.y) if (isFinite(y)) ymax = Math.max(ymax, y);
    ax3.ylim = [0, ymax > -Infinity ? Math.max(1.1 * ymax, 1e-12) : 1];
    drawCfg(A[1], tick('tR') ? 4 : 0, false, sepOn, cumOn);
    drawCfg(A[2], tick('tMn') ? 2 : 1, true, sepOn, cumOn);
    drawCfg(A[3], 3, true, sepOn, cumOn);
    let vl = inputVlim() || vlim;
    if (vl) { if (vl[1] <= vl[0]) vl = [vl[0] - 1, vl[1] + 1]; A.forEach(ax => { ax.xlim = vl; }); }
}

// one result plot: q = 0 d_max, 1 M_max, 2 M_min, 3 V_max, 4 phi_max (rotation)
function drawCfg(ax, q, sfMode, sepOn, cumOn) {
    const phiOn = tick('tPhi');
    const qraw = ['d', 'M', 'Mn', 'V', 'r'], qphi = ['pd', 'pM', 'pMn', 'pV', 'pr'];
    const ylab = ['$d_\\mathrm{max}$ (mm)', '$M$ (MNm)', '$M$ (MNm)', '$V$ (MN)', '$\\phi_\\mathrm{max}$ (mrad)'];
    const sc = [1e3, 1e-6, 1e-6, 1e-6, 1e3];
    const k = phiOn ? 1 : sc[q], fld = phiOn ? qphi[q] : qraw[q];
    const key = phiOn ? 'p' : qraw[q], useMin = q === 2 && !phiOn;     // legend text
    let yr = [Infinity, -Infinity];
    const names = [], cols = [];
    const line = (v, y, c, train) => {
        const yy = applyCum(y, v, cumOn).map(x => x * k);
        ax.plot(v, yy, { color: c, train });
        for (const x of yy) if (isFinite(x)) { yr[0] = Math.min(yr[0], x); yr[1] = Math.max(yr[1], x); }
    };
    if (envMode() && envData) {
        line(envData.v, envData['h' + fld], PL.color(0));
        line(envData.v, envData['o' + fld], PL.color(1));
        names.push('HSLM-A: ' + peakText(key, peakOf(envData['h' + fld], useMin)), 'A1-E5: ' + peakText(key, peakOf(envData['o' + fld], useMin)));
        cols.push(PL.color(0), PL.color(1));
    } else {
        runs.forEach((rr, r) => {
            const X = rr[fld];
            if (!X) return;
            names.push(runText(rr, key, peakOf(X, useMin)));
            if (sepOn && X.length > 1) {            // the ten HSLM-A trains; the peak in the legend, the train in the tip
                X.forEach((row, t) => line(rr.v, row, PL.color(t), 'A' + (t + 1)));
                cols.push(null);
            } else { line(rr.v, useMin ? minrow(X) : maxrow(X), PL.color(r)); cols.push(PL.color(r)); }
        });
    }
    const yphi = ['\\mathrm{dmax}', '\\mathrm{Mmax}', '\\mathrm{Mmin}', '\\mathrm{Vmax}', '\\phi\\mathrm{max}'];
    ax.ylabel = phiOn ? `$\\varphi'_{${yphi[q]}}$ (-)` : ylab[q];
    setLimFrom(ax, yr, sfMode && !phiOn);
    if (ax.series.length) ax.legend = { names, colors: cols, size: names.length > 4 ? 10.7 : undefined };
}

// legend text of a peak: a m/s^2, d m, rotation rad, M Nm, V N, phi' -
function peakText(key, val) {
    switch (key) {
        case 'a':  return `$a_\\mathrm{max}$ = ${val.toFixed(1)} m/s$^2$`;
        case 'd':  return `$d_\\mathrm{max}$ = ${(val * 1e3).toFixed(1)} mm`;
        case 'r':  return `$\\phi_\\mathrm{max}$ = ${(val * 1e3).toFixed(1)} mrad`;
        case 'M':  return `$M_\\mathrm{max}$ = ${(val / 1e3).toFixed(0)} kNm`;
        case 'Mn': return `$M_\\mathrm{min}$ = ${(val / 1e3).toFixed(0)} kNm`;
        case 'V':  return `$V_\\mathrm{max}$ = ${(val / 1e3).toFixed(0)} kN`;
        default:   return `$\\varphi'$ = ${val.toFixed(2)}`;
    }
}
// the same for a run, named when there are several
function runText(rr, key, val) { return (runs.length === 1 ? '' : rr.load + ': ') + peakText(key, val); }
// peak of one row or of rows [ntr][nv], finite values only
function peakOf(X, useMin) {
    let b = NaN;
    for (const row of (X[0] && X[0].length !== undefined ? X : [X]))
        for (const y of row) if (isFinite(y) && !(useMin ? y >= b : y <= b)) b = y;
    return b;
}

// ============================================================ phi''
function ppFactor(v) { return Float64Array.from(v, x => 1 + 0.5 * BD.EN.phiBis(x, mdl.f[0], BD.EN.LPhi(mdl.L))); }
function pbisLabel() {
    let s = '';
    if (tick('tPbis') && mdl && inp.v)
        s = PL.texHTML(`$1 + \\varphi''/2$ = ${(1 + 0.5 * BD.EN.phiBis(Math.max(...inp.v), mdl.f[0], BD.EN.LPhi(mdl.L))).toFixed(2)}`);
    $('rfac').innerHTML = s;
}
// fcn on results including phi'' (EN 1991-2:2024 8.4.6.5(3)-(4)): load
// effects y_dyn + phi''/2*y_stat, with y_stat = y_dyn/(1 + phi'); the
// acceleration times 1 + phi''/2; phi' itself not scaled
let inPbis = false;
function withPbis(fcn) {
    if (!tick('tPbis') || inPbis) return fcn();     // already scaled (report -> plots)
    const r0 = runs, e0 = envData;
    [runs, envData] = scaledResults(runs, envData);
    inPbis = true;
    try { return fcn(); } finally { runs = r0; envData = e0; inPbis = false; }
}
function addPbis(Y, Pp, k) {
    if (!Y) return Y;
    return Y.map((row, t) => Float64Array.from(row, (y, j) => {
        let ys = y / (1 + Pp[t][j]); if (!isFinite(ys)) ys = 0;
        return y + (k[j] - 1) * ys;
    }));
}
function scaledResults(rs, e) {
    const fr = ['d', 'r', 'M', 'Mn', 'V'], fp = ['pd', 'pr', 'pM', 'pMn', 'pV'];
    rs = rs.map(r0 => {
        const r = Object.assign({}, r0), k = r.fpp;
        r.a = r.a.map(row => Float64Array.from(row, (y, j) => y * k[j]));
        fr.forEach((f, j) => { r[f] = addPbis(r[f], r[fp[j]], k); });
        return r;
    });
    if (!e) return [rs, e];
    const k = e.fpp, E = Object.assign({}, e);
    const hd = Object.assign({}, e.hdyn);
    hd.amax = hd.amax.map(row => Float64Array.from(row, (y, j) => y * k[j]));
    for (const [f, p] of [['dmax', 'dmax'], ['rmax', 'rmax'], ['Mmax', 'Mmax'], ['Mmin', 'Mmin'], ['Vmax', 'Vmax']]) hd[f] = addPbis(hd[f], e.hphi[p], k);
    E.hdyn = hd;
    E.cat = e.cat.map(c0 => {
        const c = Object.assign({}, c0);
        c.a = Float64Array.from(c.a, (y, j) => y * k[j]);
        fr.forEach((f, j) => { c[f] = addPbis([c[f]], [c[fp[j]]], k)[0]; });
        return c;
    });
    E.ha = maxrow(hd.amax); E.hd = maxrow(hd.dmax); E.hr = maxrow(hd.rmax);
    E.hM = maxrow(hd.Mmax); E.hMn = minrow(hd.Mmin); E.hV = maxrow(hd.Vmax);
    const cm = (f, mn) => Float64Array.from(E.cat[0][f], (_, j) => (mn ? Math.min : Math.max)(...E.cat.map(c => c[f][j])));
    E.oa = cm('a'); E.od = cm('d'); E.or = cm('r'); E.oM = cm('M'); E.oMn = cm('Mn', true); E.oV = cm('V');
    return [rs, E];
}

// ============================================================ the run
// worker pool: one job per train, or per group of speeds for one train
const nWorkers = Math.max(1, Math.min(8, (navigator.hardwareConcurrency || 2) - 1));
let pool = null, useWorkers = true;
function makePool() {
    try {
        pool = Array.from({ length: nWorkers }, () => new Worker('solver.worker.js?v=2026-11-02d'));
    } catch (e) { pool = null; useWorkers = false; }
}
function killPool() { if (pool) pool.forEach(w => w.terminate()); pool = null; }
function runJobs(jobs, onDone) {
    return new Promise((resolve, reject) => {
        if (!pool) makePool();
        if (!pool) { reject(new Error('noworkers')); return; }
        const res = new Array(jobs.length);
        let next = 0, done = 0, failed = false;
        ui.abort = () => { failed = true; killPool(); reject(new Error('stopped')); };
        const feed = (w) => {
            if (failed || next >= jobs.length) return;
            const j = next++;
            w.onmessage = ev => {
                if (failed) return;
                const r = ev.data;
                if (!r.ok) { failed = true; killPool(); reject(new Error(r.message)); return; }
                res[j] = r.T; done++; onDone(done, jobs.length);
                if (done === jobs.length) resolve(res); else feed(w);
            };
            w.onerror = e => { if (!failed) { failed = true; killPool(); reject(new Error(e.message || 'worker error')); } };
            w.postMessage(jobs[j]);
        };
        pool.forEach(feed);
    });
}
// serial fallback (no workers, e.g. a page opened from file://): chunks of
// speeds with a pause in between, so the page stays responsive and Stop works
async function runSerial(jobs, onDone) {
    const P = BD.dynModel(jobs[0].spec), res = [];
    for (let j = 0; j < jobs.length; j++) {
        const Q = BD.trainPrep(jobs[j].train, P), iv = jobs[j].iv;
        let T = null;
        for (let c = 0; c < iv.length; c += 4) {
            const Tc = BD.speedPeaks(Q, P, iv.slice(c, c + 4));
            if (!T) T = Tc; else for (const k of BD.FIELDS) if (T[k]) for (const i of iv.slice(c, c + 4)) T[k][i] = Tc[k][i];
            await new Promise(r => setTimeout(r, 0));
            if (ui.stopReq) throw new Error('stopped');
        }
        res.push(T); onDone(j + 1, jobs.length);
    }
    return res;
}

// dynamic and static response of one load model (dyn [ntr][nv], stat, phi)
async function runCase(loadName) {
    const trains = BD.trainLoads(loadName, loadName === 'HSLM-A' ? [] : inp.N_car);
    const spec = { L: mdl.L, dL: mdl.dL, EI: inp.EI, m: inp.m, clamped: mdl.clamped, f: Array.from(mdl.f), Phi: mdl.Phi,
                   xi: inp.xi, fs: inp.fs, t_fr: inp.t_fr, v: inp.v, loadDist: inp.loadDist, b_Q: inp.b_Q };
    const key = JSON.stringify(Object.assign({}, spec, { Phi: mdl.Phi.length }));
    const nv = inp.v.length, ntr = trains.length, all = inp.v.map((_, i) => i);
    let jobs;
    if (ntr > 1) jobs = trains.map((t, k) => ({ id: k, key, spec, train: t, iv: all }));
    else {
        const ng = useWorkers ? Math.min(nWorkers, nv) : 1;
        jobs = Array.from({ length: ng }, (_, g) => ({ id: g, key, spec, train: trains[0], iv: all.filter(i => i % ng === g) }));
    }
    const onDone = (k, n) => {
        if (ntr > 1) setMsg(`${loadName}: ${k}/${ntr}`); else setMsg(loadName);
        setProg((prog.done + (ntr > 1 ? k : k / n)) / prog.total);
    };
    let R;
    try { R = useWorkers ? await runJobs(jobs, onDone) : await runSerial(jobs, onDone); }
    catch (e) {
        if (e.message !== 'noworkers') throw e;
        useWorkers = false; R = await runSerial(jobs.length > 1 && ntr === 1 ? [Object.assign({}, jobs[0], { iv: all })] : jobs, onDone);
    }
    // collect [ntr][nv]; one train in groups of speeds: merge the groups
    const dyn = {}, stat = {};
    for (const k of BD.FIELDS) dyn[k] = [];
    for (const k of ['dmax', 'rmax', 'Mmax', 'Mmin', 'Vmax', 'RFmax', 'RFmin', 'RMmax']) stat[k] = [];
    if (ntr > 1 || R.length === 1) R.forEach(T => {
        for (const k of BD.FIELDS) dyn[k].push(T[k]);
        for (const k in stat) stat[k].push(T.stat[k]);
    });
    else {
        const T = R[0];
        R.slice(1).forEach((Tg, g) => { for (const k of BD.FIELDS) if (T[k]) for (const i of jobs[g + 1].iv) T[k][i] = Tg[k][i]; });
        for (const k of BD.FIELDS) dyn[k].push(T[k]);
        for (const k in stat) stat[k].push(T.stat[k]);
    }
    if (!dyn.RMmax[0]) dyn.RMmax = null;
    prog.done += ntr;
    setProg(prog.done / prog.total);
    return { dyn, stat, phi: BD.dynamicFactor(stat, dyn) };
}

let stopTimer = null;
function runStart() {
    ui.busy = true; ui.stopReq = false;
    ui.expOn = enabled('sExp');
    setEnable(['sRst', 'sOpen', 'sExp', 'sRep', 'sPrp'], false);
    $('fig').classList.add('busy');
    const t0 = performance.now();
    stopTimer = setInterval(() => {                 // Run becomes Stop after 5 s
        if (performance.now() - t0 > 5000 && $('sRun').dataset.state !== 'stop') {
            setIcon('sRun', 'stop'); $('sRun').dataset.state = 'stop';
            $('sRun').dataset.tip = 'Stop: ends the run, keeps the input'; hideTip();
        }
    }, 200);
}
function runEnd() {
    clearInterval(stopTimer);
    ui.busy = false; ui.abort = null;
    setIcon('sRun', 'run'); $('sRun').dataset.state = ''; $('sRun').dataset.tip = SYMS.sRun[1];
    $('fig').classList.remove('busy');
    setEnable(['sRun', 'sRst', 'sOpen'], true);
}
function runFailed(e) {
    setProg(null);
    if (e.message === 'stopped') {
        setMsg('stopped');
        if (ui.expOn) setEnable(['sExp', 'sRep', 'sPrp'], true);
    } else { setMsg('error', '#b30000'); alert('Bridge dynamics: ' + e.message); }
}
function onStop() {
    if ($('sRun').dataset.state !== 'stop') return;
    ui.stopReq = true; setMsg('stopping ...');
    if (ui.abort) ui.abort();
}
const nextFrame = () => new Promise(r => setTimeout(r, 0));

async function onRun() {
    if (!updateModel()) return;
    if (envMode()) return runEnvelope();
    runStart();
    try {
        const t0 = performance.now();
        prog.done = 0; prog.total = inp.load === 'HSLM-A' ? 10 : 1;
        setProg(0); await nextFrame();
        const { dyn, phi } = await runCase(inp.load);
        if (isEnv) { runs = []; isEnv = false; }
        runs.push({ load: inp.load, v: inp.v.slice(), d: dyn.dmax, a: dyn.amax, M: dyn.Mmax, Mn: dyn.Mmin, V: dyn.Vmax,
                    pd: phi.dmax, pM: phi.Mmax, pMn: phi.Mmin, pV: phi.Vmax, r: dyn.rmax, pr: phi.rmax,
                    fpp: ppFactor(inp.v), dyn, phi, f: mdl.f });
        redrawResults();
        setMsg(`done in ${((performance.now() - t0) / 1000).toFixed(1)} sec`);
        setEnable(['sExp', 'sRep', 'sPrp'], true);
        window.out = { runs, mdl };
    } catch (e) { runFailed(e); }
    runEnd();
}

// HSLM-A against all A1-E5 line categories, built on Run with train 'All'
async function runEnvelope() {
    runStart();
    try {
        const t0 = performance.now();
        prog.done = 0; prog.total = 10 + catList.length;
        setProg(0); await nextFrame();
        let { dyn, phi } = await runCase('HSLM-A');
        const nv = inp.v.length;
        const E = { v: inp.v.slice(), hdyn: dyn, hphi: phi, cat: [] };
        Object.assign(E, { hd: maxrow(dyn.dmax), ha: maxrow(dyn.amax), hM: maxrow(dyn.Mmax), hMn: minrow(dyn.Mmin), hV: maxrow(dyn.Vmax),
            hpd: maxrow(phi.dmax), hpM: maxrow(phi.Mmax), hpMn: maxrow(phi.Mmin), hpV: maxrow(phi.Vmax), hr: maxrow(dyn.rmax), hpr: maxrow(phi.rmax) });
        const ninf = () => new Float64Array(nv).fill(-Infinity);
        Object.assign(E, { od: ninf(), oa: ninf(), oM: ninf(), oMn: new Float64Array(nv).fill(Infinity), oV: ninf(),
                           opd: ninf(), opM: ninf(), opMn: ninf(), opV: ninf(), or: ninf(), opr: ninf() });
        const mx = (a, b) => Float64Array.from(a, (x, j) => Math.max(x, b[j])), mn = (a, b) => Float64Array.from(a, (x, j) => Math.min(x, b[j]));
        for (const c of catList) {
            ({ dyn, phi } = await runCase(c));
            const r = { name: c, d: maxrow(dyn.dmax), a: maxrow(dyn.amax), M: maxrow(dyn.Mmax), Mn: minrow(dyn.Mmin), V: maxrow(dyn.Vmax),
                        pd: maxrow(phi.dmax), pM: maxrow(phi.Mmax), pMn: maxrow(phi.Mmin), pV: maxrow(phi.Vmax), r: maxrow(dyn.rmax), pr: maxrow(phi.rmax) };
            E.od = mx(E.od, r.d); E.oa = mx(E.oa, r.a); E.or = mx(E.or, r.r); E.opr = mx(E.opr, r.pr);
            E.oM = mx(E.oM, r.M); E.oMn = mn(E.oMn, r.Mn); E.oV = mx(E.oV, r.V);
            E.opd = mx(E.opd, r.pd); E.opM = mx(E.opM, r.pM); E.opMn = mx(E.opMn, r.pMn); E.opV = mx(E.opV, r.pV);
            E.cat.push(r);
        }
        E.fpp = ppFactor(inp.v);
        E.f = mdl.f;
        envData = E; isEnv = true;
        redrawResults();
        setMsg(`done in ${((performance.now() - t0) / 1000).toFixed(1)} sec`);
        setEnable(['sExp', 'sRep', 'sPrp'], true);
        window.out = { envData, mdl };
    } catch (e) { runFailed(e); }
    runEnd();
}

function onClear() {
    runs = []; envData = null; isEnv = false;
    setMsg(''); setProg(null);
    setEnable(['sExp', 'sRep', 'sPrp'], false);
    redrawResults();
}
function onReset() {
    setDefaults(); sup = [];
    ui.cumOn = false; setIcon('sCum', 'cum');
    $('tSep').checked = false; $('tD').checked = true; $('tR').checked = false; $('tPhi').checked = false;
    $('tPbis').checked = true; $('tMx').checked = true; $('tMn').checked = false;
    onClear();
    updateModel();
}

// ============================================================ open / save
// window.__saveHook(blob, name) / window.__printHook(frame, title): test hooks
// (test/gui_test_web.html) that replace the file and print dialogs
async function saveFile(blob, name, desc, ext, mime) {
    if (window.__saveHook) return window.__saveHook(blob, name);
    if (window.showSaveFilePicker) {
        try {
            const hd = await window.showSaveFilePicker({ suggestedName: name, types: [{ description: desc, accept: { [mime]: [ext] } }] });
            const ws = await hd.createWritable(); await ws.write(blob); await ws.close();
            return hd.name;
        } catch (e) { if (e && e.name === 'AbortError') return null; }
    }
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob); a.download = name; a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 5000);
    return name;
}
async function onSaveInput() {
    const S = { L: $('iL').value, m: $('im').value, EI: $('iEI').value, xi: $('ixi').value, xiType: ui.xiType,
                Nmod: $('iNmod').value, NmodAuto: ui.autoN, fmax: $('ifmax').value, fmaxAuto: ui.autoF, dL: $('idL').value,
                load: $('iLoad').value, N_car: $('iNcar').value, v: $('iv').value, bQType: ui.bQType, b_Q: ui.bQVal,
                clamped: sup.slice() };
    const fn = await saveFile(new Blob([JSON.stringify(S, null, 2)], { type: 'application/json' }),
                              'BeamDyn_input.json', 'JSON input', '.json', 'application/json');
    if (fn) setMsg('input saved: ' + fn);
}
$('fileIn').addEventListener('change', async ev => {
    const file = ev.target.files[0];
    ev.target.value = '';
    if (!file) return;
    let S;
    try { S = JSON.parse(await file.text()); } catch (e) { setMsg('open failed', '#b30000'); alert(e.message); return; }
    const str = v => (Array.isArray(v) ? v.join(' ') : String(v));
    const flds = { L: 'iL', m: 'im', EI: 'iEI', xi: 'ixi', Nmod: 'iNmod', fmax: 'ifmax', dL: 'idL', N_car: 'iNcar', v: 'iv' };
    if ('W' in S) { S.bQType = bQTypes[+S.W > 0 ? 3 : 0]; S.b_Q = S.W; }     // older files: W
    if ('b_Q' in S) ui.bQVal = str(S.b_Q);
    if (bQTypes.includes(S.bQType)) { ui.bQType = S.bQType; if (S.bQType !== 'Triangular') ui.bQLast = S.bQType; }
    if (xiTypes.includes(S.xiType)) { ui.xiType = S.xiType; if (S.xiType !== 'Other') ui.xiLast = S.xiType; }
    else if ('xi' in S) ui.xiType = 'Other';
    ui.autoN = 'NmodAuto' in S ? !!S.NmodAuto : !('Nmod' in S);
    ui.autoF = 'fmaxAuto' in S ? !!S.fmaxAuto : !('fmax' in S);
    for (const k in flds) if (k in S) $(flds[k]).value = str(S[k]);
    if (loadList.includes(S.load)) $('iLoad').value = S.load;
    sup = Array.isArray(S.clamped) ? S.clamped.map(Boolean) : [];
    onClear();
    if (updateModel()) setMsg('input opened: ' + file.name);
});

const defaultName = () => (isEnv ? 'Envelope' : runs.length ? runs[runs.length - 1].load : 'results');

// ============================================================ export (xlsx)
function inputRows(model) {
    return [['Quantity', 'Value', 'Unit'], ['Span lengths L', mdl.L.map(gfmt).join(' + '), 'm'], ['Mass m', inp.m, 'kg/m'],
            ['Flexural rigidity EI', inp.EI, 'Nm^2'], ['Damping ratio', inp.xi, '%'], ['Damping, EN 1991-2 Table 8.6', inp.xiType, '-'],
            ['Modes kept', inp.Nmod, '-'], ['Max frequency fmax', inp.fmax, 'Hz'], ['Element size dL', inp.dL, 'm'],
            ['Load model', model, '-'], ['Cars', inp.N_car, '-'], ['Speed range', $('iv').value.trim(), 'km/h'],
            ['Load distribution', inp.loadDist, '-'], ['Triangle base b_Q', inp.b_Q, 'm']];
}
function modeRows() {
    const x = Array.from({ length: mdl.EL + 1 }, (_, i) => i * mdl.dL), nm = mdl.f.length;
    const ph = Array.from({ length: nm }, (_, n) => BD.modeShape(mdl, mdl.Phi, n, x));
    const rows = [['Natural frequencies (Hz):', ...Array.from(mdl.f, f => Math.round(f * 100) / 100)], [],
                  ['x (m)', ...ph.map((_, n) => 'mode ' + (n + 1))]];
    x.forEach((xv, i) => rows.push([xv, ...ph.map(p => p[i])]));
    return rows;
}
function onExport() { withPbis(() => exportXlsx()); }
function exportXlsx() {
    const v = inp.v, nv = v.length, sheets = [];
    if (isEnv) {
        const E = envData;
        sheets.push({ name: 'Input', rows: inputRows('All: HSLM-A and A1-E5') }, { name: 'ModeShapes', rows: modeRows() });
        const cols = [['a_max', 'a', 1], ['d_max', 'd', 1e3], ['M_max', 'M', 1e-6], ['M_min', 'Mn', 1e-6], ['V_max', 'V', 1e-6],
                      ['phi d_max', 'pd', 1], ['phi M_max', 'pM', 1], ['phi M_min', 'pMn', 1], ['phi V_max', 'pV', 1]];
        const hdr = ['v (km/h)'], D = v.map(x => [x]);
        for (const [lab, f, s] of cols) for (const [p, t] of [['h', ' HSLM-A'], ['o', ' A1-E5']]) {
            const y = E[p + f]; if (!y) continue;
            hdr.push(lab + t); y.forEach((val, i) => D[i].push(val * s));
        }
        sheets.push({ name: 'Results', rows: [hdr, ...D] });
        const DH = ['train', 'd_max (mm)', 'v (dmax)', 'a_max (m/s^2)', 'v (amax)', 'M_max (MNm)', 'v (Mmax)', 'M_min (MNm)', 'v (Mmin)', 'V_max (MN)', 'v (Vmax)'];
        const rows = [DH];
        for (const c of E.cat) {
            const row = [c.name];
            for (const [f, s, useMin] of [['d', 1e3, false], ['a', 1, false], ['M', 1e-6, false], ['Mn', 1e-6, true], ['V', 1e-6, false]]) {
                const y = Array.from(c[f]), pk = useMin ? Math.min(...y) : Math.max(...y);
                row.push(Math.round(pk * s * 100) / 100, v[y.indexOf(pk)]);
            }
            rows.push(row);
        }
        sheets.push({ name: 'TrainDetails', rows });
    } else {
        const rr = runs[runs.length - 1], sep = sepShown();
        sheets.push({ name: 'Input', rows: inputRows(rr.load) }, { name: 'ModeShapes', rows: modeRows() });
        const hdr = ['v (km/h)'], D = v.map(x => [x]);
        const qcols = (Y, label, isMin, s) => {
            if (!Y) return;
            if (sep && Y.length > 1) Y.forEach((row, t) => { hdr.push(`${label} A${t + 1}`); row.forEach((y, i) => D[i].push(y * s)); });
            else { const e = isMin ? minrow(Y) : maxrow(Y); hdr.push(label); e.forEach((y, i) => D[i].push(y * s)); }
        };
        qcols(rr.a, 'a_max', false, 1); qcols(rr.d, 'd_max', false, 1e3); qcols(rr.M, 'M_max', false, 1e-6);
        qcols(rr.Mn, 'M_min', true, 1e-6); qcols(rr.V, 'V_max', false, 1e-6);
        qcols(rr.pd, 'phi d_max', false, 1); qcols(rr.pM, 'phi M_max', false, 1); qcols(rr.pMn, 'phi M_min', false, 1); qcols(rr.pV, 'phi V_max', false, 1);
        const rows = [hdr, ...D, []];
        rows.push(['Quantity', 'Peak', 'at v (km/h)', "peak phi'", 'at v (km/h)']);
        const spec = [['d_max (mm)', 'd', 1e3, false, 'pd'], ['a_max (m/s^2)', 'a', 1, false, null], ['M_max (MNm)', 'M', 1e-6, false, 'pM'],
                      ['M_min (MNm)', 'Mn', 1e-6, true, 'pMn'], ['V_max (MN)', 'V', 1e-6, false, 'pV']];
        for (const [qn, f, s, useMin, pf] of spec) {
            const env = Array.from(useMin ? minrow(rr[f]) : maxrow(rr[f])), pk = useMin ? Math.min(...env) : Math.max(...env);
            let pc = '-', vc = '-';
            if (pf && rr[pf]) { const pe = Array.from(maxrow(rr[pf])), pp = Math.max(...pe.filter(isFinite)); pc = Math.round(pp * 1000) / 1000; vc = v[pe.indexOf(pp)]; }
            rows.push([qn, Math.round(pk * s * 100) / 100, v[env.indexOf(pk)], pc, vc]);
        }
        sheets.push({ name: 'Results', rows });
    }
    saveFile(XLSX.write(sheets), defaultName() + '.xlsx', 'Excel workbook', '.xlsx',
             'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet').then(fn => { if (fn) setMsg('saved'); });
}

// ============================================================ print
// an A4 page in a hidden frame, printed with the browser (Save as PDF); page
// margin boxes for the date (top right) and the page number "N (M)"
function printHTML(title, body) {
    const today = new Date().toISOString().slice(0, 10);
    const css = `@page { size: A4; margin: 25mm 25mm 25mm 25mm;
        @top-right { content: "${today}"; font: 11pt "Times New Roman", Times, serif; }
        @bottom-center { content: counter(page) " (" counter(pages) ")"; font: 11pt "Times New Roman", Times, serif; } }
      body { margin: 0; font: 11pt/1.3 "Times New Roman", Times, serif; color: #000; }
      h1 { font-size: 16pt; font-weight: normal; margin: 0 0 6pt; }
      h2 { font-size: 14pt; margin: 14pt 0 4pt; }
      p { margin: 0 0 6pt; text-align: justify; hyphens: auto; }
      table { border-collapse: collapse; margin: 4pt auto 10pt; }
      caption { caption-side: top; padding-bottom: 6pt; }
      th, td { padding: 1pt 7pt; text-align: left; font-weight: normal; }
      td.r, th.r { text-align: right; } td.c, th.c { text-align: center; }
      thead tr { border-top: 0.8pt solid #000; border-bottom: 0.8pt solid #000; }
      tbody tr:last-child { border-bottom: 0.8pt solid #000; }
      figure { margin: 6pt 0 10pt; text-align: center; break-inside: avoid; }
      figcaption { margin-top: 6pt; }
      .np { break-before: page; }
      svg { display: block; margin: 0 auto; }`;
    const fr = document.createElement('iframe');
    fr.style.cssText = 'position:fixed; right:0; bottom:0; width:0; height:0; border:0';
    document.body.appendChild(fr);
    const doc = fr.contentDocument;
    doc.open();
    doc.write(`<!doctype html><html><head><meta charset="utf-8"><title>${title}</title><style>${css}</style></head><body>${body}</body></html>`);
    doc.close();
    if (window.__printHook) { window.__printHook(fr, title); return; }
    const t0 = document.title;
    document.title = title;
    setTimeout(() => {
        fr.contentWindow.focus(); fr.contentWindow.print();
        document.title = t0;
        setTimeout(() => fr.remove(), 1000);
    }, 300);
}
// the four result plots in an SVG of w x h (px), plot height hp
function resultsSVG(w, xL, wp, hp, gap, top) {
    const h = top + 4 * hp + 3 * gap + 50;
    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.setAttribute('viewBox', `0 0 ${w} ${h}`); svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    const A = [0, 1, 2, 3].map(k => new PL.Axes(svg, { x: xL, y: top + k * (hp + gap), w: wp, h: hp }, { xTickLabels: k === 3, pick: false }));
    A[0].ylabel = AX[0].ylabel; A[3].xlabel = AX[3].xlabel;
    redrawResults(A);
    A.forEach((ax, k) => {                          // keep the GUI's datatips
        ax.tips = AX[k].tips.map(t => ({ s: ax.series[AX[k].series.indexOf(t.s)], i: t.i })).filter(t => t.s);
        ax.onTipText = AX[k].onTipText; ax.draw();
    });
    return { svg, h };
}
function onPrintPlots() {
    const { svg, h } = resultsSVG(560, 64, 480, 140, 14, 8);
    svg.setAttribute('width', '160mm'); svg.setAttribute('height', (160 * h / 560).toFixed(1) + 'mm');
    printHTML(defaultName() + '_plots', `<div style="margin-top:4mm">${svg.outerHTML}</div>`);
    setMsg('print dialog: Save as PDF');
}

function onReport() { withPbis(() => buildReport()); }
function buildReport() {
    const H = PL.texHTML, fr = mdl.f, Lphi = BD.EN.LPhi(mdl.L), clamped = sup.map((c, i) => (c ? i + 1 : 0)).filter(Boolean);
    const num = x => gfmt(Number(x.toPrecision(3)));
    const sci = x => { const p = Math.floor(Math.log10(Math.abs(x))); return p >= -2 && p <= 3 ? num(x) : `${gfmt(Number((x / 10 ** p).toFixed(4)))}&thinsp;·&thinsp;10<sup>${p}</sup>`; };
    const table = (cap, head, rows, al) => `<table><caption>Table ${++tn}:&ensp;${cap}</caption><thead><tr>${head.map((h, i) => `<th class="${al[i]}">${h}</th>`).join('')}</tr></thead>` +
        `<tbody>${rows.map(r => `<tr>${r.map((c, i) => `<td class="${al[i]}">${c}</td>`).join('')}</tr>`).join('')}</tbody></table>`;
    let tn = 0, b = '<h1>Railway Bridge Dynamics</h1>';

    // basis
    b += '<h2>Basis</h2><p>Dynamic analysis of a railway bridge for moving axle loads according to EN 1991-2:2024, 8.4.6: ' +
         'peak deck acceleration, displacement, end rotation and section forces as functions of the train speed, and the dynamic factor ' +
         H("$\\varphi' = y_\\mathrm{dyn}/y_\\mathrm{stat} - 1$") + ' of the same train. Vehicle/structure mass interaction is neglected (8.4.6.4(2)). ';
    if (tick('tPbis')) {
        b += `All results except ${H("$\\varphi'$")} include the increase for track defects and vehicle imperfections for carefully maintained track ` +
             `(8.4.6.5(3)–(4)): load effects ${H("$(1 + \\varphi'_\\mathrm{dyn} + \\varphi''/2)\\,y_\\mathrm{stat}$")}, accelerations times ${H("$1 + \\varphi''/2$")}, ` +
             `with ${H("$\\varphi''$")} = ${BD.EN.phiBis(Math.max(...inp.v), fr[0], Lphi).toFixed(2)} from Annex C, Formula (C.6) ` +
             `(${H('$n_0$')} = ${fr[0].toFixed(2)} Hz, ${H('$L_\\Phi$')} = ${num(Lphi)} m, ${H('$\\alpha_\\mathrm{v}$')} = 1 for ${H('$v$')} &gt; 22 m/s). `;
    } else b += `The increase for track defects and vehicle imperfections, ${H("$\\varphi''$")} according to Annex C (8.4.6.5(3)), is not included. `;
    b += inp.xiType === 'Other' ? `Damping ${H('$\\zeta$')} = ${gfmt(inp.xi)} %, given. `
        : `Damping ${H('$\\zeta$')} = ${gfmt(inp.xi)} % from Table 8.6 (${inp.xiType.toLowerCase()} bridge) with ${H('$L$ = $L_\\Phi$')} = ${num(BD.EN.LPhi(inp.L))} m (Table 8.2, 5.2). `;
    b += `Modes up to ${H('$f_\\mathrm{max}$')} = ${gfmt(inp.fmax)} Hz, ` + (ui.autoF
        ? 'at least the greatest of 30 Hz, 1.5 times the fundamental frequency and the third mode frequency (EN 1990:2023, A.2.9.4.2.1(3)).</p>'
        : 'given.</p>');

    // structural model
    const stx = clamped.length ? `support${clamped.length > 1 ? 's' : ''} ${clamped.join(', ')} clamped, the others pinned` : 'all supports pinned';
    const D = BD.loadDistribution(inp.loadDist, inp.b_Q, inp.dL), w = D.w.map(x => (100 * x).toFixed(1)).join('/');
    const ld = { '1-point': ' at the axles',
                 '3-point': `, each axle load distributed by the rail over three sleepers at 0.6 m, ${w} % (EN 1991-2:2024, 8.3.6.2)`,
                 '5-point': `, each axle load distributed over five sleepers at 0.6 m as the deflection of the rail on an elastic (Winkler) foundation, rail 60E1, track stiffness 100 MN/m per rail, ${w} %`,
                 'Triangular': `, each axle load distributed as a triangle of base ${H('$b_Q$')} = ${num(inp.b_Q)} m` }[inp.loadDist];
    b += `<h2>Structural model</h2><p>2D Bernoulli-Euler beam with uniform mass and stiffness, ${mdl.EL} elements of ${num(mdl.dL)} m, ${stx}. ` +
         `Modal superposition with ${inp.Nmod} modes, each modal equation integrated exactly for a piecewise-linear load (Nigam-Jennings), ` +
         'with the static correction of the truncated modes (mode acceleration) for displacement, rotation and section forces. ' +
         `Axle loads as moving forces${ld}. Model and mode shapes in figure 1, natural frequencies in table 3.</p>`;
    const msvg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    msvg.setAttribute('viewBox', '0 0 700 210'); msvg.setAttribute('width', '160mm'); msvg.setAttribute('height', '48mm');
    drawModel(msvg, { x: 20, y: 5, w: 660, h: 200 }, { top: true, legendSize: 13 });
    b += `<figure>${msvg.outerHTML}<figcaption>Figure 1:&ensp;Model with mode shapes, dimensions in m.</figcaption></figure>`;

    // input
    b += '<h2 class="np">Input</h2>';
    b += table('Bridge and analysis.', ['', 'Value', 'Unit'], [
        [`Spans, ${H('$L$')}`, mdl.L.map(gfmt).join(' + '), 'm'], [`Mass, ${H('$m$')}`, num(inp.m), 'kg/m'],
        [`Flexural rigidity, ${H('$EI$')}`, sci(inp.EI), 'Nm<sup>2</sup>'], [`Damping, ${H('$\\zeta$')}`, gfmt(inp.xi), '%'],
        [`Number of modes, ${H('$N_\\mathrm{mod}$')}`, inp.Nmod, '–'], [`Maximum frequency, ${H('$f_\\mathrm{max}$')}`, gfmt(inp.fmax), 'Hz'],
        [`Element length, ${H('$\\mathrm{d}L$')}`, num(mdl.dL), 'm']], ['', 'r', '']);
    const lm = isEnv ? 'HSLM-A (A1–A10) and EN 15528 A1–E5 (envelope)' : inp.load === 'HSLM-A' ? 'HSLM-A (A1–A10)' : `EN 15528, line category ${inp.load}`;
    const nc = isEnv ? `${inp.N_car} (A1–E5)` : inp.load === 'HSLM-A' ? '–' : String(inp.N_car);
    const trows = [['Load model', lm, '–'], ['Cars', nc, '–'], [`Speeds, ${H('$v$')}`, $('iv').value.trim(), 'km/h'], ['Load distribution', inp.loadDist.toLowerCase(), '–']];
    if (inp.loadDist === 'Triangular') trows.push([`Triangle base, ${H('$b_Q$')}`, num(inp.b_Q), 'm']);
    b += table('Train load.', ['', 'Value', 'Unit'], trows, ['', 'r', '']);
    b += table('Natural frequencies.', ['Mode', H('$f$ (Hz)')], Array.from(fr, (f, n) => [n + 1, f.toFixed(2)]), ['c', 'c']);

    // results
    b += '<h2 class="np">Results</h2><p>Peak values over the bridge and the speed range, with the speed and the governing train; ' +
         `${H("$\\varphi'$")} is the largest dynamic factor of the quantity. Diagrams in figure 2.</p>`;
    const head = ['', 'Peak', H('$v$ (km/h)'), 'train', H("$\\varphi'$"), H('$v_{\\varphi}$ (km/h)'), 'train'], al = ['', 'r', 'r', 'r', 'r', 'r', 'r'];
    const dash = s => (s === '-' || s === '' || s === 'NaN' ? '—' : String(s));   // negative: hyphen, as MATLAB
    if (isEnv) {
        b += table('Peak values, HSLM-A.', head, peakRows(envData.hdyn, envData.hphi, true).map(r => r.map(dash)), al);
        b += table('Peak values, EN 15528 A1–E5.', head, rowsA1E5().map(r => r.map(dash)), al);
        const dh = ['train', H('$d_\\mathrm{max}$'), H('$v_d$'), H('$a_\\mathrm{max}$'), H('$v_a$'), H('$M_\\mathrm{max}$'), H('$v_M$'),
                    H('$M_\\mathrm{min}$'), H('$v_{Mn}$'), H('$V_\\mathrm{max}$'), H('$v_V$')];
        b += table(`Peak values per line category (mm, m/s<sup>2</sup>, kNm, kN; ${H('$v$')} in km/h).`, dh, rowsDetail().map(r => r.map(dash)), ['', ...Array(10).fill('r')]);
    } else {
        const rr = runs[runs.length - 1];
        b += table('Peak values.', head, peakRows(rr.dyn, rr.phi, rr.d.length > 1, rr).map(r => r.map(dash)), al);
    }
    const { svg, h } = resultsSVG(700, 80, 590, 130, 12, 6);
    svg.setAttribute('width', '160mm'); svg.setAttribute('height', (160 * h / 700).toFixed(1) + 'mm');
    b += `<figure>${svg.outerHTML}<figcaption>Figure 2:&ensp;Peak acceleration, displacement, bending moment and shear force as functions of the speed.</figcaption></figure>`;
    printHTML(defaultName(), b);
    setMsg('print dialog: Save as PDF');
}

// peak rows: per quantity the peak, its speed and the governing train, then
// phi'.  With the phi'' tick the stored run (rr) holds the scaled values.
function peakRows(dyn, phi, isHSLM, rr) {
    const v = inp.v, H = PL.texHTML;
    const qn = ['$d_\\mathrm{max}$', '$a_\\mathrm{max}$', '$M_\\mathrm{max}$', '$M_\\mathrm{min}$', '$V_\\mathrm{max}$'];
    const us = ['mm', 'm/s<sup>2</sup>', 'kNm', 'kNm', 'kN'], dec = [1, 1, 0, 0, 0], sc = [1e3, 1, 1e-3, 1e-3, 1e-3];
    const useMin = [false, false, false, true, false];
    const raw = rr ? [rr.d, rr.a, rr.M, rr.Mn, rr.V] : [dyn.dmax, dyn.amax, dyn.Mmax, dyn.Mmin, dyn.Vmax];
    const pf = [phi.dmax, null, phi.Mmax, phi.Mmin, phi.Vmax];
    const rows = [];
    for (let j = 0; j < 5; j++) {
        const Y = raw[j]; if (!Y) continue;
        const pk = argExt(Y, useMin[j]);
        let pC = '-', vC = '-', pt = '';
        if (pf[j]) { const pp = argExt(pf[j], false); if (pp) { pC = pp.val.toFixed(2); vC = gfmt(v[pp.j]); pt = isHSLM ? 'A' + (pp.t + 1) : ''; } }
        rows.push([`${H(qn[j])} (${us[j]})`, (pk.val * sc[j]).toFixed(dec[j]), gfmt(v[pk.j]), isHSLM ? 'A' + (pk.t + 1) : '', pC, vC, pt]);
    }
    return rows;
}
// extreme over [ntr][nv], finite values only: {val, t, j}
function argExt(Y, useMin) {
    let best = null;
    Y.forEach((row, t) => row.forEach((y, j) => {
        if (!isFinite(y)) return;
        if (!best || (useMin ? y < best.val : y > best.val)) best = { val: y, t, j };
    }));
    return best;
}
function rowsA1E5() {
    const v = inp.v, H = PL.texHTML, C = envData.cat;
    const qn = ['$d_\\mathrm{max}$', '$a_\\mathrm{max}$', '$M_\\mathrm{max}$', '$M_\\mathrm{min}$', '$V_\\mathrm{max}$'];
    const us = ['mm', 'm/s<sup>2</sup>', 'kNm', 'kNm', 'kN'], dec = [1, 1, 0, 0, 0], sc = [1e3, 1, 1e-3, 1e-3, 1e-3];
    const rf = ['d', 'a', 'M', 'Mn', 'V'], pf = ['pd', '', 'pM', 'pMn', 'pV'], useMin = [false, false, false, true, false];
    const catPeak = (f, mn) => {
        let best = null;
        C.forEach(c => Array.from(c[f]).forEach((y, j) => { if (isFinite(y) && (!best || (mn ? y < best.val : y > best.val))) best = { val: y, j, name: c.name }; }));
        return best;
    };
    return rf.map((f, j) => {
        const pk = catPeak(f, useMin[j]);
        let pC = '-', vC = '-', pn = '';
        if (pf[j]) { const pp = catPeak(pf[j], false); if (pp) { pC = pp.val.toFixed(2); vC = gfmt(v[pp.j]); pn = pp.name; } }
        return [`${H(qn[j])} (${us[j]})`, (pk.val * sc[j]).toFixed(dec[j]), gfmt(v[pk.j]), pk.name, pC, vC, pn];
    });
}
function rowsDetail() {
    const v = inp.v;
    return envData.cat.map(c => {
        const row = [c.name];
        for (const [f, s, mn, d] of [['d', 1e3, false, 1], ['a', 1, false, 1], ['M', 1e-3, false, 0], ['Mn', 1e-3, true, 0], ['V', 1e-3, false, 0]]) {
            const y = Array.from(c[f]), pk = mn ? Math.min(...y) : Math.max(...y);
            row.push((pk * s).toFixed(d), v[y.indexOf(pk)].toFixed(0));
        }
        return row;
    });
}

// ============================================================ start
setDefaults();
setEnable(['sExp', 'sRep', 'sPrp'], false);
updateModel();
})();

/* plot.js -------------------------------------------------------------------
   SVG plotting for the BeamDyn page: line plots in a box with inward ticks,
   a legend at the 'best' place, draggable datatips, and a small LaTeX subset
   for labels ($a_\mathrm{max}$ (m/s$^2$), $\varphi'$, ...) rendered as SVG
   tspans or HTML.  Vector output, so the printed report stays vector.

   Andreas Andersson, 2026-11-01 (assisted by Claude)
---------------------------------------------------------------------------*/
(function (root) {
'use strict';

const PL = {};
const NS = 'http://www.w3.org/2000/svg';
PL.FONT = '"Times New Roman", Times, serif';
// MATLAB's default colour order
PL.COLORS = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'];
PL.color = i => PL.COLORS[((i % 7) + 7) % 7];

PL.el = function (tag, attrs, parent) {
    const e = document.createElementNS(NS, tag);
    for (const k in attrs || {}) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
};

// ------------------------------------------------------------ LaTeX subset
// Runs {t, it, lvl}: italic inside $...$ for letters, \mathrm{} upright,
// _ and ^ (one character or {...}) one level down/up, a few symbols.
const SYM = { '\\varphi': 'φ', '\\phi': 'ϕ', '\\zeta': 'ζ', '\\alpha': 'α', '\\beta': 'β',
              '\\Phi': 'Φ', '\\cdot': '·', '\\%': '%', '\\,': ' ' };
PL.texRuns = function (s) {
    const runs = [];
    let math = false, i = 0;
    const push = (t, it, lvl) => {
        if (!t) return;
        const last = runs[runs.length - 1];
        if (last && last.it === it && last.lvl === lvl) last.t += t; else runs.push({ t, it, lvl });
    };
    const group = (j) => {                       // {...} or one character at j
        if (s[j] === '{') { let d = 1, k = j + 1; while (k < s.length && d) { if (s[k] === '{') d++; else if (s[k] === '}') d--; k++; } return [s.slice(j + 1, k - 1), k]; }
        if (s[j] === '\\') { const m = /^\\[a-zA-Z]+/.exec(s.slice(j)); return [m[0], j + m[0].length]; }
        return [s[j], j + 1];
    };
    const emit = (txt, lvl, upright) => {        // text inside math
        let k = 0;
        while (k < txt.length) {
            const c = txt[k];
            if (c === '\\') {
                const m = /^\\([a-zA-Z]+|.)/.exec(txt.slice(k));
                const cmd = m[0];
                if (cmd === '\\mathrm') { const [g, k2] = groupIn(txt, k + 7); emit(g, lvl, true); k = k2; continue; }
                push(SYM[cmd] !== undefined ? SYM[cmd] : cmd.slice(1), !upright && /^\\[a-zA-Z]/.test(cmd) && cmd !== '\\cdot', lvl);
                k += cmd.length; continue;
            }
            if (c === '_' || c === '^') {
                const [g, k2] = groupIn(txt, k + 1);
                emit(g, lvl + (c === '_' ? -1 : 1), upright); k = k2; continue;
            }
            if (c === "'") { push('′', false, lvl); k++; continue; }
            if (c === '-') { push('−', false, lvl); k++; continue; }
            if (c === '{' || c === '}') { k++; continue; }
            push(c, !upright && /[a-zA-Z]/.test(c), lvl); k++;
        }
    };
    function groupIn(txt, j) {
        const brace = (j0) => { let d = 1, k = j0 + 1; while (k < txt.length && d) { if (txt[k] === '{') d++; else if (txt[k] === '}') d--; k++; } return k; };
        if (txt[j] === '{') { const k = brace(j); return [txt.slice(j + 1, k - 1), k]; }
        if (txt[j] === '\\') {
            const m = /^\\([a-zA-Z]+|.)/.exec(txt.slice(j)), e = j + m[0].length;
            if (txt[e] === '{') { const k = brace(e); return [txt.slice(j, k), k]; }   // \mathrm{...}
            return [m[0], e];
        }
        return [txt[j] || '', j + 1];
    }
    void group;
    let buf = '';
    while (i < s.length) {
        if (s[i] === '$') {
            if (math) { emit(buf, 0, false); } else push(buf, false, 0);
            buf = ''; math = !math; i++; continue;
        }
        buf += s[i]; i++;
    }
    if (math) emit(buf, 0, false); else push(buf, false, 0);
    return runs;
};

const escH = t => t.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
PL.texHTML = function (s) {
    return PL.texRuns(s).map(r => {
        let h = escH(r.t);
        if (r.it) h = '<i>' + h + '</i>';
        if (r.lvl < 0) h = '<sub>' + h + '</sub>';
        if (r.lvl > 0) h = '<sup>' + h + '</sup>';
        return h;
    }).join('');
};

// SVG text from a LaTeX string; opt: size (px), anchor, rotate (deg), fill
PL.text = function (parent, s, x, y, opt) {
    opt = opt || {};
    const size = opt.size || 14.7;
    const t = PL.el('text', {
        x, y, 'font-family': PL.FONT, 'font-size': size, fill: opt.fill || '#000',
        'text-anchor': opt.anchor || 'start', 'dominant-baseline': opt.baseline || 'middle'
    }, parent);
    if (opt.rotate) t.setAttribute('transform', `rotate(${opt.rotate} ${x} ${y})`);
    let lvl = 0;
    for (const r of PL.texRuns(s)) {
        const ts = PL.el('tspan', {}, t);
        if (r.it) ts.setAttribute('font-style', 'italic');
        if (r.lvl !== lvl) {
            const dy = (lvl - r.lvl) * size * 0.3 * (r.lvl === 0 || lvl === 0 ? 1 : 1);
            ts.setAttribute('dy', dy.toFixed(2));
            lvl = r.lvl;
        }
        if (r.lvl) ts.setAttribute('font-size', (size * 0.72).toFixed(1));
        ts.textContent = r.t;
    }
    if (lvl) { const ts = PL.el('tspan', { dy: (lvl * size * 0.3).toFixed(2) }, t); ts.textContent = '​'; }
    return t;
};

// width of a LaTeX string in px (canvas measurement)
let mctx = null;
PL.textWidth = function (s, size) {
    size = size || 14.7;
    if (!mctx) mctx = document.createElement('canvas').getContext('2d');
    let w = 0;
    for (const r of PL.texRuns(s)) {
        mctx.font = (r.it ? 'italic ' : '') + (r.lvl ? size * 0.72 : size) + 'px ' + PL.FONT;
        w += mctx.measureText(r.t).width;
    }
    return w;
};

// ------------------------------------------------------------------- ticks
PL.ticks = function (lo, hi, n) {
    if (!(isFinite(lo) && isFinite(hi)) || hi <= lo) return { t: [lo], step: 1 };
    const raw = (hi - lo) / Math.max(n, 1);
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const nr = raw / mag;
    const step = (nr <= 1 ? 1 : nr <= 2 ? 2 : nr <= 2.5 ? 2.5 : nr <= 5 ? 5 : 10) * mag;
    const t = [];
    for (let k = Math.ceil(lo / step - 1e-9); k * step <= hi + step * 1e-9; k++) t.push(Math.abs(k * step) < step * 1e-9 ? 0 : k * step);
    return { t, step };
};
PL.fmtTick = function (v, step) {
    const dec = Math.max(0, -Math.floor(Math.log10(step) + 1e-9)) + 1;
    const s = String(Number(v.toFixed(Math.min(dec, 8))));   // 0, 0.5, 1 as MATLAB
    return v < 0 ? '−' + s.slice(1) : s;
};

// ------------------------------------------------------------------- Axes
// A plot box at (x, y, w, h) inside an <svg> group.
PL.Axes = class {
    constructor(parent, box, opt) {
        this.g = PL.el('g', {}, parent);
        this.box = box;
        this.opt = Object.assign({ size: 14.7, xTickLabels: true }, opt || {});
        this.xlim = [0, 1]; this.ylim = [0, 1];
        this.series = []; this.tips = [];
        this.xlabel = ''; this.ylabel = ''; this.legend = null;
        this.onTipText = null;
        const id = 'clip' + Math.random().toString(36).slice(2);
        const cp = PL.el('clipPath', { id }, this.g);
        PL.el('rect', { x: box.x, y: box.y, width: box.w, height: box.h }, cp);
        this.clipId = id;
    }
    X(x) { return this.box.x + (x - this.xlim[0]) / (this.xlim[1] - this.xlim[0]) * this.box.w; }
    Y(y) { return this.box.y + this.box.h - (y - this.ylim[0]) / (this.ylim[1] - this.ylim[0]) * this.box.h; }
    clear() { this.series = []; this.tips = []; this.legend = null; }
    plot(x, y, opt) { this.series.push(Object.assign({ x: Array.from(x), y: Array.from(y), color: PL.color(0), width: 1, name: '' }, opt || {})); }

    draw() {
        const g = this.g, b = this.box, sz = this.opt.size;
        while (g.childNodes.length > 1) g.removeChild(g.lastChild);   // keep the clipPath
        PL.el('rect', { x: b.x, y: b.y, width: b.w, height: b.h, fill: '#fff' }, g);
        // ticks inward, labels outside
        const tx = PL.ticks(this.xlim[0], this.xlim[1], Math.max(2, Math.round(b.w / 45)));
        const ty = PL.ticks(this.ylim[0], this.ylim[1], Math.max(2, Math.round(b.h / 32)));
        const tl = Math.min(b.w, b.h) * 0.025 + 2;
        const tick = (x1, y1, x2, y2) => PL.el('line', { x1, y1, x2, y2, stroke: '#000', 'stroke-width': 0.5 }, g);
        for (const v of tx.t) {
            const px = this.X(v);
            if (px < b.x - 0.5 || px > b.x + b.w + 0.5) continue;
            tick(px, b.y + b.h, px, b.y + b.h - tl); tick(px, b.y, px, b.y + tl);
            if (this.opt.xTickLabels) PL.text(g, PL.fmtTick(v, tx.step), px, b.y + b.h + sz * 0.75, { size: sz, anchor: 'middle' });
        }
        let ylw = 0;
        for (const v of ty.t) {
            const py = this.Y(v);
            if (py < b.y - 0.5 || py > b.y + b.h + 0.5) continue;
            tick(b.x, py, b.x + tl, py); tick(b.x + b.w, py, b.x + b.w - tl, py);
            const s = PL.fmtTick(v, ty.step);
            PL.text(g, s, b.x - 4, py, { size: sz, anchor: 'end' });
            ylw = Math.max(ylw, PL.textWidth(s, sz));
        }
        if (this.ylabel) PL.text(g, this.ylabel, b.x - ylw - 4 - sz * 0.75, b.y + b.h / 2, { size: sz, anchor: 'middle', rotate: -90 });
        if (this.xlabel) PL.text(g, this.xlabel, b.x + b.w / 2, b.y + b.h + sz * 1.9, { size: sz, anchor: 'middle' });
        // curves
        const cg = PL.el('g', { 'clip-path': `url(#${this.clipId})` }, g);
        this.series.forEach((s, k) => {
            let d = '', pen = false;
            for (let i = 0; i < s.x.length; i++) {
                const yv = s.y[i];
                if (!isFinite(yv)) { pen = false; continue; }
                d += (pen ? 'L' : 'M') + this.X(s.x[i]).toFixed(2) + ' ' + this.Y(yv).toFixed(2);
                pen = true;
            }
            if (s.x.length === 1 && isFinite(s.y[0]))                // one speed: a marker
                PL.el('circle', { cx: this.X(s.x[0]), cy: this.Y(s.y[0]), r: 2.5, fill: s.color }, cg);
            PL.el('path', { d, fill: 'none', stroke: s.color, 'stroke-width': s.width, 'stroke-linejoin': 'round' }, cg);
            if (this.opt.pick !== false) {
                const hit = PL.el('path', { d, fill: 'none', stroke: 'transparent', 'stroke-width': 8, style: 'cursor:pointer' }, cg);
                hit.addEventListener('mousedown', ev => this.pick(ev, k));
            }
        });
        PL.el('rect', { x: b.x, y: b.y, width: b.w, height: b.h, fill: 'none', stroke: '#000', 'stroke-width': 0.5 }, g);
        if (this.legend) this.drawLegend();
        this.tipLayer = PL.el('g', {}, g);
        this.tips.forEach(t => this.drawTip(t));
    }

    // legend inside the box at the corner with the fewest data points
    drawLegend() {
        const L = this.legend, sz = L.size || this.opt.size * 0.9, b = this.box;
        const names = L.names, nc = L.cols || 1, lw = 18, pad = 4;
        const rowH = sz * 1.15;
        const colW = [];
        names.forEach((n, i) => { const c = i % nc; colW[c] = Math.max(colW[c] || 0, PL.textWidth(n, sz) + lw + 6); });
        const W = colW.reduce((a, c) => a + c, 0) + 2 * pad, H = Math.ceil(names.length / nc) * rowH + 2 * pad;
        const m = 4, cand = [[b.x + b.w - W - m, b.y + m], [b.x + m, b.y + m], [b.x + m, b.y + b.h - H - m],
                             [b.x + b.w - W - m, b.y + b.h - H - m], [b.x + (b.w - W) / 2, b.y + m], [b.x + (b.w - W) / 2, b.y + b.h - H - m]];
        let best = cand[0], bestN = Infinity;
        for (const [cx, cy] of cand) {
            let n = 0;
            for (const s of this.series) for (let i = 0; i < s.x.length; i++) {
                const px = this.X(s.x[i]), py = this.Y(s.y[i]);
                if (px >= cx && px <= cx + W && py >= cy && py <= cy + H) n++;
            }
            if (n < bestN) { bestN = n; best = [cx, cy]; }
        }
        const g = PL.el('g', {}, this.g);
        PL.el('rect', { x: best[0], y: best[1], width: W, height: H, fill: '#fff', stroke: L.box ? '#000' : 'none', 'stroke-width': 0.5 }, g);
        names.forEach((n, i) => {
            const c = i % nc, r = Math.floor(i / nc);
            let x0 = best[0] + pad; for (let k = 0; k < c; k++) x0 += colW[k];
            const y0 = best[1] + pad + (r + 0.5) * rowH;
            PL.el('line', { x1: x0, y1: y0, x2: x0 + lw, y2: y0, stroke: L.colors[i], 'stroke-width': 1 }, g);
            PL.text(g, n, x0 + lw + 4, y0, { size: sz });
        });
    }

    // ---- datatips: press on a curve (and drag along it); press on a tip and
    // drag to move it, double-click to remove; shift adds another tip
    pick(ev, k) {
        ev.preventDefault(); ev.stopPropagation();
        const near = !ev.shiftKey && this.tips.find(t => Math.abs(this.X(t.s.x[t.i]) - this.mouseX(ev)) < 0.02 * this.box.w && t.s === this.series[k]);
        const tip = near || { s: this.series[k], i: 0 };
        if (!near) this.tips.push(tip);
        this.dragTip(tip, ev);
    }
    mouseX(ev) {
        const svg = this.g.ownerSVGElement, p = svg.createSVGPoint();
        p.x = ev.clientX; p.y = ev.clientY;
        return p.matrixTransform(svg.getScreenCTM().inverse()).x;
    }
    dragTip(tip, ev) {
        const move = e => {
            const xm = this.mouseX(e);
            let bi = 0, bd = Infinity;
            tip.s.x.forEach((xv, i) => { const d = Math.abs(this.X(xv) - xm); if (d < bd && isFinite(tip.s.y[i])) { bd = d; bi = i; } });
            tip.i = bi;
            this.redrawTips();
        };
        const up = () => { window.removeEventListener('mousemove', move); window.removeEventListener('mouseup', up); };
        window.addEventListener('mousemove', move); window.addEventListener('mouseup', up);
        move(ev);
    }
    redrawTips() {
        while (this.tipLayer.firstChild) this.tipLayer.removeChild(this.tipLayer.firstChild);
        this.tips.forEach(t => this.drawTip(t));
    }
    drawTip(t) {
        const x = t.s.x[t.i], y = t.s.y[t.i], px = this.X(x), py = this.Y(y), b = this.box;
        const lines = this.onTipText ? this.onTipText(x, y) : [x.toFixed(0), y.toPrecision(3)];
        const g = PL.el('g', { style: 'cursor:move' }, this.tipLayer);
        PL.el('circle', { cx: px, cy: py, r: 3, fill: '#000' }, g);
        const fs = 11, w = Math.max(...lines.map(s => { if (!mctx) PL.textWidth('', fs); mctx.font = fs + 'px system-ui, sans-serif'; return mctx.measureText(s).width; })) + 8;
        const h = lines.length * fs * 1.2 + 6;
        const right = px > b.x + 0.75 * b.w;
        let bx = right ? px - 6 - w : px + 6, by = py - h - 4;
        if (by < b.y) by = Math.min(py + 6, b.y + b.h - h);
        PL.el('rect', { x: bx, y: by, width: w, height: h, fill: '#fff', stroke: '#b3b3b3', 'stroke-width': 0.7 }, g);
        lines.forEach((s, i) => {
            const tt = PL.el('text', { x: bx + 4, y: by + 3 + (i + 0.5) * fs * 1.2, 'font-size': fs, 'font-family': 'system-ui, sans-serif', 'dominant-baseline': 'middle' }, g);
            tt.textContent = s;
        });
        g.addEventListener('mousedown', ev => { ev.preventDefault(); ev.stopPropagation(); this.dragTip(t, ev); });
        g.addEventListener('dblclick', ev => { ev.stopPropagation(); this.tips.splice(this.tips.indexOf(t), 1); this.redrawTips(); });
    }
};

root.PL = PL;
})(typeof self !== 'undefined' ? self : this);

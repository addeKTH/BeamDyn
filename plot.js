/* plot.js -------------------------------------------------------------------
   A small canvas plotting layer, sized to do exactly what the MATLAB figure
   did and nothing more: line plots with a box, ticks, legend and draggable
   datatips, plus the beam sketch of PlotBeam2D.m.

   Text supports a light TeX-like markup, _{...} and ^{...}, so axis labels
   read like the MATLAB ones (a_{max} (m/s^2), \varphi') without pulling in a
   maths typesetting library.
---------------------------------------------------------------------------*/
(function (root) {
'use strict';

const P = {};

// MATLAB's default ColorOrder, so the curves come out the familiar colours.
P.COLORS = [
    '#0072BD', '#D95319', '#EDB120', '#7E2F8E',
    '#77AC30', '#4DBEEE', '#A2142F'
];
P.color = (i) => P.COLORS[((i % P.COLORS.length) + P.COLORS.length) % P.COLORS.length];

const FONT = '11px "Segoe UI", system-ui, sans-serif';
const FONT_SM = '9px "Segoe UI", system-ui, sans-serif';

// ------------------------------------------------------------- text markup
// Split "M_{max} (kNm)" into runs of {t, lvl} where lvl is 0, -1 (sub) or +1.
function parseRuns(s) {
    const runs = []; let i = 0, buf = '';
    const push = (lvl) => { if (buf) { runs.push({ t: buf, lvl }); buf = ''; } };
    while (i < s.length) {
        const c = s[i];
        if ((c === '_' || c === '^') && s[i + 1] === '{') {
            push(0);
            const j = s.indexOf('}', i + 2);
            runs.push({ t: s.slice(i + 2, j), lvl: c === '_' ? -1 : 1 });
            i = j + 1;
        } else if ((c === '_' || c === '^') && i + 1 < s.length) {
            push(0);
            runs.push({ t: s[i + 1], lvl: c === '_' ? -1 : 1 });
            i += 2;
        } else { buf += c; i++; }
    }
    push(0);
    return runs;
}

P.textWidth = function (ctx, s, size) {
    size = size || 11;
    const runs = parseRuns(s); let w = 0;
    for (const r of runs) {
        ctx.font = (r.lvl ? Math.round(size * 0.78) : size) + 'px "Segoe UI", system-ui, sans-serif';
        w += ctx.measureText(r.t).width;
    }
    return w;
};

// align: 'left' | 'center' | 'right'; rot in radians (for y labels)
P.text = function (ctx, s, x, y, opt) {
    opt = opt || {};
    const size = opt.size || 11;
    const runs = parseRuns(s);
    const w = P.textWidth(ctx, s, size);
    ctx.save();
    ctx.translate(x, y);
    if (opt.rot) ctx.rotate(opt.rot);
    let dx = 0;
    if (opt.align === 'center') dx = -w / 2;
    else if (opt.align === 'right') dx = -w;
    ctx.fillStyle = opt.color || '#000';
    ctx.textBaseline = opt.baseline || 'middle';
    for (const r of runs) {
        const fs = r.lvl ? Math.round(size * 0.78) : size;
        ctx.font = fs + 'px "Segoe UI", system-ui, sans-serif';
        const dy = r.lvl === -1 ? size * 0.28 : (r.lvl === 1 ? -size * 0.35 : 0);
        ctx.fillText(r.t, dx, dy);
        dx += ctx.measureText(r.t).width;
    }
    ctx.restore();
    return w;
};

// ------------------------------------------------------------------- ticks
// 1-2-5 decade steps, about n intervals, as MATLAB would pick them.
P.ticks = function (lo, hi, n) {
    if (!(isFinite(lo) && isFinite(hi)) || hi <= lo) return [lo];
    const raw = (hi - lo) / (n || 5);
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const norm = raw / mag;
    const step = (norm <= 1.5 ? 1 : norm <= 3 ? 2 : norm <= 7 ? 5 : 10) * mag;
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) {
        out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
    }
    return out;
};

P.fmtTick = function (v) {
    if (v === 0) return '0';
    const a = Math.abs(v);
    if (a >= 1e5 || a < 1e-3) return v.toExponential(0).replace('e+', 'e');
    const s = (Math.round(v * 1e6) / 1e6).toString();
    return s;
};

// ------------------------------------------------------------------- Axes
// box is the plot rectangle in CSS pixels inside the canvas.
P.Axes = class Axes {
    constructor(canvas, box) {
        this.cv = canvas;
        this.box = box;
        this.xlim = [0, 1];
        this.ylim = [0, 1];
        this.rlim = null;               // right-hand ruler, when in use
        this.series = [];
        this.marks = [];                // datatips
        this.xlabel = ''; this.ylabel = ''; this.rlabel = '';
        this.showXTickLabels = true;
        this.legendEntries = null;
        this.visible = true;
        this.fit();
    }

    // Redraw into an offscreen canvas at a higher pixel ratio and return a PNG
    // data URL - used by the report, where the screen resolution is too coarse.
    snapshot(scale) {
        const off = document.createElement('canvas');
        off.style.width = this.cv.clientWidth + 'px';
        off.style.height = this.cv.clientHeight + 'px';
        Object.defineProperty(off, 'clientWidth', { value: this.cv.clientWidth });
        Object.defineProperty(off, 'clientHeight', { value: this.cv.clientHeight });
        const keep = this.cv, keepD = this.dprOverride;
        this.cv = off; this.dprOverride = scale || 3;
        this.draw();
        this.cv = keep; this.dprOverride = keepD;
        return off.toDataURL('image/png');
    }

    fit() {
        const dpr = this.dprOverride || window.devicePixelRatio || 1;
        const w = this.cv.clientWidth, h = this.cv.clientHeight;
        if (this.cv.width !== Math.round(w * dpr) || this.cv.height !== Math.round(h * dpr)) {
            this.cv.width = Math.round(w * dpr);
            this.cv.height = Math.round(h * dpr);
        }
        const ctx = this.cv.getContext('2d');
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        this.ctx = ctx;
    }

    clear() { this.series = []; this.marks = []; this.rlim = null; this.rlabel = ''; this.legendEntries = null; }

    add(x, y, opt) {
        opt = opt || {};
        this.series.push({
            x, y,
            color: opt.color || P.color(this.series.length),
            dash: opt.dash || null,
            width: opt.width || 1,
            name: opt.name || '',
            right: !!opt.right,
            noLegend: !!opt.noLegend
        });
        return this.series[this.series.length - 1];
    }

    // data -> pixel
    px(x) { const b = this.box; return b.x + (x - this.xlim[0]) / (this.xlim[1] - this.xlim[0]) * b.w; }
    py(y, right) {
        const b = this.box, l = right && this.rlim ? this.rlim : this.ylim;
        return b.y + b.h - (y - l[0]) / (l[1] - l[0]) * b.h;
    }
    xData(px) { const b = this.box; return this.xlim[0] + (px - b.x) / b.w * (this.xlim[1] - this.xlim[0]); }

    autoY(pad) {
        let lo = Infinity, hi = -Infinity, rlo = Infinity, rhi = -Infinity;
        for (const s of this.series) {
            for (let i = 0; i < s.y.length; i++) {
                const v = s.y[i];
                if (!isFinite(v)) continue;
                if (s.right) { if (v < rlo) rlo = v; if (v > rhi) rhi = v; }
                else { if (v < lo) lo = v; if (v > hi) hi = v; }
            }
        }
        this.ylim = padLim(lo, hi, pad);
        this.rlim = isFinite(rlo) ? padLim(rlo, rhi, pad) : null;
    }

    draw() {
        this.fit();
        const ctx = this.ctx, b = this.box;
        ctx.clearRect(0, 0, this.cv.clientWidth, this.cv.clientHeight);
        if (!this.visible) return;

        // curves, clipped to the box
        ctx.save();
        ctx.beginPath(); ctx.rect(b.x, b.y, b.w, b.h); ctx.clip();
        for (const s of this.series) {
            ctx.beginPath();
            ctx.strokeStyle = s.color; ctx.lineWidth = s.width;
            ctx.setLineDash(s.dash || []);
            let started = false;
            for (let i = 0; i < s.x.length; i++) {
                const v = s.y[i];
                if (!isFinite(v)) { started = false; continue; }
                const X = this.px(s.x[i]), Y = this.py(v, s.right);
                if (!started) { ctx.moveTo(X, Y); started = true; } else ctx.lineTo(X, Y);
            }
            ctx.stroke();
        }
        ctx.setLineDash([]);
        ctx.restore();

        // frame
        ctx.strokeStyle = '#262626'; ctx.lineWidth = 0.8;
        ctx.strokeRect(b.x + 0.5, b.y + 0.5, b.w, b.h);

        // ticks
        const xt = P.ticks(this.xlim[0], this.xlim[1], 5);
        ctx.beginPath();
        for (const t of xt) {
            const X = Math.round(this.px(t)) + 0.5;
            if (X < b.x - 1 || X > b.x + b.w + 1) continue;
            ctx.moveTo(X, b.y + b.h); ctx.lineTo(X, b.y + b.h - 4);
            ctx.moveTo(X, b.y); ctx.lineTo(X, b.y + 4);
            if (this.showXTickLabels)
                P.text(ctx, P.fmtTick(t), X, b.y + b.h + 9, { align: 'center', size: 10 });
        }
        const yt = P.ticks(this.ylim[0], this.ylim[1], 4);
        for (const t of yt) {
            const Y = Math.round(this.py(t)) + 0.5;
            if (Y < b.y - 1 || Y > b.y + b.h + 1) continue;
            ctx.moveTo(b.x, Y); ctx.lineTo(b.x + 4, Y);
            if (!this.rlim) { ctx.moveTo(b.x + b.w, Y); ctx.lineTo(b.x + b.w - 4, Y); }
            P.text(ctx, P.fmtTick(t), b.x - 5, Y, { align: 'right', size: 10 });
        }
        if (this.rlim) {
            const rt = P.ticks(this.rlim[0], this.rlim[1], 4);
            for (const t of rt) {
                const Y = Math.round(this.py(t, true)) + 0.5;
                if (Y < b.y - 1 || Y > b.y + b.h + 1) continue;
                ctx.moveTo(b.x + b.w, Y); ctx.lineTo(b.x + b.w - 4, Y);
                P.text(ctx, P.fmtTick(t), b.x + b.w + 5, Y, { align: 'left', size: 10 });
            }
        }
        ctx.stroke();

        // labels
        if (this.xlabel) P.text(ctx, this.xlabel, b.x + b.w / 2, b.y + b.h + 23, { align: 'center' });
        if (this.ylabel) P.text(ctx, this.ylabel, b.x - 40, b.y + b.h / 2, { align: 'center', rot: -Math.PI / 2 });
        if (this.rlabel) P.text(ctx, this.rlabel, b.x + b.w + 42, b.y + b.h / 2, { align: 'center', rot: -Math.PI / 2 });

        this.drawLegend();
        this.drawMarks();
    }

    drawLegend() {
        const ents = this.series.filter(s => s.name && !s.noLegend);
        if (!ents.length) return;
        const ctx = this.ctx, b = this.box;
        const rows = ents.length;
        let wmax = 0;
        for (const e of ents) wmax = Math.max(wmax, P.textWidth(ctx, e.name, 9));
        const lw = wmax + 30, lh = rows * 11 + 6;
        // pick the emptiest corner: top right unless the data crowds it
        let x = b.x + b.w - lw - 4, y = b.y + 4;
        if (this.legendLoc === 'topleft') x = b.x + 4;
        ctx.fillStyle = 'rgba(255,255,255,0.85)';
        ctx.fillRect(x, y, lw, lh);
        let yy = y + 8;
        for (const e of ents) {
            ctx.strokeStyle = e.color; ctx.lineWidth = e.width;
            ctx.setLineDash(e.dash || []);
            ctx.beginPath(); ctx.moveTo(x + 4, yy); ctx.lineTo(x + 22, yy); ctx.stroke();
            ctx.setLineDash([]);
            P.text(ctx, e.name, x + 26, yy, { size: 9 });
            yy += 11;
        }
    }

    drawMarks() {
        const ctx = this.ctx, b = this.box;
        for (const m of this.marks) {
            const s = this.series[m.si];
            if (!s) continue;
            const X = this.px(s.x[m.i]), Y = this.py(s.y[m.i], s.right);
            ctx.fillStyle = '#000';
            ctx.beginPath(); ctx.arc(X, Y, 3, 0, 2 * Math.PI); ctx.fill();
            const lines = m.text;
            const w = Math.max(...lines.map(t => P.textWidth(ctx, t, 9))) + 8;
            const h = lines.length * 11 + 6;
            let tx = X + 8, ty = Y - h - 6;
            if (tx + w > b.x + b.w) tx = X - w - 8;
            if (ty < b.y) ty = Y + 8;
            ctx.fillStyle = '#fff'; ctx.strokeStyle = '#b3b3b3'; ctx.lineWidth = 0.8;
            ctx.fillRect(tx, ty, w, h); ctx.strokeRect(tx + 0.5, ty + 0.5, w, h);
            let yy = ty + 9;
            for (const t of lines) { P.text(ctx, t, tx + 4, yy, { size: 9 }); yy += 11; }
        }
    }

    // nearest sample of any series to a pixel position
    hit(mx, my, maxPx) {
        let best = null, bd = Infinity;
        for (let si = 0; si < this.series.length; si++) {
            const s = this.series[si];
            for (let i = 0; i < s.x.length; i++) {
                if (!isFinite(s.y[i])) continue;
                const dx = this.px(s.x[i]) - mx, dy = this.py(s.y[i], s.right) - my;
                const d = dx * dx + dy * dy;
                if (d < bd) { bd = d; best = { si, i }; }
            }
        }
        if (!best || Math.sqrt(bd) > (maxPx || 14)) return null;
        return best;
    }
};

function padLim(lo, hi, pad) {
    if (!isFinite(lo) || !isFinite(hi)) return [0, 1];
    if (hi === lo) { hi = lo + Math.max(1e-12, Math.abs(lo) * 0.1 + 1e-12); }
    const p = (pad === undefined ? 0.1 : pad) * (hi - lo);
    return [lo < 0 ? lo - p : Math.min(lo, 0), hi + p];
}
P.padLim = padLim;

// ------------------------------------------------------------- beam sketch
// Port of PlotBeam2D.m / plotBeam: beam line, supports and dimension lines.
// The 1:1 data aspect comes from the y-limits, so the sketch lines up with the
// mode shape plot below it.  Returns the support x positions in pixels for
// hit-testing.
P.beam = function (canvas, box, L, bc, opt) {
    opt = opt || {};
    const dpr = window.devicePixelRatio || 1;
    if (canvas.width !== Math.round(canvas.clientWidth * dpr)) {
        canvas.width = Math.round(canvas.clientWidth * dpr);
        canvas.height = Math.round(canvas.clientHeight * dpr);
    }
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, canvas.clientWidth, canvas.clientHeight);

    const Ltot = L.reduce((a, b) => a + b, 0);
    const xl = [-0.02 * Ltot, 1.02 * Ltot];
    const ysp = (xl[1] - xl[0]) * box.h / box.w;          // 1:1 scale
    const yl = [-0.043 * Ltot - 0.5 * ysp, -0.043 * Ltot + 0.5 * ysp];
    const X = v => box.x + (v - xl[0]) / (xl[1] - xl[0]) * box.w;
    const Y = v => box.y + box.h - (v - yl[0]) / (yl[1] - yl[0]) * box.h;

    ctx.lineWidth = 4; ctx.strokeStyle = '#999';
    ctx.beginPath(); ctx.moveTo(X(0), Y(0)); ctx.lineTo(X(Ltot), Y(0)); ctx.stroke();

    const pos = Array.from(new Set(bc.map(b => b[0]))).sort((a, b) => a - b);
    ctx.strokeStyle = '#000'; ctx.fillStyle = '#000';
    for (const p of pos) {
        const dirs = bc.filter(b => b[0] === p).map(b => b[1]);
        const hasV = dirs.includes(1), hasR = dirs.includes(2);
        ctx.lineWidth = 1;
        if (hasV && hasR) {                                   // clamped
            ctx.beginPath();
            ctx.moveTo(X(p), Y(-0.03 * Ltot)); ctx.lineTo(X(p), Y(0.03 * Ltot)); ctx.stroke();
        } else if (hasV) {                                    // pinned
            ctx.beginPath(); ctx.arc(X(p), Y(0), 2, 0, 2 * Math.PI); ctx.fill();
            ctx.beginPath();
            ctx.moveTo(X(p - 0.02 * Ltot), Y(-0.03 * Ltot));
            ctx.lineTo(X(p + 0.02 * Ltot), Y(-0.03 * Ltot));
            ctx.lineTo(X(p), Y(0));
            ctx.closePath(); ctx.stroke();
        } else if (hasR) {                                    // rotation only
            ctx.beginPath();
            ctx.moveTo(X(p), Y(-0.03 * Ltot)); ctx.lineTo(X(p), Y(0.03 * Ltot)); ctx.stroke();
            const r = 0.01 * Ltot, s = (p === Ltot) ? 0 : 1;
            for (const yc of [-1.5 * r, 1.5 * r]) {
                ctx.beginPath();
                ctx.ellipse(X(p - 2 * r * s + r), Y(yc), Math.abs(X(r) - X(0)),
                            Math.abs(Y(r) - Y(0)), 0, 0, 2 * Math.PI);
                ctx.stroke();
            }
        }
    }

    // dimension lines
    const dim = Array.from(new Set([0, ...pos, Ltot])).sort((a, b) => a - b);
    ctx.lineWidth = 0.5;
    ctx.beginPath();
    ctx.moveTo(X(0), Y(-0.1 * Ltot)); ctx.lineTo(X(Ltot), Y(-0.1 * Ltot));
    for (const p of dim) {
        ctx.moveTo(X(p), Y(-0.1 * Ltot - 0.8 * 0.02 * Ltot));
        ctx.lineTo(X(p), Y(-0.1 * Ltot + 2 * 0.02 * Ltot));
        ctx.moveTo(X(p - 0.015 * Ltot), Y(-0.1 * Ltot - 0.015 * Ltot));
        ctx.lineTo(X(p + 0.015 * Ltot), Y(-0.1 * Ltot + 0.015 * Ltot));
    }
    ctx.stroke();
    for (let i = 0; i < dim.length - 1; i++) {
        const Li = dim[i + 1] - dim[i];
        P.text(ctx, Li.toFixed(1), X(dim[i] + 0.5 * Li), Y(-0.075 * Ltot), { align: 'center', size: 9 });
    }

    if (opt.marks) {
        for (const m of opt.marks) {
            ctx.strokeStyle = m.color || '#D95319'; ctx.lineWidth = 1;
            ctx.beginPath();
            ctx.moveTo(X(m.x), Y(0.005 * Ltot)); ctx.lineTo(X(m.x), Y(0.035 * Ltot));
            ctx.stroke();
            P.text(ctx, m.label, X(m.x), Y(0.045 * Ltot), { align: 'center', size: 8, color: m.color || '#D95319' });
        }
    }
    return { X, Y, pos, supPix: pos.map(p => X(p)) };
};

if (typeof module !== 'undefined' && module.exports) module.exports = P;
root.Plot = P;
})(typeof self !== 'undefined' ? self : this);

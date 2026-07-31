/* dynamicfactor.js ----------------------------------------------------------
   EN 1991-2 dynamic amplification factor, phi' = y_dyn/y_stat - 1.
   Port of DynamicFactor.m.

   For each response quantity, phi' is formed per train and per speed as the
   peak dynamic response divided by the peak static response of the SAME train,
   minus one.  The two peaks need not occur at the same position or instant.
   The static peak does not depend on speed, so it is broadcast across the
   speed axis.

   Acceleration is omitted: it has no static counterpart.  phi' is formed
   like-for-like - dynamic max over static max (Mmax, RFmax), dynamic min over
   static min (Mmin, RFmin), magnitude over magnitude (dmax, rmax, Vmax,
   RMmax).  Where a static peak is zero - hogging on a single simply supported
   span, so Mmin = 0 - phi' is undefined and comes back as NaN.
---------------------------------------------------------------------------*/
(function (root) {
'use strict';

const FIELDS = ['dmax', 'rmax', 'Mmax', 'Mmin', 'Vmax', 'RFmax', 'RFmin', 'RMmax'];

// stat: {field: [per train scalar]}, dyn: {field: [per train Float64Array(nv)]}
function dynamicFactor(stat, dyn) {
    const phi = {};
    for (const f of FIELDS) {
        const ys = stat[f], yd = dyn[f];
        if (!ys || !yd || !ys.length || !yd.length) { phi[f] = []; continue; }
        if (ys.length !== yd.length)
            throw new Error('dynamicFactor: different train counts for "' + f + '"');
        phi[f] = yd.map((row, k) => {
            const s = ys[k];
            return Float64Array.from(row, v => {
                const r = v / s;
                return isFinite(r) ? r - 1 : NaN;
            });
        });
    }
    return phi;
}

dynamicFactor.FIELDS = FIELDS;

if (typeof module !== 'undefined' && module.exports) module.exports = dynamicFactor;
root.dynamicFactor = dynamicFactor;
})(typeof self !== 'undefined' ? self : this);

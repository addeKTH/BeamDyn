/* trainloads.js -------------------------------------------------------------
   Axle definitions for the two load families the GUI offers.

   HSLM-A        EN 1991-2:2024 figure 8.12 and table 8.4, ported from
                 TrainLoadHSR.m ('HSLMA').  Ten Universal Trains.
   A1 ... E5     EN 15528:2021 table A.1 reference wagons, ported from
                 TrainLoadSOU.m (the EN 15528 branch).  Reference wagon A is
                 called 'A1' so as not to clash with SOU 1938 type A.

   A train is {x, P}: axle position from the leading axle (m) and axle load
   (N), which is what the solvers consume.  Vertical loads only, no load
   distribution here - the solvers apply inp.W themselves.
---------------------------------------------------------------------------*/
(function (root) {
'use strict';

const TL = {};
const G = 9.81;

// EN 1991-2 table 8.4: N coaches, coach length D, bogie axle spacing dBA,
// point force P (kN), for Universal Trains A1 to A10.
TL.HSLMA_TABLE = [
    [18, 18, 2.0, 170], [17, 19, 3.5, 200], [16, 20, 2.0, 180], [15, 21, 3.0, 190],
    [14, 22, 2.0, 170], [13, 23, 2.0, 180], [13, 24, 2.0, 190], [12, 25, 2.5, 190],
    [11, 26, 2.0, 210], [11, 27, 2.0, 210]
];

// EN 15528 table A.1: name, axle load (t), front overhang, internal axle
// spacings, rear overhang (m).
TL.EN15528_TABLE = [
    ['A1', 16.0, 1.5, [1.8, 6.20, 1.8], 1.5],
    ['B1', 18.0, 1.5, [1.8, 7.80, 1.8], 1.5],
    ['B2', 18.0, 1.5, [1.8, 4.65, 1.8], 1.5],
    ['C2', 20.0, 1.5, [1.8, 5.90, 1.8], 1.5],
    ['C3', 20.0, 1.5, [1.8, 4.50, 1.8], 1.5],
    ['C4', 20.0, 1.5, [1.8, 3.40, 1.8], 1.5],
    ['D2', 22.5, 1.5, [1.8, 7.45, 1.8], 1.5],
    ['D3', 22.5, 1.5, [1.8, 5.90, 1.8], 1.5],
    ['D4', 22.5, 1.5, [1.8, 4.65, 1.8], 1.5],
    ['D5', 22.5, 1.5, [1.8, 3.62, 1.8], 1.5],
    ['E4', 25.0, 1.5, [1.8, 5.90, 1.8], 1.5],
    ['E5', 25.0, 1.5, [1.8, 4.75, 1.8], 1.5]
];

TL.CATEGORIES = TL.EN15528_TABLE.map(r => r[0]);

// ------------------------------------------------------------------- HSLM-A
// Power car (4 axles), end coach, N+1 articulation bogies, end coach, mirrored
// power car.  2N + 14 axles in all.
TL.hslmA = function (k) {
    const row = TL.HSLMA_TABLE[k - 1];
    if (!row) throw new Error('hslmA: train must be 1 to 10');
    const N = row[0], D = row[1], dBA = row[2], P = row[3];

    const pc = [0, 3, 14, 17];                     // power car
    const c0 = 17 + 3.525 / 2;                     // coupling, power car to end coach
    const xL = [17 + 3.525, 17 + 3.525 + dBA];     // outer bogie, leading end coach
    const Lt = c0 + (N + 2) * D + 3.525 / 2 + 17;  // total train length
    const xT = [Lt - 17 - 3.525 - dBA, Lt - 17 - 3.525];

    const x = [];
    pc.forEach(v => x.push(v));
    xL.forEach(v => x.push(v));
    for (let i = 1; i <= N + 1; i++) {             // articulation bogies
        const sb = c0 + i * D;
        x.push(sb - dBA / 2);
        x.push(sb + dBA / 2);
    }
    xT.forEach(v => x.push(v));
    for (let i = pc.length - 1; i >= 0; i--) x.push(Lt - pc[i]);

    const x0 = x[0];
    return {
        name: 'HSLM-A' + k,
        x: Float64Array.from(x, v => v - x0),
        P: Float64Array.from(x, () => P * 1e3)
    };
};

// ------------------------------------------------------- EN 15528 category
// nWagon identical reference wagons, coupled end to end.
TL.en15528 = function (name, nWagon) {
    const row = TL.EN15528_TABLE.find(r => r[0] === name);
    if (!row) throw new Error('en15528: unknown line category ' + name);
    const P0 = row[1], ohF = row[2], sp = row[3], ohR = row[4];
    const nw = Math.max(1, Math.round(nWagon || 1));

    const local = [ohF];                           // axle positions in one wagon
    for (let i = 0; i < sp.length; i++) local.push(local[local.length - 1] + sp[i]);
    const Lw = ohF + sp.reduce((a, b) => a + b, 0) + ohR;

    const x = [];
    for (let w = 0; w < nw; w++) for (const v of local) x.push(w * Lw + v);
    const x0 = x[0];
    return {
        name: name,
        x: Float64Array.from(x, v => v - x0),
        P: Float64Array.from(x, () => P0 * G * 1e3)
    };
};

// --------------------------------------------------------------- GUI helper
// The load pull-down maps to a set of trains: 'HSLM-A' is the ten Universal
// Trains, a line category is one train of nCar wagons.
TL.buildTrains = function (loadName, nCar) {
    if (loadName === 'HSLM-A') {
        const F = [];
        for (let k = 1; k <= 10; k++) F.push(TL.hslmA(k));
        return F;
    }
    return [TL.en15528(loadName, nCar)];
};

TL.LOAD_LIST = ['HSLM-A'].concat(TL.CATEGORIES);

if (typeof module !== 'undefined' && module.exports) module.exports = TL;
root.TL = TL;
})(typeof self !== 'undefined' ? self : this);

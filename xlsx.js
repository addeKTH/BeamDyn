/* xlsx.js -------------------------------------------------------------------
   A minimal Excel writer: enough to put several sheets of numbers and text
   into a real .xlsx, with no dependencies and no build step.

   Zip entries are stored uncompressed, strings are written inline (no shared
   string table), and numbers keep full precision.  Excel, LibreOffice and
   Numbers all open the result.

   Usage:  XLSX.write([{name:'Input', rows:[['a',1],['b',2]]}, ...])  -> Blob
---------------------------------------------------------------------------*/
(function (root) {
'use strict';

const XLSX = {};

// ------------------------------------------------------------------- crc32
const CRC = (() => {
    const t = new Uint32Array(256);
    for (let n = 0; n < 256; n++) {
        let c = n;
        for (let k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
        t[n] = c >>> 0;
    }
    return t;
})();

function crc32(buf) {
    let c = 0xFFFFFFFF;
    for (let i = 0; i < buf.length; i++) c = CRC[(c ^ buf[i]) & 0xFF] ^ (c >>> 8);
    return (c ^ 0xFFFFFFFF) >>> 0;
}

const enc = new TextEncoder();

// ---------------------------------------------------------------- zip, store
function zip(files) {
    const chunks = [], central = [];
    let offset = 0;
    for (const f of files) {
        const name = enc.encode(f.name);
        const data = f.data;
        const crc = crc32(data);
        const lh = new Uint8Array(30 + name.length);
        const dv = new DataView(lh.buffer);
        dv.setUint32(0, 0x04034b50, true);
        dv.setUint16(4, 20, true); dv.setUint16(6, 0, true); dv.setUint16(8, 0, true);
        dv.setUint16(10, 0, true); dv.setUint16(12, 0, true);
        dv.setUint32(14, crc, true);
        dv.setUint32(18, data.length, true); dv.setUint32(22, data.length, true);
        dv.setUint16(26, name.length, true); dv.setUint16(28, 0, true);
        lh.set(name, 30);
        chunks.push(lh, data);

        const ch = new Uint8Array(46 + name.length);
        const cv = new DataView(ch.buffer);
        cv.setUint32(0, 0x02014b50, true);
        cv.setUint16(4, 20, true); cv.setUint16(6, 20, true);
        cv.setUint16(8, 0, true); cv.setUint16(10, 0, true);
        cv.setUint16(12, 0, true); cv.setUint16(14, 0, true);
        cv.setUint32(16, crc, true);
        cv.setUint32(20, data.length, true); cv.setUint32(24, data.length, true);
        cv.setUint16(28, name.length, true);
        cv.setUint32(42, offset, true);
        ch.set(name, 46);
        central.push(ch);
        offset += lh.length + data.length;
    }
    let clen = 0; for (const c of central) clen += c.length;
    const end = new Uint8Array(22);
    const ev = new DataView(end.buffer);
    ev.setUint32(0, 0x06054b50, true);
    ev.setUint16(8, files.length, true); ev.setUint16(10, files.length, true);
    ev.setUint32(12, clen, true); ev.setUint32(16, offset, true);
    return new Blob([...chunks, ...central, end], { type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' });
}

// ---------------------------------------------------------------- sheet xml
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;')
                          .replace(/>/g, '&gt;').replace(/"/g, '&quot;');

function colName(n) {                                    // 1 -> A
    let s = '';
    while (n > 0) { const r = (n - 1) % 26; s = String.fromCharCode(65 + r) + s; n = Math.floor((n - 1) / 26); }
    return s;
}

function sheetXml(rows) {
    let x = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>';
    for (let r = 0; r < rows.length; r++) {
        const row = rows[r] || [];
        x += '<row r="' + (r + 1) + '">';
        for (let c = 0; c < row.length; c++) {
            const v = row[c];
            if (v === null || v === undefined || v === '') continue;
            const ref = colName(c + 1) + (r + 1);
            if (typeof v === 'number' && isFinite(v)) {
                x += '<c r="' + ref + '"><v>' + v + '</v></c>';
            } else {
                x += '<c r="' + ref + '" t="inlineStr"><is><t xml:space="preserve">' +
                     esc(v) + '</t></is></c>';
            }
        }
        x += '</row>';
    }
    return x + '</sheetData></worksheet>';
}

// sheets: [{name, rows}]
XLSX.write = function (sheets) {
    const files = [];
    const add = (name, str) => files.push({ name, data: enc.encode(str) });

    let types = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">' +
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>' +
        '<Default Extension="xml" ContentType="application/xml"/>' +
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>';
    let wb = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" ' +
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>';
    let rels = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">';

    sheets.forEach((s, i) => {
        const n = i + 1;
        add('xl/worksheets/sheet' + n + '.xml', sheetXml(s.rows));
        types += '<Override PartName="/xl/worksheets/sheet' + n + '.xml" ' +
                 'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>';
        wb += '<sheet name="' + esc(s.name) + '" sheetId="' + n + '" r:id="rId' + n + '"/>';
        rels += '<Relationship Id="rId' + n + '" ' +
                'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" ' +
                'Target="worksheets/sheet' + n + '.xml"/>';
    });
    types += '</Types>';
    wb += '</sheets></workbook>';
    rels += '</Relationships>';

    add('[Content_Types].xml', types);
    add('_rels/.rels', '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' +
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">' +
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>' +
        '</Relationships>');
    add('xl/workbook.xml', wb);
    add('xl/_rels/workbook.xml.rels', rels);
    return zip(files);
};

// Ask where to put the file.  Chrome and Edge on the desktop have the File
// System Access API, which gives a real Save-as dialog with an editable name;
// everywhere else this falls back to a download, which lands in the download
// folder under the suggested name.  Returns the name written, or null if the
// dialog was cancelled.
XLSX.save = async function (sheets, filename) {
    const blob = XLSX.write(sheets);
    if (typeof window !== 'undefined' && window.showSaveFilePicker) {
        let handle = null;
        try {
            handle = await window.showSaveFilePicker({
                suggestedName: filename,
                types: [{
                    description: 'Excel workbook',
                    accept: { 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet': ['.xlsx'] }
                }]
            });
        } catch (e) {
            if (e && e.name === 'AbortError') return null;      // cancelled
            handle = null;                                       // blocked, fall back
        }
        if (handle) {
            const ws = await handle.createWritable();
            await ws.write(blob);
            await ws.close();
            return handle.name;
        }
    }
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = filename;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 5000);
    return filename;
};

XLSX.download = function (sheets, filename) {
    const blob = XLSX.write(sheets);
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = filename;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 5000);
};

if (typeof module !== 'undefined' && module.exports) module.exports = XLSX;
root.XLSX = XLSX;
})(typeof self !== 'undefined' ? self : this);

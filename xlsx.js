/*--------------------------------------------------------------------------
  Minimal .xlsx writer with scatter charts.  No dependencies.
  Works in the browser and in Node.

  An .xlsx is a ZIP of XML parts.  Everything here is hand-built: a STORE-only
  ZIP (no compression, so no deflate implementation is needed) plus the
  SpreadsheetML and DrawingML parts for a worksheet with an embedded chart.

  writeXlsx(sheets) -> Uint8Array
    sheets: [{ name, headers:[...], cols:[[...], ...],
               chart:{title,xlab,ylab}, series:[colIdx, ...] }]
    cols[0] is the x-axis.  Without `series`, every other column becomes a
    line; with it, only the listed column indices are plotted.

  Andreas Andersson, 2026-07-16
  assisted by Claude Opus 4.8
--------------------------------------------------------------------------*/
(function (root) {
'use strict';

/* ===================== ZIP (store only) ================================ */

const CRC = (() => {
  const t = new Uint32Array(256);
  for (let i = 0; i < 256; i++) {
    let c = i;
    for (let k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
    t[i] = c >>> 0;
  }
  return t;
})();

function crc32(u8) {
  let c = 0xFFFFFFFF;
  for (let i = 0; i < u8.length; i++) c = CRC[(c ^ u8[i]) & 0xFF] ^ (c >>> 8);
  return (c ^ 0xFFFFFFFF) >>> 0;
}

function utf8(s) {
  if (typeof TextEncoder !== 'undefined') return new TextEncoder().encode(s);
  return Uint8Array.from(Buffer.from(s, 'utf8'));           // Node fallback
}

/** files: [{name, data:Uint8Array}] -> Uint8Array of a STORE-only zip */
function zip(files) {
  const parts = [], central = [];
  let off = 0;
  // DOS timestamp: fixed, so the output is reproducible
  const time = 0, date = ((2026 - 1980) << 9) | (7 << 5) | 16;

  for (const f of files) {
    const nm = utf8(f.name), d = f.data, c = crc32(d);
    const lh = new Uint8Array(30 + nm.length);
    const v = new DataView(lh.buffer);
    v.setUint32(0, 0x04034b50, true);
    v.setUint16(4, 20, true); v.setUint16(6, 0, true); v.setUint16(8, 0, true);
    v.setUint16(10, time, true); v.setUint16(12, date, true);
    v.setUint32(14, c, true);
    v.setUint32(18, d.length, true); v.setUint32(22, d.length, true);
    v.setUint16(26, nm.length, true); v.setUint16(28, 0, true);
    lh.set(nm, 30);
    parts.push(lh, d);

    const ch = new Uint8Array(46 + nm.length);
    const w = new DataView(ch.buffer);
    w.setUint32(0, 0x02014b50, true);
    w.setUint16(4, 20, true); w.setUint16(6, 20, true);
    w.setUint16(8, 0, true); w.setUint16(10, 0, true);
    w.setUint16(12, time, true); w.setUint16(14, date, true);
    w.setUint32(16, c, true);
    w.setUint32(20, d.length, true); w.setUint32(24, d.length, true);
    w.setUint16(28, nm.length, true);
    w.setUint16(30, 0, true); w.setUint16(32, 0, true);
    w.setUint16(34, 0, true); w.setUint16(36, 0, true);
    w.setUint32(38, 0, true); w.setUint32(42, off, true);
    ch.set(nm, 46);
    central.push(ch);
    off += lh.length + d.length;
  }

  const cdSize = central.reduce((a, b) => a + b.length, 0);
  const eo = new Uint8Array(22);
  const e = new DataView(eo.buffer);
  e.setUint32(0, 0x06054b50, true);
  e.setUint16(4, 0, true); e.setUint16(6, 0, true);
  e.setUint16(8, files.length, true); e.setUint16(10, files.length, true);
  e.setUint32(12, cdSize, true); e.setUint32(16, off, true);
  e.setUint16(20, 0, true);

  const all = parts.concat(central, [eo]);
  const total = all.reduce((a, b) => a + b.length, 0);
  const out = new Uint8Array(total);
  let p = 0;
  for (const a of all) { out.set(a, p); p += a.length; }
  return out;
}

/* ===================== XML helpers ==================================== */

const XE = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;')
  .replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const HDR = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n';
const R_NS = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships';

function colName(n) {                       // 0 -> A, 26 -> AA
  let s = ''; n++;
  while (n > 0) { const r = (n - 1) % 26; s = String.fromCharCode(65 + r) + s; n = (n - r - 1) / 26; }
  return s;
}
const num = v => (v === null || v === undefined || !isFinite(v)) ? '' : String(v);

/* Okabe-Ito palette, matching the page */
const PAL = ['0072B2', 'D55E00', '009E73', 'CC79A7', '56B4E9', 'E69F00'];

/* ===================== parts ========================================== */

function sheetXml(sh, hasChart) {
  const nr = sh.cols[0].length;
  let x = HDR + '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
        + ` xmlns:r="${R_NS}"><sheetData>`;
  // header row
  x += '<row r="1">';
  sh.headers.forEach((h, j) => {
    x += `<c r="${colName(j)}1" s="1" t="inlineStr"><is><t>${XE(h)}</t></is></c>`;
  });
  x += '</row>';
  // data
  for (let i = 0; i < nr; i++) {
    x += `<row r="${i + 2}">`;
    for (let j = 0; j < sh.cols.length; j++) {
      const v = num(sh.cols[j][i]);
      if (v !== '') x += `<c r="${colName(j)}${i + 2}"><v>${v}</v></c>`;
    }
    x += '</row>';
  }
  x += '</sheetData>';
  if (hasChart) x += '<drawing r:id="rId1"/>';
  x += '</worksheet>';
  return x;
}

function chartXml(sh, si) {
  const nr = sh.cols[0].length;
  // sh.series: column indices to plot (0 = x).  Default: every column but x.
  const ser = sh.series || Array.from({ length: sh.cols.length - 1 }, (_, i) => i + 1);
  const nser = ser.length;
  const q = `'${sh.name.replace(/'/g, "''")}'`;
  const xRef = `${q}!$A$2:$A$${nr + 1}`;
  const ax1 = 100000000 + si * 10, ax2 = ax1 + 1;

  const title = t => '<c:title><c:tx><c:rich><a:bodyPr/><a:lstStyle/><a:p><a:pPr>'
    + '<a:defRPr sz="1000"/></a:pPr><a:r><a:rPr lang="en-US" sz="1000"/>'
    + `<a:t>${XE(t)}</a:t></a:r></a:p></c:rich></c:tx><c:overlay val="0"/></c:title>`;

  let s = HDR
    + '<c:chartSpace xmlns:c="http://schemas.openxmlformats.org/drawingml/2006/chart"'
    + ' xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"'
    + ` xmlns:r="${R_NS}"><c:chart>`
    + title(sh.chart.title)
    + '<c:autoTitleDeleted val="0"/><c:plotArea><c:layout/>'
    + '<c:scatterChart><c:scatterStyle val="lineMarker"/><c:varyColors val="0"/>';

  for (let k = 0; k < nser; k++) {
    const col = colName(ser[k]);
    s += '<c:ser>'
      + `<c:idx val="${k}"/><c:order val="${k}"/>`
      + `<c:tx><c:strRef><c:f>${q}!$${col}$1</c:f></c:strRef></c:tx>`
      + `<c:spPr><a:ln w="19050" cap="rnd"><a:solidFill><a:srgbClr val="${PAL[k % PAL.length]}"/>`
      + '</a:solidFill><a:round/></a:ln></c:spPr>'
      + '<c:marker><c:symbol val="none"/></c:marker>'
      + `<c:xVal><c:numRef><c:f>${xRef}</c:f></c:numRef></c:xVal>`
      + `<c:yVal><c:numRef><c:f>${q}!$${col}$2:$${col}$${nr + 1}</c:f></c:numRef></c:yVal>`
      + '<c:smooth val="0"/></c:ser>';
  }

  const axis = (id, cross, pos, lab) =>
    `<c:valAx><c:axId val="${id}"/><c:scaling><c:orientation val="minMax"/></c:scaling>`
    + `<c:delete val="0"/><c:axPos val="${pos}"/>` + title(lab)
    + '<c:numFmt formatCode="General" sourceLinked="1"/>'
    + '<c:majorTickMark val="out"/><c:minorTickMark val="none"/><c:tickLblPos val="nextTo"/>'
    + `<c:crossAx val="${cross}"/><c:crosses val="autoZero"/></c:valAx>`;

  s += `<c:axId val="${ax1}"/><c:axId val="${ax2}"/></c:scatterChart>`
    + axis(ax1, ax2, 'b', sh.chart.xlab)
    + axis(ax2, ax1, 'l', sh.chart.ylab)
    + '</c:plotArea>'
    + (nser > 1 ? '<c:legend><c:legendPos val="r"/><c:overlay val="0"/></c:legend>' : '')
    + '<c:plotVisOnly val="1"/><c:dispBlanksAs val="gap"/></c:chart></c:chartSpace>';
  return s;
}

function drawingXml(si) {
  const c0 = 1 + 8;                      // place the chart clear of the data
  return HDR
    + '<xdr:wsDr xmlns:xdr="http://schemas.openxmlformats.org/drawingml/2006/spreadsheetDrawing"'
    + ' xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
    + '<xdr:twoCellAnchor>'
    + `<xdr:from><xdr:col>${c0}</xdr:col><xdr:colOff>0</xdr:colOff>`
    + '<xdr:row>1</xdr:row><xdr:rowOff>0</xdr:rowOff></xdr:from>'
    + `<xdr:to><xdr:col>${c0 + 9}</xdr:col><xdr:colOff>0</xdr:colOff>`
    + '<xdr:row>22</xdr:row><xdr:rowOff>0</xdr:rowOff></xdr:to>'
    + '<xdr:graphicFrame macro=""><xdr:nvGraphicFramePr>'
    + `<xdr:cNvPr id="${si + 2}" name="Chart ${si + 1}"/><xdr:cNvGraphicFramePr/>`
    + '</xdr:nvGraphicFramePr>'
    + '<xdr:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/></xdr:xfrm>'
    + '<a:graphic><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/chart">'
    + '<c:chart xmlns:c="http://schemas.openxmlformats.org/drawingml/2006/chart"'
    + ` xmlns:r="${R_NS}" r:id="rId1"/>`
    + '</a:graphicData></a:graphic></xdr:graphicFrame>'
    + '<xdr:clientData/></xdr:twoCellAnchor></xdr:wsDr>';
}

const STYLES = HDR
  + '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
  + '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>'
  + '<font><b/><sz val="11"/><name val="Calibri"/></font></fonts>'
  + '<fills count="2"><fill><patternFill patternType="none"/></fill>'
  + '<fill><patternFill patternType="gray125"/></fill></fills>'
  + '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
  + '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
  + '<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
  + '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/></cellXfs>'
  + '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
  + '</styleSheet>';

/* ===================== assemble ======================================= */

function writeXlsx(sheets) {
  const files = [];
  const N = sheets.length;

  let ct = HDR + '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
    + '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
    + '<Default Extension="xml" ContentType="application/xml"/>'
    + '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
    + '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>';
  for (let i = 0; i < N; i++) {
    ct += `<Override PartName="/xl/worksheets/sheet${i + 1}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>`;
    if (sheets[i].chart) {
      ct += `<Override PartName="/xl/drawings/drawing${i + 1}.xml" ContentType="application/vnd.openxmlformats-officedocument.drawing+xml"/>`;
      ct += `<Override PartName="/xl/charts/chart${i + 1}.xml" ContentType="application/vnd.openxmlformats-officedocument.drawingml.chart+xml"/>`;
    }
  }
  ct += '</Types>';
  files.push({ name: '[Content_Types].xml', data: utf8(ct) });

  files.push({ name: '_rels/.rels', data: utf8(HDR
    + `<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">`
    + `<Relationship Id="rId1" Type="${R_NS}/officeDocument" Target="xl/workbook.xml"/>`
    + '</Relationships>') });

  let wb = HDR + '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    + ` xmlns:r="${R_NS}"><sheets>`;
  sheets.forEach((s, i) => {
    wb += `<sheet name="${XE(s.name)}" sheetId="${i + 1}" r:id="rId${i + 1}"/>`;
  });
  wb += '</sheets></workbook>';
  files.push({ name: 'xl/workbook.xml', data: utf8(wb) });

  let wr = HDR + '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">';
  sheets.forEach((s, i) => {
    wr += `<Relationship Id="rId${i + 1}" Type="${R_NS}/worksheet" Target="worksheets/sheet${i + 1}.xml"/>`;
  });
  wr += `<Relationship Id="rId${N + 1}" Type="${R_NS}/styles" Target="styles.xml"/></Relationships>`;
  files.push({ name: 'xl/_rels/workbook.xml.rels', data: utf8(wr) });
  files.push({ name: 'xl/styles.xml', data: utf8(STYLES) });

  sheets.forEach((sh, i) => {
    files.push({ name: `xl/worksheets/sheet${i + 1}.xml`, data: utf8(sheetXml(sh, !!sh.chart)) });
    if (!sh.chart) return;
    files.push({ name: `xl/worksheets/_rels/sheet${i + 1}.xml.rels`, data: utf8(HDR
      + '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
      + `<Relationship Id="rId1" Type="${R_NS}/drawing" Target="../drawings/drawing${i + 1}.xml"/>`
      + '</Relationships>') });
    files.push({ name: `xl/drawings/drawing${i + 1}.xml`, data: utf8(drawingXml(i)) });
    files.push({ name: `xl/drawings/_rels/drawing${i + 1}.xml.rels`, data: utf8(HDR
      + '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
      + `<Relationship Id="rId1" Type="${R_NS}/chart" Target="../charts/chart${i + 1}.xml"/>`
      + '</Relationships>') });
    files.push({ name: `xl/charts/chart${i + 1}.xml`, data: utf8(chartXml(sh, i)) });
  });

  return zip(files);
}

const API = { writeXlsx, zip, crc32, colName };
if (typeof module !== 'undefined' && module.exports) module.exports = API;
else root.XlsxOut = API;

})(typeof self !== 'undefined' ? self : this);

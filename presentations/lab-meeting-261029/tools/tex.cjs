// Typeset LaTeX to MathML at build time: stdin JSON [[src, display], ...] -> stdout JSON [mathml, ...]
const katex = require('./vendor/katex.min.js');

// Chrome's MathML ignores mathvariant other than "normal", so \mathbf, \boldsymbol, \mathcal,
// \mathbb and \mathsf would render as plain letters. Swap in the Unicode math alphanumerics instead.
const SCRIPT_HOLES = { B: 0x212C, E: 0x2130, F: 0x2131, H: 0x210B, I: 0x2110, L: 0x2112, M: 0x2133, R: 0x211B };
const DS_HOLES = { C: 0x2102, H: 0x210D, N: 0x2115, P: 0x2119, Q: 0x211A, R: 0x211D, Z: 0x2124 };
function mapChar(c, v) {
  const A = c >= 'A' && c <= 'Z', a = c >= 'a' && c <= 'z', d = c >= '0' && c <= '9';
  const L = A ? c.charCodeAt(0) - 65 : a ? 26 + c.charCodeAt(0) - 97 : -1;
  const greek = 'αβγδεζηθικλμνξοπρςστυφχψω'.indexOf(c);
  switch (v) {
    case 'bold':
      if (L >= 0) return String.fromCodePoint(0x1D400 + L);
      if (d) return String.fromCodePoint(0x1D7CE + +c);
      if (greek >= 0) return String.fromCodePoint(0x1D6C2 + greek);
      break;
    case 'bold-italic':
      if (L >= 0) return String.fromCodePoint(0x1D468 + L);
      if (greek >= 0) return String.fromCodePoint(0x1D736 + greek);
      break;
    case 'script':
      if (A) return String.fromCodePoint(SCRIPT_HOLES[c] || 0x1D49C + L);
      break;
    case 'double-struck':
      if (A) return String.fromCodePoint(DS_HOLES[c] || 0x1D538 + L);
      if (d) return String.fromCodePoint(0x1D7D8 + +c);
      break;
    case 'sans-serif':
      if (L >= 0) return String.fromCodePoint(0x1D5A0 + L);
      if (d) return String.fromCodePoint(0x1D7E2 + +c);
      break;
  }
  return null;
}
function unicodeVariants(ml) {
  return ml.replace(/<(mi|mn) mathvariant="(bold|bold-italic|script|double-struck|sans-serif)">([^<]*)<\/\1>/g, (m, tag, v, txt) => {
    const chars = [...txt].map(c => mapChar(c, v));
    if (chars.some(c => c === null)) return m;
    return `<mi${chars.length === 1 ? ' mathvariant="normal"' : ''}>${chars.join('')}</mi>`;
  });
}

let buf = '';
process.stdin.on('data', d => buf += d).on('end', () => {
  // Chrome stretches an accent <mo>^</mo> into a tall diamond; pin it to its natural size
  const out = JSON.parse(buf).map(([src, display]) => unicodeVariants(katex.renderToString(src, { displayMode: display, output: 'mathml', throwOnError: true })).replace(/<mo>\^<\/mo>/g, '<mo stretchy="false">ˆ</mo>'));
  process.stdout.write(JSON.stringify(out));
});

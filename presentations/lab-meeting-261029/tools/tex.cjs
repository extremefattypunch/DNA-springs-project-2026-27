// Typeset LaTeX to MathML at build time: stdin JSON [[src, display], ...] -> stdout JSON [mathml, ...]
const katex = require('./vendor/katex.min.js');
let buf = '';
process.stdin.on('data', d => buf += d).on('end', () => {
  const out = JSON.parse(buf).map(([src, display]) => katex.renderToString(src, { displayMode: display, output: 'mathml', throwOnError: true }));
  process.stdout.write(JSON.stringify(out));
});

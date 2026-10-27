// Drive the deck in headless Chrome and capture frames mid-animation.
//   node tools/play.mjs URL OUTDIR action...
// actions: key:<Key>  wait:<ms>  shot:<name>  eval:<js>
import { spawn } from 'node:child_process';
import { writeFileSync, mkdirSync } from 'node:fs';

const [url, outDir, ...actions] = process.argv.slice(2);
mkdirSync(outDir, { recursive: true });
const port = 9300 + Math.floor(Math.random() * 500);
const chrome = spawn('google-chrome-stable', ['--headless=new', '--no-sandbox', '--disable-gpu', '--hide-scrollbars',
  `--remote-debugging-port=${port}`, `--user-data-dir=${outDir}/.chrome-play`, '--window-size=1920,1080', 'about:blank'],
  { stdio: 'ignore' });
const sleep = ms => new Promise(r => setTimeout(r, ms));
let ws, id = 0; const pending = new Map(); const logs = [];
async function connect() {
  for (let i = 0; i < 50; i++) {
    try {
      const list = await (await fetch(`http://127.0.0.1:${port}/json/list`)).json();
      const pg = list.find(t => t.type === 'page');
      if (pg) return pg.webSocketDebuggerUrl;
    } catch (e) {}
    await sleep(200);
  }
  throw new Error('chrome did not start');
}
function send(method, params = {}) {
  return new Promise((res, rej) => { const i = ++id; pending.set(i, { res, rej }); ws.send(JSON.stringify({ id: i, method, params })); });
}
const KEYS = { ArrowRight: 39, ArrowLeft: 37, ' ': 32, PageDown: 34, PageUp: 33, Home: 36, End: 35 };
try {
  ws = new WebSocket(await connect());
  await new Promise(r => ws.addEventListener('open', r));
  ws.addEventListener('message', ev => {
    const m = JSON.parse(ev.data);
    if (m.id && pending.has(m.id)) { const p = pending.get(m.id); pending.delete(m.id); m.error ? p.rej(new Error(m.error.message)) : p.res(m.result); }
    if (m.method === 'Runtime.consoleAPICalled') logs.push(m.params.type + ': ' + m.params.args.map(a => a.value ?? a.description).join(' '));
    if (m.method === 'Runtime.exceptionThrown') logs.push('EXCEPTION: ' + (m.params.exceptionDetails.exception?.description || m.params.exceptionDetails.text));
  });
  await send('Runtime.enable'); await send('Page.enable');
  await send('Emulation.setDeviceMetricsOverride', { width: 1920, height: 1080, deviceScaleFactor: 1, mobile: false });
  await send('Page.navigate', { url });
  await sleep(2500);
  for (const a of actions) {
    const [k, ...rest] = a.split(':'); const v = rest.join(':');
    if (k === 'wait') await sleep(+v);
    else if (k === 'key') {
      const code = KEYS[v] || v.charCodeAt(0);
      await send('Input.dispatchKeyEvent', { type: 'keyDown', key: v, windowsVirtualKeyCode: code, code: v });
      await send('Input.dispatchKeyEvent', { type: 'keyUp', key: v, windowsVirtualKeyCode: code, code: v });
    } else if (k === 'shot') {
      const r = await send('Page.captureScreenshot', { format: 'png' });
      writeFileSync(`${outDir}/${v}.png`, Buffer.from(r.data, 'base64'));
    } else if (k === 'eval') {
      const r = await send('Runtime.evaluate', { expression: v, returnByValue: true });
      console.log('eval', JSON.stringify(r.result.value));
    }
  }
  if (logs.length) console.log(logs.join('\n'));
} finally {
  try { ws && ws.close(); } catch (e) {}
  chrome.kill('SIGKILL');
}

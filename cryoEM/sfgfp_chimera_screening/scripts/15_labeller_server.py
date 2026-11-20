"""Local labelling server for labeler/index.html (stdlib only; binds 127.0.0.1, so only this node's browser sees it).

  GET  /               the labeller page (wrapped in a doctype/head)
  GET  /manifest.json  the 40-micrograph manifest      ($SFGFP_DATA/label_images/manifest.json)
  GET  /img/<file>.jpg denoised / low-pass images       ($SFGFP_DATA/label_images/)
  GET  /labels         {mic_id: label doc} for every saved micrograph
  PUT  /labels/<mic>   save one micrograph's labels -> selections/labels/<mic>.json (atomic)

Run on the desktop node and open http://127.0.0.1:<port>/ in its Firefox:
  $DNASPRING_PY -I scripts/15_labeller_server.py [port]
Label coordinates are display pixels at 2.84 A/px (x2 = CryoSPARC micrograph pixels at 1.42 A/px).
"""
import json
import os
import re
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

REPO = os.environ.get("SFGFP_REPO", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.environ.get("SFGFP_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112")
IMG = os.path.join(DATA, "label_images")
OUT = os.path.join(REPO, "selections", "labels")
PAGE = os.path.join(REPO, "labeler", "index.html")
HEAD = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover"></head><body>')
MICS = {m["mic_id"] for m in json.load(open(os.path.join(IMG, "manifest.json")))}
LABELS = {"R", "L", "S", "P", "D", "T"}


class H(BaseHTTPRequestHandler):
    def _send(self, code, body, ctype="application/json"):
        if isinstance(body, str):
            body = body.encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        p = self.path.split("?")[0]
        if p in ("/", "/index.html"):
            return self._send(200, HEAD + open(PAGE).read() + "</body></html>", "text/html; charset=utf-8")
        if p == "/manifest.json":
            return self._send(200, open(os.path.join(IMG, "manifest.json")).read())
        m = re.fullmatch(r"/img/(m\d\d_(?:dn|lp)\.jpg)", p)
        if m:
            return self._send(200, open(os.path.join(IMG, m.group(1)), "rb").read(), "image/jpeg")
        if p == "/guide.png":
            return self._send(200, open(os.path.join(REPO, "figures", "qc", "label_guide.png"), "rb").read(), "image/png")
        if p == "/labels":
            docs = {}
            if os.path.isdir(OUT):
                for f in sorted(os.listdir(OUT)):
                    if f.endswith(".json"):
                        docs[f[:-5]] = json.load(open(os.path.join(OUT, f)))
            return self._send(200, json.dumps(docs))
        self._send(404, '{"error":"not found"}')

    def do_PUT(self):
        m = re.fullmatch(r"/labels/(m\d\d)", self.path)
        if not m or m.group(1) not in MICS:
            return self._send(404, '{"error":"unknown micrograph"}')
        try:
            doc = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))))
            pts = [{"x": float(q["x"]), "y": float(q["y"]), "l": q["l"]} for q in doc.get("points", [])]
            assert all(q["l"] in LABELS and 0 <= q["x"] <= 1024 and 0 <= q["y"] <= 1024 for q in pts)
        except Exception:
            return self._send(400, '{"error":"bad label document"}')
        clean = {"mic_id": m.group(1), "points": pts, "done": bool(doc.get("done")), "updated": str(doc.get("updated", ""))}
        os.makedirs(OUT, exist_ok=True)
        tmp = os.path.join(OUT, f".{m.group(1)}.tmp")
        json.dump(clean, open(tmp, "w"), indent=0)
        os.replace(tmp, os.path.join(OUT, f"{m.group(1)}.json"))
        self._send(200, '{"ok":true}')

    def log_message(self, fmt, *args):  # quiet: only saves
        if self.command == "PUT":
            sys.stderr.write("saved %s\n" % self.path)


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8766
    print(f"labeller on http://127.0.0.1:{port}/  ({len(MICS)} micrographs; labels -> {OUT})", flush=True)
    ThreadingHTTPServer(("127.0.0.1", port), H).serve_forever()

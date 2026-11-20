"""Summarise the facility rsync listing (raw_inventory/listing.txt) before/after download.

Reads the `rsync --list-only -r` output and writes
  results/raw_inventory.csv          one row per grid square (movies, bytes, time span)
  results/raw_inventory_summary.json totals by file kind + session-level files
Run:  $DNASPRING_PY -I scripts/00_inventory.py
"""
import csv
import json
import os
import re
import sys
from collections import defaultdict

DATA = os.environ.get("SFGFP_DATA", "/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112")
REPO = os.environ.get("SFGFP_REPO", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LISTING = os.path.join(DATA, "raw_inventory", "listing.txt")

LINE = re.compile(r"^(?P<perm>[-dl]\S+)\s+(?P<size>[\d,]+)\s+(?P<date>\S+)\s+(?P<time>\S+)\s+(?P<path>.+)$")


def kind(path):
    p = path.lower()
    if p.endswith("_fractions.tiff") or p.endswith("_fractions.tif"):
        return "movie_tiff"
    if p.endswith(".eer"):
        return "movie_eer"
    if p.endswith("_fractions.mrc"):
        return "movie_mrc"
    if p.endswith(".gain") or "gain" in os.path.basename(p):
        return "gain"
    if "/data/" in p and p.endswith(".mrc"):
        return "data_micrograph_preview"   # EPU single-frame sum written next to each movie
    if "/data/" in p and p.endswith(".xml"):
        return "data_xml"
    if "/data/" in p and p.endswith(".jpg"):
        return "data_jpg"
    if "foilholes" in p:
        return "foilhole_image"
    if "gridsquare_" in os.path.basename(p):
        return "gridsquare_image"
    return "other"


def main():
    rows = []
    with open(LISTING) as fh:
        for line in fh:
            m = LINE.match(line.rstrip("\n"))
            if not m or m["perm"].startswith("d"):
                continue
            rows.append((m["path"], int(m["size"].replace(",", "")), f'{m["date"]} {m["time"]}'))
    by_kind = defaultdict(lambda: [0, 0])
    squares = defaultdict(lambda: {"movies": 0, "movie_bytes": 0, "first": None, "last": None, "disc": ""})
    session_files = []
    for path, size, ts in rows:
        k = kind(path)
        by_kind[k][0] += 1
        by_kind[k][1] += size
        if "/" not in path:
            session_files.append({"path": path, "bytes": size, "mtime": ts})
        sq = re.search(r"(Images-Disc\d+)/(GridSquare_\d+)", path)
        if sq and k.startswith("movie"):
            s = squares[sq.group(2)]
            s["disc"] = sq.group(1)
            s["movies"] += 1
            s["movie_bytes"] += size
            s["first"] = min(filter(None, [s["first"], ts]))
            s["last"] = max(filter(None, [s["last"], ts]))
    os.makedirs(os.path.join(REPO, "results"), exist_ok=True)
    with open(os.path.join(REPO, "results", "raw_inventory.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["disc", "grid_square", "movies", "movie_GB", "first_movie_mtime", "last_movie_mtime"])
        for name, s in sorted(squares.items(), key=lambda kv: kv[1]["first"] or ""):
            w.writerow([s["disc"], name, s["movies"], round(s["movie_bytes"] / 1e9, 2), s["first"], s["last"]])
    summary = {
        "listing": LISTING,
        "files": len(rows),
        "total_TB": round(sum(r[1] for r in rows) / 1e12, 4),
        "by_kind": {k: {"files": v[0], "GB": round(v[1] / 1e9, 2)} for k, v in sorted(by_kind.items())},
        "grid_squares": len(squares),
        "movies": sum(s["movies"] for s in squares.values()),
        "session_level_files": session_files,
    }
    with open(os.path.join(REPO, "results", "raw_inventory_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    json.dump({k: summary[k] for k in ("files", "total_TB", "grid_squares", "movies", "by_kind")}, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()

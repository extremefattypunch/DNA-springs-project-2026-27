# Where the data lives (and when netscratch deletes it)

| What | Path | Size | Purge (90 d after mtime) |
|---|---|---|---|
| Raw EPU session (movies, gain, xml, square/hole images; **no** per-movie preview .mrc) | `/n/netscratch/hekstra_lab/Lab/ian_poon/sfgfp_chimera_screening_20261112/raw/` | ~1.11 TB | **~2027-02-10** (rsync keeps the facility mtimes, 2026-11-12/13) |
| Facility listing + inventory | `.../sfgfp_chimera_screening_20261112/raw_inventory/listing.txt` | 5 MB | ~2027-02-11 |
| CryoSPARC project | `/n/netscratch/hekstra_lab/Lab/ian_poon/cryosparc_projects/CS-sfgfp-chimera-screening/` | grows | 90 d after each job's last write |
| Exports for figures | `.../sfgfp_chimera_screening_20261112/exports/` | small | 90 d |

Source: `rsync://em00000013@10.245.3.142/em00000013` (password in `~/.config/rsync/em00000013.pass`, mode 600,
never in git). **Only Cannon login nodes can reach it** — compute nodes are firewalled (probe job 60306535) — so the
transfer runs via `scripts/download_login.sh` on holylogin05. All movie-dependent work (motion, RBMC) must finish
by **2027-01-31**. Ask the facility how long they keep the session before relying on a re-download.

"""Connect to the running CryoSPARC master.

Host and port are read from the master's config.sh, which `launchcryosparc` rewrites on
every launch, so this keeps working across restarts without hard-coding anything.
"""

import os
import re

CRYODIR = "/n/holylabs/hekstra_lab/Everyone/ianpoon/cryosparc/cryosparc_master"
SETUP = "/n/home05/ianpoon/cryosparc-setup/configure.sh"


def _scan(path, keys):
    out = {}
    with open(path) as f:
        for line in f:
            m = re.match(r'\s*(?:export\s+)?([A-Z_]+)=["\']?([^"\'\n#]+)', line)
            if m and m.group(1) in keys:
                out[m.group(1)] = m.group(2).strip()
    return out


def master_config():
    cfg = _scan(os.path.join(CRYODIR, "config.sh"),
                {"CRYOSPARC_MASTER_HOSTNAME", "CRYOSPARC_BASE_PORT", "CRYOSPARC_LICENSE_ID"})
    creds = _scan(SETUP, {"EMAIL", "PASSWORD"})
    return {
        "host": cfg["CRYOSPARC_MASTER_HOSTNAME"],
        "port": int(cfg["CRYOSPARC_BASE_PORT"]),
        "license": cfg.get("CRYOSPARC_LICENSE_ID"),
        "email": creds.get("EMAIL", "ian_poon@fas.harvard.edu"),
        "password": creds.get("PASSWORD"),
    }


def connect():
    from cryosparc.tools import CryoSPARC
    c = master_config()
    cs = CryoSPARC(host=c["host"], base_port=c["port"], license=c["license"],
                   email=c["email"], password=c["password"])
    assert cs.test_connection(), f"cannot reach CryoSPARC at {c['host']}:{c['port']}"
    return cs, c


if __name__ == "__main__":
    cs, c = connect()
    print(f"connected to http://{c['host']}:{c['port']} as {c['email']}")
    print("existing projects:", [(p["uid"], p["title"]) for p in cs.cli.list_projects()]
          if hasattr(cs, "cli") else "n/a")

#!/usr/bin/env python3
# Role: Run selection by run length: per-participant fMRIPrep --bids-filter-file JSONs selecting the
#     task-rest run with the most volumes (>= 50).
# Manuscript item: Methods and Supplementary Methods (analysed run).
"""Build per-participant fMRIPrep --bids-filter-file JSONs selecting the task-rest run by run length.

Runs on the cluster. For each participant in the roster it scans the BIDS func/ directory, counts
the volumes of every task-rest BOLD run from the gzipped NIfTI header (dim[4]) and selects the run
with the most volumes (at least MIN_BOLD_VOLS = 50; shorter series are calibration scans).
Selection is never by the dir- entity. The manifest carries an expected-direction column (see
EXPECTED_AP_SUBJECTS below) so that the selection can be checked before a batch submission.

Deterministic tie-break (ties are not expected): 1. more volumes wins; 2. larger file size wins;
3. lexicographically first filename wins.

Outputs (next to this script):
  filters/<subject>_filter.json   fMRIPrep --bids-filter-file, bold query only
  run_selection_manifest.csv      every run seen and what was picked

Exit status: 0 only if every roster participant got a valid, unambiguous filter.

Usage (login node):
  python3 select_runs_by_length.py --roster <roster.txt> [--bids-root <BIDS>] [--session BL]
"""
import argparse
import csv
import gzip
import json
import os
import re
import struct
import sys

GROUP_DIRS = ["Apathy_BL", "Both_BL", "ICB_BL", "Neither_BL"]
MIN_BOLD_VOLS = 50  # config.py MIN_BOLD_VOLS: below this = calibration scan

# Participants whose full run is dir-AP (everyone else PA): supply them as
# EXPECTED_AP_SUBJECTS=sub-...,sub-... . The list feeds only the manifest's direction check
# and the exit status, never the run selection (run length).
EXPECTED_DIR = {s: "AP" for s in os.environ.get("EXPECTED_AP_SUBJECTS", "").split(",") if s}
DEFAULT_EXPECTED_DIR = "PA"

# BIDS entities we mirror from the selected filename into the filter query.
# Maps filename key -> pybids query key.
ENTITY_MAP = {"dir": "direction", "run": "run", "acq": "acquisition",
              "echo": "echo"}


def nvols(path):
    """#volumes (dim[4]) from a gzipped NIfTI-1/2 header, or None."""
    try:
        with gzip.open(path, "rb") as f:
            hdr = f.read(544)
    except OSError:
        return None
    if len(hdr) < 348:
        return None
    sizeof_hdr = struct.unpack("<i", hdr[0:4])[0]
    if sizeof_hdr == 348:  # NIfTI-1 little-endian
        dim = struct.unpack("<8h", hdr[40:56])
    elif sizeof_hdr == 540 and len(hdr) >= 544:  # NIfTI-2 little-endian
        dim = struct.unpack("<8q", hdr[16:80])
    else:
        # try big-endian NIfTI-1
        if struct.unpack(">i", hdr[0:4])[0] == 348:
            dim = struct.unpack(">8h", hdr[40:56])
        else:
            return None
    if not (1 <= dim[0] <= 7):
        return None
    return int(dim[4]) if dim[0] >= 4 else 1


def parse_entities(fname):
    """key -> value for the BIDS entities in a filename."""
    ents = {}
    for key in ENTITY_MAP:
        m = re.search(r"_%s-([A-Za-z0-9]+)" % key, fname)
        if m:
            ents[key] = m.group(1)
    return ents


def matches_filter(fname, query):
    """Would this bold filename match the entity query? (local re-check)"""
    ents = parse_entities(fname)
    for fkey, qkey in ENTITY_MAP.items():
        if qkey in query:
            want = str(query[qkey])
            have = ents.get(fkey)
            if fkey == "run":
                # pybids treats run as an index: run-01 == run-1
                if have is None or int(have) != int(want):
                    return False
            elif have != want:
                return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bids-root",
                    default=os.path.expanduser("~/PPMI/MRI_Data/BIDS"))
    ap.add_argument("--session", default="BL")
    ap.add_argument("--roster", default=None)
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    roster_path = args.roster or os.path.join(here, "roster.txt")
    filt_dir = os.path.join(here, "filters")
    os.makedirs(filt_dir, exist_ok=True)

    subjects = []
    with open(roster_path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                subjects.append(line)
    print("Roster: %d subjects" % len(subjects))

    rows = []
    n_bad = 0
    for subj in subjects:
        row = {"subject": subj, "group": "", "n_bold_runs": 0,
               "runs_found": "", "selected_file": "", "selected_nvols": "",
               "selected_dir": "", "expected_dir": EXPECTED_DIR.get(
                   subj, DEFAULT_EXPECTED_DIR),
               "dir_check": "", "filter_file": "", "status": ""}

        # locate the group dir holding this subject's baseline func
        func_dir = None
        for g in GROUP_DIRS:
            cand = os.path.join(args.bids_root, g, subj,
                                "ses-" + args.session, "func")
            if os.path.isdir(cand):
                func_dir = cand
                row["group"] = g
                break
        if func_dir is None:
            row["status"] = "ERROR: no ses-%s func dir in any group" % \
                args.session
            n_bad += 1
            rows.append(row)
            continue

        bolds = sorted(f for f in os.listdir(func_dir)
                       if "task-rest" in f and f.endswith("_bold.nii.gz"))
        row["n_bold_runs"] = len(bolds)
        if not bolds:
            row["status"] = "ERROR: no task-rest bold runs"
            n_bad += 1
            rows.append(row)
            continue

        scored = []
        for f in bolds:
            full = os.path.join(func_dir, f)
            nv = nvols(full)
            scored.append((f, nv, os.path.getsize(full)))
        row["runs_found"] = ";".join(
            "%s:%s" % (f, "?" if nv is None else nv) for f, nv, _ in scored)

        usable = [s for s in scored if s[1] is not None
                  and s[1] >= MIN_BOLD_VOLS]
        if not usable:
            row["status"] = "ERROR: no run with >= %d volumes" % MIN_BOLD_VOLS
            n_bad += 1
            rows.append(row)
            continue

        # RUN-LENGTH selection with deterministic tie-break
        usable.sort(key=lambda s: (-s[1], -s[2], s[0]))
        sel_file, sel_nv, _ = usable[0]
        row["selected_file"] = sel_file
        row["selected_nvols"] = sel_nv

        ents = parse_entities(sel_file)
        row["selected_dir"] = ents.get("dir", "(none)")
        if "dir" in ents:
            row["dir_check"] = ("MATCH" if ents["dir"] == row["expected_dir"]
                                else "MISMATCH")
        else:
            row["dir_check"] = "NO_DIR_ENTITY"

        # Deliberately NO "session" or "extension" keys here: fMRIPrep's
        # collect_data shares those as layout.get kwargs across ALL
        # modality queries and DELETES any that appear in a filter entry
        # (its "avoid clobbering" block). Including them stripped the
        # extension guard from the t1w query and fed a .json sidecar into
        # the anatomical workflow (pilot crash #2, pilot participant). Session
        # and extension scoping come from fMRIPrep's own defaults.
        query = {"datatype": "func", "suffix": "bold", "task": "rest"}
        for fkey, qkey in ENTITY_MAP.items():
            if fkey in ents:
                # Values must be STRINGS: fMRIPrep's filter loader probes
                # `'Query' in value` on each one and raises TypeError on
                # ints (proven by the pilot crash). pybids casts
                # the raw token ("01") to the run entity's int dtype itself.
                query[qkey] = ents[fkey]

        # the filter must single out exactly the selected run
        hits = [f for f in bolds if matches_filter(f, query)]
        if hits != [sel_file]:
            row["status"] = "ERROR: filter ambiguous, matches %s" % hits
            n_bad += 1
            rows.append(row)
            continue

        fpath = os.path.join(filt_dir, "%s_filter.json" % subj)
        with open(fpath, "w") as fh:
            json.dump({"bold": query}, fh, indent=2, sort_keys=True)
            fh.write("\n")
        row["filter_file"] = fpath
        row["status"] = "OK"
        rows.append(row)

    manifest = os.path.join(here, "run_selection_manifest.csv")
    with open(manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()),
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    n_ok = sum(1 for r in rows if r["status"] == "OK")
    n_mism = sum(1 for r in rows if r["dir_check"] == "MISMATCH")
    print("Manifest: %s" % manifest)
    print("OK: %d / %d   dir MISMATCH flags: %d" % (n_ok, len(rows), n_mism))
    for r in rows:
        if r["status"] != "OK" or r["dir_check"] == "MISMATCH":
            print("  ATTN %s: status=%s dir_check=%s" %
                  (r["subject"], r["status"], r["dir_check"]))
    sys.exit(0 if (n_ok == len(rows) and n_mism == 0) else 1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# Role: NBS surfaces + chord diagrams (atlas-derived labels read with --atlas-labels, printed-width
#     layout); imports plot_nbs_circos_revised.py from this folder.
# Manuscript item: Figure 3 (N = 305 run) and Supplementary Figure S4 (N = 343 run).
"""plot_nbs_combined_revised.py: manuscript Figure 3 (and Supplementary Figure S4 for the extended sample).

One row per significant NBS contrast, five columns: LH lateral, LH medial, RH lateral, RH medial
surface maps of signed node strength, and a chord diagram of the component edges. The figure is
authored at its printed width (--width-mm, default 160 mm), so font sizes are true printed points
(all text at least 6 pt). Parcel labels, ring segments and colours come from the 4S256Parcels atlas
order (--atlas-labels, required); nilearn's Schaefer-200 labels place the surface maps. Writes
<out-name>.png/.pdf and <out-name>_manifest.json (components, drawn nodes and labels, printed
strings, colour range, text checks). Node strength is the sum, over the significant components a
parcel belongs to, of -log10(p_FWE) * sqrt(n_edges), signed by the direction of the component.

Usage:
    python3 plot_nbs_combined_revised.py --nbs-dir <run>/fc_connectome_continuous \
        --atlas-labels atlas_4S256Parcels_labels_xcpd.tsv --fig-dir <out> \
        --contrasts main_apathy diff_apathy_gt_quip interaction_apathyxquip --out-name Figure3
    (--sample-label "extended sample" for the N = 343 run; --svg also writes an SVG)
"""
from __future__ import annotations

import sys
sys.dont_write_bytecode = True   # FIX leave no __pycache__ behind

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np

# circos helpers (atlas-derived labels and layout; draw_circos refuses an edge list that does not
# match the component's n_edges)
from plot_nbs_circos_revised import (
    CONTRAST_TITLES,
    NEG_EDGE_COLOR,
    NET_SHORT,
    POS_EDGE_COLOR,
    SKIP_CONTRASTS,
    _build_parcel_layout,
    audit_texts,
    check_labels_against_decoded,
    draw_circos,
    get_labels,
    load_significant_components,
    network_of,
    segment_at_angle,
)

# Pre-fix canvas (col ratios 1,1,1,1,1.8 x 3.1 in; rows x 3.8 in + 1 in at
# 300 dpi). The PNG keeps exactly these pixel dimensions.
ORIG_W_IN = 5.8 * 3.1
ORIG_ROW_H_IN, ORIG_EXTRA_H_IN, ORIG_DPI = 3.8, 1.0, 300

# Printed font sizes (pt) at --width-mm
FS = dict(title=7.5, row=7.5, header=6.5, summary=6.5, ring=6.0,
          cbar_tick=6.0, cbar_label=6.5, legend=6.5)
# Layout (mm at the printed size)
MM = dict(row_label=6.0, surf_pitch=25.0, surf_ax_w=27.0, title=4.4,
          header=3.2, bottom=8.2, chord_gap=1.5, right=0.8, ring_r=9.5,
          cbar_w=56.0, cbar_h=1.3, cbar_y=5.6)
MINUS = "−"


# ── Copied verbatim from plot_nbs_surface.py ──
def compute_node_strength(nbs_json: dict, p_thresh: float = 0.05,
                          n_parcels: int = 256) -> np.ndarray:  # noqa: C901
    """Compute signed node strength for each parcel.

    For each significant component:
      - Positive component: each node gets +weight
      - Negative component: each node gets -weight
      where weight = -log10(p) * n_edges (captures both significance and size)

    Returns (n_parcels,) array with signed strength values.
    """
    strength = np.zeros(n_parcels)

    for sign_label, key in [("pos", "obs_components_pos"),
                            ("neg", "obs_components_neg")]:
        for comp in nbs_json.get(key, []):
            if comp["p"] >= p_thresh:
                continue
            weight = -np.log10(comp["p"]) * np.sqrt(comp["n_edges"])
            sign = 1.0 if sign_label == "pos" else -1.0
            for node in comp["nodes"]:
                if node < n_parcels:
                    strength[node] += sign * weight

    return strength


def make_cortical_volume(strength_cortical: np.ndarray,
                          schaefer_img) -> "nibabel.Nifti1Image":
    """Map cortical parcel strengths onto the Schaefer atlas volume."""
    import nibabel as nib

    atlas_data = np.asarray(schaefer_img.dataobj, dtype=np.int32)
    out_data = np.zeros_like(atlas_data, dtype=np.float64)

    for parcel_idx in range(len(strength_cortical)):
        val = strength_cortical[parcel_idx]
        if val != 0:
            mask = atlas_data == (parcel_idx + 1)  # parcels are 1-indexed
            out_data[mask] = val

    return nib.Nifti1Image(out_data, schaefer_img.affine, schaefer_img.header)
# ── end of verbatim copy (note: the weight is -log10(p) * sqrt(n_edges);
#    the docstring's "* n_edges" is the original's wording) ──


def plot_surface_views(fig, rects, strength, schaefer_img, fsaverage, vmax):
    """Render the 4 cortical surface views into the given axes rectangles."""
    from nilearn import plotting, surface

    cortical_strength = strength[:200]
    stat_vol = make_cortical_volume(cortical_strength, schaefer_img)

    views = [
        ("left", "lateral"),
        ("left", "medial"),
        ("right", "lateral"),
        ("right", "medial"),
    ]
    axes = []
    for rect, (hemi, view) in zip(rects, views):
        ax = fig.add_axes(rect, projection="3d")
        mesh = fsaverage[f"pial_{hemi}"]
        bg_map = fsaverage[f"sulc_{hemi}"]

        try:
            texture = surface.vol_to_surf(
                stat_vol, mesh, interpolation="nearest", radius=3.0)
        except Exception:
            texture = surface.vol_to_surf(stat_vol, mesh)

        plotting.plot_surf_stat_map(
            mesh, texture,
            hemi=hemi, view=view,
            bg_map=bg_map,
            threshold=0.01,
            cmap="RdBu_r",
            vmax=vmax,
            symmetric_cbar=True,
            colorbar=False,
            axes=ax,
        )
        # FIX transparent 3D-axes background, so the (slightly
        # overlapping) axes never hide a neighbouring brain.
        ax.set_facecolor((1, 1, 1, 0))
        axes.append(ax)
    return axes


def plot_circos_panel(fig, cell, fig_mm, angles, net_arcs, net_of_idx,
                      labels, comps, lines, ring_r_mm):
    """Chord cell: component summary lines on top (centred, as before), ring
    (1 data unit = ring_r_mm) below, centred in the cell.
    cell = (x0, y0, w, h) in mm."""
    x0, y0, w, h = cell
    FW, FH = fig_mm
    R = ring_r_mm
    ax = fig.add_axes([x0 / FW, y0 / FH, w / FW, h / FH])
    ax.set_xlim(-w / 2 / R, w / 2 / R)          # equal mm scale on x and y
    ax.set_ylim(-h / 2 / R, h / 2 / R)
    info = draw_circos(ax, angles, net_arcs, net_of_idx, labels, comps, "",
                       fontsize=FS["ring"])

    renderer = fig.canvas.get_renderer()
    px_per_mm = fig.dpi / 25.4
    stexts = [ax.text(0, 0, txt, fontsize=FS["summary"], fontweight="bold",
                      ha="center", va="center", zorder=6) for txt in lines]
    heights = [t.get_window_extent(renderer).height / px_per_mm for t in stexts]
    widths = [t.get_window_extent(renderer).width / px_per_mm for t in stexts]
    pitch = max(heights) * 1.12
    blk_h = pitch * len(stexts)

    bx0, by0, bx1, by1 = info["block_bbox"]          # data units, ring centre = 0
    top_pad, gap, bot_pad = 0.3, 0.8, 0.3
    band_lo, band_hi = bot_pad, h - top_pad - blk_h - gap
    need = (by1 - by0) * R
    if need > band_hi - band_lo + 1e-6:
        raise RuntimeError(f"chord block ({need:.1f} mm) does not fit under "
                           f"the summary ({band_hi - band_lo:.1f} mm)")
    cx = w / 2 - (bx0 + bx1) / 2 * R                 # ring centre in the cell
    cy = (band_lo + band_hi) / 2 - (by0 + by1) / 2 * R
    if cx + bx0 * R < -1e-6 or cx + bx1 * R > w + 1e-6:
        raise RuntimeError("chord block wider than its cell")
    if max(widths) > w:
        raise RuntimeError("component summary wider than its cell")
    ax.set_xlim(-cx / R, (w - cx) / R)
    ax.set_ylim(-cy / R, (h - cy) / R)
    for k, t in enumerate(stexts):
        ymm = h - top_pad - pitch * (k + 0.5)
        t.set_position(((w / 2 - cx) / R, (ymm - cy) / R))
    info["ring_centre_mm"] = (round(x0 + cx, 3), round(y0 + cy, 3))
    info["ring_radius_mm"] = R
    return ax, info, stexts


def _md5(path: Path) -> str:
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--nbs-dir", required=True, type=Path)
    p.add_argument("--fig-dir", required=True, type=Path)
    p.add_argument("--atlas-labels", required=True,
                   help="4S256Parcels labels: index/label table "
                        "(atlas_4S256Parcels_labels_xcpd.tsv), XCP-D "
                        "timeseries/relmat TSV, or a directory of them")
    p.add_argument("--p-thresh", type=float, default=0.05)
    p.add_argument("--contrasts", nargs="*", default=None,
                   help="Ordered contrast names (default: auto-discover, "
                        "skip redundant mirrors)")
    p.add_argument("--out-name", default="fig_main_combined")
    p.add_argument("--vmax", type=float, default=0,
                   help="Override surface colourbar range (default: auto)")
    p.add_argument("--title", default="",
                   help="Override figure suptitle (default: auto from N/nperm)")
    p.add_argument("--sample-label", default="",
                   help='Inserted in the automatic title, e.g. "extended sample"')
    p.add_argument("--width-mm", type=float, default=160.0,
                   help="Printed width the figure is authored at (default 160)")
    p.add_argument("--svg", action="store_true", help="also write an SVG")
    args = p.parse_args()

    nbs_dir: Path = args.nbs_dir
    out_dir: Path = args.fig_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"],
                         "pdf.fonttype": 42, "ps.fonttype": 42,
                         "axes.unicode_minus": True})

    # ── Atlas + surface ──
    from nilearn import datasets
    import nibabel as nib
    print("Fetching Schaefer 200 atlas...")
    schaefer = datasets.fetch_atlas_schaefer_2018(
        n_rois=200, resolution_mm=2, verbose=0)
    schaefer_img = nib.load(schaefer.maps)
    print("Fetching fsaverage surface...")
    fsaverage = datasets.fetch_surf_fsaverage(mesh="fsaverage5")

    # ── FIX circos layout from the atlas labels ──
    label_src = Path(args.atlas_labels)
    labels = get_labels(label_src)
    n_dec = check_labels_against_decoded(nbs_dir, labels)
    sch = [lb.decode() if isinstance(lb, bytes) else lb for lb in schaefer.labels]
    sch = [lb.replace("7Networks_", "") for lb in sch]
    if sch and "background" in sch[0].lower():
        sch = sch[1:]
    if sch != labels[:200]:
        raise ValueError("nilearn Schaefer-200 labels (surface maps) differ "
                         "from atlas parcels 0-199")
    print(f"Atlas labels: {len(labels)} parcels ({label_src.name}); "
          f"{n_dec} decoder index:label pairs agree; Schaefer 0-199 agree")
    angles, net_arcs, net_of_idx = _build_parcel_layout(labels)

    # ── Discover NBS outputs ──
    all_contrast_files = {}
    for jp in sorted(nbs_dir.glob("nbs_*.json")):
        data_peek = json.load(open(jp))
        cname = data_peek.get("contrast")
        if cname is not None:
            all_contrast_files[cname] = jp

    if args.contrasts:
        contrast_files = []
        for c in args.contrasts:
            if c in all_contrast_files:
                contrast_files.append((c, all_contrast_files[c]))
            else:
                print(f"  WARNING: contrast '{c}' not found, skipping")
    else:
        contrast_files = [(c, p) for c, p in all_contrast_files.items()
                          if c not in SKIP_CONTRASTS]

    # ── Per-contrast: node strength (for surface) + sig components (for circos) ──
    strengths, components = {}, {}
    global_max = 0.0
    keep_order = []
    for cname, cpath in contrast_files:
        _, comps, data = load_significant_components(cpath, args.p_thresh)
        s = compute_node_strength(data, args.p_thresh, n_parcels=256)
        if not comps:
            print(f"  {cname}: no significant components at p<{args.p_thresh}, skipping")
            continue
        strengths[cname] = s
        components[cname] = comps
        global_max = max(global_max, np.abs(s).max())
        keep_order.append((cname, cpath))
        n_pos = sum(1 for c in comps if c[0] == "pos")
        n_neg = sum(1 for c in comps if c[0] == "neg")
        print(f"  {cname}: {n_pos} pos + {n_neg} neg sig components")

    if not keep_order:
        print("ERROR: no significant components at the chosen threshold")
        sys.exit(1)

    vmax = args.vmax if args.vmax > 0 else global_max * 1.1

    # ── Figure construction at the printed size ──
    n_rows = len(keep_order)
    W_px = int(round(ORIG_W_IN * ORIG_DPI))
    H_px = int(round((ORIG_ROW_H_IN * n_rows + ORIG_EXTRA_H_IN) * ORIG_DPI))
    dpi = W_px / (args.width_mm / 25.4)
    fig = plt.figure(figsize=((W_px + 0.3) / dpi, (H_px + 0.3) / dpi), dpi=dpi)
    FW, FH = [v * 25.4 for v in fig.get_size_inches()]

    def rect(x, y, w, h):
        return [x / FW, y / FH, w / FW, h / FH]

    row_h = (FH - MM["title"] - MM["header"] - MM["bottom"]) / n_rows
    chord_x0 = MM["row_label"] + 4 * MM["surf_pitch"] + MM["chord_gap"]
    chord_w = FW - MM["right"] - chord_x0

    texts_to_audit = []
    manifest_rows = []
    for r, (cname, cpath) in enumerate(keep_order):
        title = CONTRAST_TITLES.get(cname, cname)
        comps = components[cname]
        y0 = MM["bottom"] + (n_rows - 1 - r) * row_h

        # Surface: 4 views
        rects = [rect(MM["row_label"] + (c + 0.5) * MM["surf_pitch"]
                      - MM["surf_ax_w"] / 2, y0, MM["surf_ax_w"], row_h)
                 for c in range(4)]
        plot_surface_views(fig, rects, strengths[cname], schaefer_img,
                           fsaverage, vmax)

        # Row label
        fig.text(MM["row_label"] / 2 / FW, (y0 + row_h / 2) / FH, title,
                 rotation=90, ha="center", va="center",
                 fontsize=FS["row"], fontweight="bold")

        # Chord cell (the first row may use the column-header band).
        # Component summaries: same strings as before the fix.
        lines = [f"{'+' if sign == 'pos' else MINUS}{ne}e / {len(nodes)}n "
                 f"(p={pv:.3f})" for sign, nodes, edges, ne, pv in comps]
        cell_h = row_h + (MM["header"] if r == 0 else 0.0)
        ax, info, stexts = plot_circos_panel(
            fig, (chord_x0, y0, chord_w, cell_h), (FW, FH), angles, net_arcs,
            net_of_idx, labels, comps, lines, MM["ring_r"])
        texts_to_audit += info["label_texts"] + stexts

        comp_rows = []
        for ci, ((sign, nodes, edges, ne, pv), txt) in enumerate(zip(comps, lines)):
            drawn = [(a, b) for k, a, b in info["edges_drawn"] if k == ci]
            comp_rows.append({
                "sign": sign, "p": pv, "summary_text": txt,
                "n_edges_json": ne, "n_edges_listed": len(edges),
                "n_edges_drawn": len(drawn),
                "edges_drawn_equal_json": sorted(map(tuple, map(sorted, drawn)))
                == sorted(map(tuple, map(sorted, edges))),
                "n_nodes_json": len(nodes),
                "n_nodes_drawn": len(set(nodes) & set(info["node_xy"])),
                "nodes": [{"index": int(i), "label": labels[i],
                           "segment": NET_SHORT[network_of(labels[i])],
                           "drawn_in_segment": NET_SHORT[info["node_segment"][i]],
                           "angle_deg": round(float(np.degrees(angles[i])) % 360, 3)}
                          for i in nodes],
            })
        manifest_rows.append({
            "contrast": cname, "row_label": title, "json": cpath.name,
            "components": comp_rows,
            "ring_centre_mm": info["ring_centre_mm"],
            "ring_radius_mm": info["ring_radius_mm"],
            "node_strength_nonzero": {str(i): round(float(v), 4) for i, v in
                                      enumerate(strengths[cname]) if v != 0}})

    # ── Column headers above the top row only ──
    for c, lab in enumerate(["LH lateral", "LH medial", "RH lateral", "RH medial"]):
        fig.text((MM["row_label"] + (c + 0.5) * MM["surf_pitch"]) / FW,
                 (FH - MM["title"] - MM["header"] / 2) / FH, lab,
                 ha="center", va="center", fontsize=FS["header"],
                 fontweight="bold")

    # ── Surface colourbar (under the surfaces, fully inside the canvas) ──
    sm = ScalarMappable(cmap="RdBu_r", norm=Normalize(-vmax, vmax))
    sm.set_array([])
    surf_mid = MM["row_label"] + 2 * MM["surf_pitch"]
    cbar_ax = fig.add_axes(rect(surf_mid - MM["cbar_w"] / 2, MM["cbar_y"],
                                MM["cbar_w"], MM["cbar_h"]))
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("Signed node strength (red, hyper-connectivity; "
                   "blue, hypo-connectivity)", fontsize=FS["cbar_label"],
                   labelpad=1.0)
    cbar.ax.tick_params(labelsize=FS["cbar_tick"], length=1.6, width=0.4,
                        pad=1.0)
    cbar.outline.set_linewidth(0.4)

    # ── Circos legend (bottom-right) ──
    legend_elements = [
        Line2D([0], [0], color=POS_EDGE_COLOR, linewidth=1.2,
               label="Positive edge (hyper-connectivity)"),
        Line2D([0], [0], color=NEG_EDGE_COLOR, linewidth=1.2,
               label="Negative edge (hypo-connectivity)"),
    ]
    leg = fig.legend(handles=legend_elements, loc="lower right", ncol=1,
                     fontsize=FS["legend"], frameon=False,
                     bbox_to_anchor=(1 - MM["right"] / FW, 0.5 / FH),
                     handlelength=1.8, handletextpad=0.6, labelspacing=0.25,
                     borderaxespad=0.0, borderpad=0.0)

    # ── Suptitle ──
    subj_file = nbs_dir / "design_subjects.tsv"
    if subj_file.is_file():
        n_subj = sum(1 for _ in open(subj_file)) - 1
    else:
        n_subj = "?"
    first_data = json.load(open(keep_order[0][1]))
    nperm = first_data.get("nperm", "?")
    t_thresh = first_data.get("t_thresh", 3.1)
    if args.title:
        suptitle = args.title
    else:
        lead = "NBS functional connectivity"
        if args.sample_label:
            lead += f", {args.sample_label}"
        nperm_s = f"{nperm:,}" if isinstance(nperm, int) else str(nperm)
        suptitle = (f"{lead} (4S256Parcels, N = {n_subj}, {nperm_s} "
                    f"permutations, |t| > {t_thresh})")
    fig.text(0.5, (FH - MM["title"] / 2) / FH, suptitle, ha="center",
             va="center", fontsize=FS["title"])

    # ── Printed-size text audit (error on any problem) ──
    fig.canvas.draw()
    # only the ticks inside the colour range are drawn
    tick_texts = [t for t, loc in zip(cbar.ax.xaxis.get_majorticklabels(),
                                      cbar.ax.xaxis.get_majorticklocs())
                  if t.get_visible() and t.get_text()
                  and -vmax - 1e-9 <= loc <= vmax + 1e-9]
    texts_to_audit += list(fig.texts) + [cbar.ax.xaxis.label] + tick_texts + \
        list(leg.get_texts())
    problems = audit_texts(fig, texts_to_audit, min_pt=6.0)
    if problems:
        raise RuntimeError("text audit failed:\n  " + "\n  ".join(problems))
    sizes = sorted({round(float(t.get_fontsize()), 2) for t in texts_to_audit
                    if t.get_text().strip()})

    png = out_dir / f"{args.out_name}.png"
    pdf = out_dir / f"{args.out_name}.pdf"
    fig.savefig(png, dpi=dpi, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    if args.svg:
        fig.savefig(out_dir / f"{args.out_name}.svg", facecolor="white")
    plt.close(fig)

    from PIL import Image
    with Image.open(png) as im:
        png_size = im.size
    if png_size != (W_px, H_px):
        raise RuntimeError(f"PNG is {png_size}, expected {(W_px, H_px)}")

    manifest = {
        "script": Path(__file__).name,
        "layout": "atlas-derived chord labels, print-size layout",
        "nbs_dir": nbs_dir.name,
        "atlas_labels": {"source": label_src.name, "md5": _md5(label_src)
                         if label_src.is_file() else None,
                         "decoder_pairs_checked": n_dec},
        "n_subjects": n_subj, "nperm": nperm, "t_thresh": t_thresh,
        "title": suptitle, "vmax": float(vmax),
        "colourbar_ticks": [t.get_text() for t in tick_texts],
        "colourbar_label": cbar.ax.xaxis.label.get_text(),
        "legend": [t.get_text() for t in leg.get_texts()],
        "png_px": list(png_size), "width_mm": args.width_mm,
        "height_mm": round(args.width_mm * H_px / W_px, 3),
        "dpi": round(dpi, 3), "font_sizes_pt": sizes,
        "min_font_pt": min(sizes), "text_audit_problems": problems,
        "rows": manifest_rows,
        "segments": [{"segment": NET_SHORT[net],
                      "n_parcels": sum(1 for x in net_of_idx if x == net),
                      "arc_deg": round(float(np.degrees(e - s)), 3),
                      "start_deg": round(float(np.degrees(s)), 3)}
                     for net, s, e in net_arcs],
    }
    man = out_dir / f"{args.out_name}_manifest.json"
    man.write_text(json.dumps(manifest, indent=1, ensure_ascii=False))
    print(f"\nSaved: {png}  ({png_size[0]} x {png_size[1]} px, "
          f"{args.width_mm:g} mm wide, min font {min(sizes)} pt)")
    print(f"Saved: {pdf}")
    print(f"Saved: {man}")


if __name__ == "__main__":
    main()

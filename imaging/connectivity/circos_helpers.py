#!/usr/bin/env python3
# Role: Circos helpers with atlas-derived parcel labels (atlas_4S256Parcels_labels.tsv, given
#     with --atlas-labels); no hard-coded label list.
# Manuscript item: Figure 3 (chord panels).
"""circos_helpers.py: chord-diagram helpers for the NBS figures.

All 256 parcels of the 4S256Parcels atlas are arranged around the ring by network (Schaefer-200
cortex) or structure (CIT168 nuclei, thalamic nuclei, hippocampus, amygdala, cerebellum) in atlas
order, read from the label table given with --atlas-labels (an index/label table extracted from
XCP-D output headers, or any 4S256Parcels timeseries or relmat TSV). Significant NBS edges are
drawn as bezier curves (red, positive; blue, negative). draw_circos checks that the drawn edge list
matches the component's recorded n_edges, that labels do not overlap or touch the ring and that
leader lines do not cross, and returns what was drawn; audit_texts checks printed font sizes.

Imported by figure3_connectivity_components.py, which draws the manuscript figure.
"""
from __future__ import annotations

import itertools
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
import matplotlib.patches as mpatches
from matplotlib.colors import to_rgb, to_hex
import numpy as np

N_PARCELS = 256
N_CORTICAL = 200

# ── Network / structure definitions (ring segments) ────────────────────
NETWORK_ORDER = [
    "Visual", "Somatomotor", "DorsalAttention", "Salience_VentAttn",
    "Limbic", "Control_FPN", "Default",
    "Hippocampus", "Amygdala", "Putamen", "Caudate",
    "NucleusAccumbens", "GlobusPallidus", "Thalamus",
    "ExtendedAmygdala", "SN_VTA", "Cerebellum", "Other",
]
CORTICAL_NETWORKS = {
    "Visual", "Somatomotor", "DorsalAttention", "Salience_VentAttn",
    "Limbic", "Control_FPN", "Default",
}

NET_COLORS = {
    "Visual":             "#781286",
    "Somatomotor":        "#4682B4",
    "DorsalAttention":    "#00760E",
    "Salience_VentAttn":  "#C43AFA",
    "Limbic":             "#C2B280",
    "Control_FPN":        "#E69422",
    "Default":            "#CD3E4E",
    "Hippocampus":        "#9370DB",
    "Amygdala":           "#FF69B4",
    "Putamen":            "#1E90FF",
    "Caudate":            "#00CED1",
    "NucleusAccumbens":   "#FFD700",
    "GlobusPallidus":     "#8FBC8F",
    "Thalamus":           "#708090",
    "ExtendedAmygdala":   "#FF6347",
    "SN_VTA":             "#2F4F4F",
    "Cerebellum":         "#A0522D",
    "Other":              "#AAAAAA",
}

NET_SHORT = {
    "Visual": "Vis", "Somatomotor": "SomMot",
    "DorsalAttention": "DorsAttn", "Salience_VentAttn": "Sal/VAttn",
    "Limbic": "Limbic", "Control_FPN": "Control",
    "Default": "Default", "Hippocampus": "Hipp", "Amygdala": "Amyg",
    "Putamen": "Put", "Caudate": "Caud",
    "NucleusAccumbens": "NAcc", "GlobusPallidus": "GP",
    "Thalamus": "Thal", "ExtendedAmygdala": "ExtAmyg",
    "SN_VTA": "SN/VTA", "Cerebellum": "Cereb", "Other": "Other",
}

POS_EDGE_COLOR = "#D32F2F"
NEG_EDGE_COLOR = "#1565C0"

CONTRAST_TITLES = {
    "main_apathy": "Apathy",
    "main_quip": "ICB (QUIP)",
    "interaction_apathyxquip": "Apathy × ICB",
    "diff_apathy_gt_quip": "Apathy > ICB",
    "shared_pos": "Shared Effect",
}
# Skip redundant contrasts (mirrors of diff/shared)
SKIP_CONTRASTS = {"diff_quip_gt_apathy", "shared_neg"}

# ── 4S256Parcels label classification ─────────────────
# Schaefer-200 7-network codes (cortex, indices 0-199: LH_<code>_..., RH_...)
_CORTICAL_CODE = {
    "Vis": "Visual", "SomMot": "Somatomotor", "DorsAttn": "DorsalAttention",
    "SalVentAttn": "Salience_VentAttn", "Limbic": "Limbic",
    "Cont": "Control_FPN", "Default": "Default",
}
# CIT168 structure codes (LH-<code> / RH-<code>, indices 200-227)
_CIT168_SEGMENT = {
    "Pu": "Putamen", "Ca": "Caudate", "NAC": "NucleusAccumbens",
    "EXA": "ExtendedAmygdala", "GPe": "GlobusPallidus",
    "GPi": "GlobusPallidus", "SNc_PBP_VTA": "SN_VTA",
    "RN": "Other", "SNr": "Other", "VeP": "Other", "HN": "Other",
    "HTH": "Other", "MN": "Other", "STH": "Other",
}
# Thalamic nuclei (LH-<name> / RH-<name>, indices 228-241)
_THALAMIC_NUCLEI = {
    "Pulvinar", "Anterior", "Medio_Dorsal", "Ventral_Latero_Dorsal",
    "Central_Lateral-Lateral_Posterior-Medial_Pulvinar",
    "Ventral_Anterior", "Ventral_Latero_Ventral",
}


def network_of(label: str) -> str:
    """Map a 4S256Parcels label (XCP-D spelling) to its ring segment.

    Explicit rules only; an unrecognised label raises ValueError so that a
    different atlas can never be drawn silently under wrong segments.
    """
    m = re.match(r"^(LH|RH)_([A-Za-z]+)_", label)
    if m and m.group(2) in _CORTICAL_CODE:
        return _CORTICAL_CODE[m.group(2)]
    m = re.match(r"^(LH|RH)[-_](.+)$", label)
    if m:
        core = m.group(2)
        if core in _CIT168_SEGMENT:
            return _CIT168_SEGMENT[core]
        if core in _THALAMIC_NUCLEI:
            return "Thalamus"
        if core == "Hippocampus":
            return "Hippocampus"
        if core == "Amygdala":
            return "Amygdala"
    if re.match(r"^Cerebellar_Region\d+$", label):
        return "Cerebellum"
    raise ValueError(f"Unrecognised 4S256Parcels label {label!r}")


# ── Atlas labels (FIX read from the atlas, never hard-coded) ─
def _header_labels(path: Path) -> list[str]:
    """Parcel labels from the header of an XCP-D timeseries/relmat TSV."""
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
    if len(head) == N_PARCELS + 1:      # relmat: first column is the row name
        head = head[1:]
    return head


def load_atlas_labels(source: str | Path) -> list[str]:
    """Return the 256 parcel labels in connectome / NBS node order.

    source: a table whose first two columns are 'index' and 'label'
    (e.g. atlas_4S256Parcels_labels.tsv); an XCP-D 4S256Parcels
    timeseries or relmat TSV; or a directory searched recursively for such
    TSVs (all headers must be identical). File names are never printed
    (XCP-D names carry participant identifiers).
    """
    if not source:
        raise ValueError("an atlas label source is required (--atlas-labels)")
    src = Path(source)
    if src.is_dir():
        files = sorted(src.rglob("*atlas-4S256Parcels*_timeseries.tsv"))
        files += sorted(src.rglob("*atlas-4S256Parcels*_relmat.tsv"))
        if not files:
            raise FileNotFoundError(
                "no XCP-D *atlas-4S256Parcels*_timeseries/relmat.tsv under "
                "the --atlas-labels directory")
        headers = {tuple(_header_labels(f)) for f in files}
        if len(headers) != 1:
            raise ValueError(f"{len(headers)} different 4S256Parcels headers "
                             f"among {len(files)} XCP-D files")
        labels = list(next(iter(headers)))
    else:
        with open(src) as fh:
            first = fh.readline().rstrip("\n").split("\t")
            if first[:2] == ["index", "label"]:
                rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
                if [int(r[0]) for r in rows] != list(range(len(rows))):
                    raise ValueError(f"{src.name}: index column is not 0..n-1")
                labels = [r[1] for r in rows]
            else:
                labels = _header_labels(src)
    if len(labels) != N_PARCELS:
        raise ValueError(f"expected {N_PARCELS} labels, got {len(labels)}")
    if len(set(labels)) != N_PARCELS:
        raise ValueError("duplicate atlas labels")
    nets = [network_of(lbl) for lbl in labels]          # raises on unknown
    if any(n not in CORTICAL_NETWORKS for n in nets[:N_CORTICAL]) or \
       any(n in CORTICAL_NETWORKS for n in nets[N_CORTICAL:]):
        raise ValueError("4S256Parcels order violated: parcels 0-199 must be "
                         "Schaefer cortex and 200-255 subcortex/cerebellum")
    return labels


def get_labels(source: str | Path) -> list[str]:
    """Backward-compatible name; see load_atlas_labels."""
    return load_atlas_labels(source)


def check_labels_against_decoded(nbs_dir: Path, labels: list[str]) -> int:
    """Compare every 'index:label' pair written by decode_nbs_components.py.

    Returns the number of pairs checked; raises on any disagreement.
    """
    n = 0
    for f in sorted(Path(nbs_dir).glob("nbs_*_decoded.txt")):
        for m in re.finditer(r"(\d+):([A-Za-z0-9_\-]+)", f.read_text()):
            i, lab = int(m.group(1)), m.group(2)
            if i >= len(labels) or labels[i] != lab:
                raise ValueError(f"{f.name}: node {i} is {lab!r} in the "
                                 f"decoder output but "
                                 f"{labels[i] if i < len(labels) else None!r} "
                                 "in the atlas labels")
            n += 1
    return n


def load_significant_components(json_path: Path, p_thresh: float = 0.05):
    """(contrast, [(sign, sorted nodes, edges, n_edges, p), ...]) of one
    NBS JSON, in JSON order (positive components first)."""
    data = json.load(open(json_path))
    comps = []
    for sign, key in [("pos", "obs_components_pos"),
                      ("neg", "obs_components_neg")]:
        for comp in data.get(key, []):
            if comp["p"] < p_thresh:
                # pre-fix NBS JSONs have no edge list;
                # fail with instructions rather than a bare KeyError.
                if "edges" not in comp:
                    raise KeyError(
                        f"{Path(json_path).name}: component lacks an 'edges' "
                        "list (older NBS export without an edge list). Re-run the "
                        "NBS extraction so each component records its "
                        "supra-threshold edges; do not reconstruct "
                        "them from the node list.")
                comps.append((sign, sorted(comp["nodes"]),
                              comp["edges"], comp["n_edges"], comp["p"]))
    return data.get("contrast"), comps, data


# ======================================================================
# Circos drawing
# ======================================================================
def _build_parcel_layout(labels: list[str], sub_centre_deg: float = 0.0,
                         gap_ctx_deg: float = 1.5, gap_sub_deg: float = 1.6,
                         gap_block_deg: float = 4.0):
    """Assign each parcel an angular position, grouped by segment.

    Segment arcs are proportional to the number of parcels (atlas order
    within a segment). The subcortical block is centred on sub_centre_deg
    (0 = 3 o'clock); segments run counter-clockwise in NETWORK_ORDER
    (subcortex bottom-to-top, then cortex from Visual to Default).

    Returns:
        angles: ndarray (N,) — angle (radians) of each parcel
        net_arcs: list of (net_name, start_angle, end_angle) in radians
        net_of_idx: list[str] — segment name of each parcel index
    """
    n = len(labels)
    net_of_idx = [network_of(lbl) for lbl in labels]
    groups = {net: [i for i in range(n) if net_of_idx[i] == net]
              for net in NETWORK_ORDER}
    ordered = [net for net in NETWORK_ORDER if groups[net]]
    ctx = [net for net in ordered if net in CORTICAL_NETWORKS]
    sub = [net for net in ordered if net not in CORTICAL_NETWORKS]

    total_gap = (gap_ctx_deg * max(len(ctx) - 1, 0)
                 + gap_sub_deg * max(len(sub) - 1, 0)
                 + (2 * gap_block_deg if ctx and sub else 0.0))
    per = (360.0 - total_gap) / n                     # degrees per parcel
    n_sub = sum(len(groups[net]) for net in sub)
    sub_span = n_sub * per + gap_sub_deg * max(len(sub) - 1, 0)

    angles = np.zeros(n)
    net_arcs = []
    cursor = sub_centre_deg - sub_span / 2.0

    def place(net):
        nonlocal cursor
        idxs = groups[net]
        start = cursor
        for k, i in enumerate(idxs):
            angles[i] = np.radians(start + (k + 0.5) * per)
        end = start + len(idxs) * per
        net_arcs.append((net, np.radians(start), np.radians(end)))
        cursor = end

    for j, net in enumerate(sub):
        place(net)
        if j < len(sub) - 1:
            cursor += gap_sub_deg
    if sub and ctx:
        cursor += gap_block_deg
    for j, net in enumerate(ctx):
        place(net)
        if j < len(ctx) - 1:
            cursor += gap_ctx_deg
    return angles, net_arcs, net_of_idx


def segment_at_angle(theta: float, net_arcs) -> str | None:
    """Ring segment whose arc contains angle theta (radians), else None."""
    for net, start, end in net_arcs:
        d = (theta - start) % (2 * np.pi)
        if d <= (end - start) + 1e-12:
            return net
    return None


def _readable(color: str, min_contrast: float = 4.5) -> str:
    """Darken a colour until its contrast ratio against white >= min_contrast
    (label text only; ring segments keep the network colour)."""
    import colorsys

    def lum(rgb):
        c = [x / 12.92 if x <= 0.03928 else ((x + 0.055) / 1.055) ** 2.4
             for x in rgb]
        return 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]

    r, g, b = to_rgb(color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    while 1.05 / (lum(colorsys.hls_to_rgb(h, l, s)) + 0.05) < min_contrast \
            and l > 0:
        l = max(0.0, l - 0.01)
    return to_hex(colorsys.hls_to_rgb(h, l, s))


def _spread(desired: list[float], pitch: float) -> list[float]:
    """Positions closest (least squares) to `desired` (ascending) with
    consecutive gaps >= pitch (isotonic regression, pool-adjacent-violators)."""
    t = [d - k * pitch for k, d in enumerate(desired)]
    blocks: list[list[float]] = []
    for v in t:
        blocks.append([v, 1.0])
        while len(blocks) > 1 and \
                blocks[-2][0] / blocks[-2][1] > blocks[-1][0] / blocks[-1][1]:
            s, c = blocks.pop()
            blocks[-1][0] += s
            blocks[-1][1] += c
    z: list[float] = []
    for s, c in blocks:
        z += [s / c] * int(c)
    return [zk + k * pitch for k, zk in enumerate(z)]


def _data_bbox(ax, artist, renderer):
    bb = artist.get_window_extent(renderer)
    (x0, y0), (x1, y1) = ax.transData.inverted().transform(
        [[bb.x0, bb.y0], [bb.x1, bb.y1]])
    return x0, y0, x1, y1


def _boxes_overlap(a, b, pad=0.0) -> bool:
    return (a[0] < b[2] - pad and b[0] < a[2] - pad and
            a[1] < b[3] - pad and b[1] < a[3] - pad)


def _segment_hits_box(p, q, box, n=80) -> bool:
    for s in np.linspace(0.0, 1.0, n):
        x = p[0] + s * (q[0] - p[0])
        y = p[1] + s * (q[1] - p[1])
        if box[0] <= x <= box[2] and box[1] <= y <= box[3]:
            return True
    return False


def _segments_cross(p1, p2, q1, q2) -> bool:
    """Proper intersection of segments p1-p2 and q1-q2."""
    def orient(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    d1, d2 = orient(q1, q2, p1), orient(q1, q2, p2)
    d3, d4 = orient(p1, p2, q1), orient(p1, p2, q2)
    return (d1 * d2 < 0) and (d3 * d4 < 0)


def _ring_distance(box) -> float:
    """Distance from the ring centre to the nearest point of a box."""
    return float(np.hypot(np.clip(0.0, box[0], box[2]),
                          np.clip(0.0, box[1], box[3])))


def draw_circos(ax, angles, net_arcs, net_of_idx, labels,
                sig_components, title_str="", *,
                fontsize: float = 6.0, fontweight: str = "bold",
                band: float = 0.07, band_alpha: float = 0.8,
                r_label: float = 1.06, ring_clearance: float = 0.035,
                column_x: float = 1.14, line_pitch: float = 1.10,
                node_size: float = 1.8, node_edge_lw: float = 0.2,
                edge_lw_scale: float = 0.35, tick_lw: float = 0.2,
                band_lw: float = 0.2, leader_lw: float = 0.3,
                check: bool = True) -> dict:
    """Draw one circos panel (ring radius = 1 data unit).

    sig_components: list of (sign, nodes, edges, n_edges, p).
    The caller fixes the axes position and data limits (1 data unit = the
    ring radius, equal scale on x and y) BEFORE calling, because label
    positions are measured. Sizes (fontsize, node_size, linewidths) are in
    points of the figure as authored; edge widths are the pre-fix widths
    times edge_lw_scale (0.35 = the 457 -> 160 mm reduction of the pre-fix
    figure, i.e. the same printed widths).
    Returns a dict describing what was drawn.
    """
    n = len(labels)
    R_outer = 1.0
    R_inner = R_outer - band
    R_node = R_inner - 0.025
    R_edge = R_node - 0.03

    ax.axis("off")
    ax.apply_aspect()
    fig = ax.figure
    renderer = fig.canvas.get_renderer()

    sig_nodes = set()
    for sign, nodes, edges, ne, p in sig_components:
        sig_nodes.update(nodes)

    # ── Ring segments (arc length proportional to parcel count) ──
    from matplotlib.patches import Polygon
    for net, start, end in net_arcs:
        theta = np.linspace(start, end, max(12, int(np.degrees(end - start) * 4)))
        outer = np.c_[np.cos(theta), np.sin(theta)] * R_outer
        inner = (np.c_[np.cos(theta), np.sin(theta)] * R_inner)[::-1]
        ax.add_patch(Polygon(np.vstack([outer, inner]), closed=True,
                             facecolor=NET_COLORS[net], alpha=band_alpha,
                             edgecolor="k", linewidth=band_lw, zorder=2))

    # ── Node ticks / significant node dots (colour = ring segment) ──
    node_xy, node_segment = {}, {}
    for i in range(n):
        a = angles[i]
        if i in sig_nodes:
            x, y = R_node * np.cos(a), R_node * np.sin(a)
            ax.plot(x, y, "o", color=NET_COLORS[net_of_idx[i]],
                    markersize=node_size, markeredgecolor="k",
                    markeredgewidth=node_edge_lw, zorder=4)
            node_xy[i] = (float(x), float(y))
            node_segment[i] = segment_at_angle(a, net_arcs)
        else:
            ax.plot([R_inner * np.cos(a), (R_inner - 0.022) * np.cos(a)],
                    [R_inner * np.sin(a), (R_inner - 0.022) * np.sin(a)],
                    color="#BBBBBB", linewidth=tick_lw, zorder=1)

    # ── Edges (bezier curves) — actual NBS supra-threshold edges only ──
    edges_drawn = []
    for ci, (sign, nodes, edges, n_edges, p) in enumerate(sig_components):
        # Guard: refuse to draw anything but the true edge
        # list. A clique over the component's nodes (the older
        # bug) or any stale edge list fails loudly here.
        if len(edges) != n_edges:
            raise ValueError(
                f"Component ({sign}, {len(nodes)} nodes) carries "
                f"{len(edges)} edges but records n_edges={n_edges}. "
                "Refusing to draw: the edge list must be the NBS "
                "supra-threshold edges, never an all-pairs clique.")
        color = POS_EDGE_COLOR if sign == "pos" else NEG_EDGE_COLOR
        alpha = min(0.85, max(0.25, 8 / max(len(edges), 1)))
        lw = edge_lw_scale * max(0.4, min(2.2, -np.log10(p)))
        for ni, nj in edges:
            ai, aj = angles[ni], angles[nj]
            xi, yi = R_edge * np.cos(ai), R_edge * np.sin(ai)
            xj, yj = R_edge * np.cos(aj), R_edge * np.sin(aj)
            # Bezier control point pulled toward center
            ang_diff = abs(ai - aj) % (2 * np.pi)
            if ang_diff > np.pi:
                ang_diff = 2 * np.pi - ang_diff
            pull = 0.1 + 0.4 * (ang_diff / np.pi)
            mx, my = pull * (xi + xj) / 2, pull * (yi + yj) / 2
            path = MPath([(xi, yi), (mx, my), (xj, yj)],
                         [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3])
            ax.add_patch(mpatches.PathPatch(
                path, facecolor="none", edgecolor=color,
                linewidth=lw, alpha=alpha, zorder=3))
            edges_drawn.append((ci, int(ni), int(nj)))

    # ── Segment labels ──
    ctx_texts, sub_items = [], []
    for net, start, end in net_arcs:
        mid = 0.5 * (start + end)
        colr = _readable(NET_COLORS[net])
        if net in CORTICAL_NETWORKS:
            c, s = np.cos(mid), np.sin(mid)
            ha = "left" if c > 0.3 else ("right" if c < -0.3 else "center")
            va = "bottom" if s > 0.3 else ("top" if s < -0.3 else "center")
            t = ax.text(r_label * c, r_label * s, NET_SHORT[net], ha=ha,
                        va=va, fontsize=fontsize, fontweight=fontweight,
                        color=colr, zorder=5)
            ctx_texts.append((net, t, mid))
        else:
            t = ax.text(column_x, r_label * np.sin(mid), NET_SHORT[net],
                        ha="left", va="center", fontsize=fontsize,
                        fontweight=fontweight, color=colr, zorder=5)
            sub_items.append((net, t, mid))

    def boxes_of(items):
        return [(net, t, _data_bbox(ax, t, renderer)) for net, t, _ in items]

    def nudge_cortical(others):
        """Move cortical labels radially outward until they clear the ring
        and every other label box."""
        for net, t, mid in ctx_texts:
            for _ in range(80):
                me = _data_bbox(ax, t, renderer)
                clash = _ring_distance(me) < R_outer + ring_clearance or any(
                    _boxes_overlap(me, bb) for nn, tt, bb in
                    boxes_of(ctx_texts) + others if tt is not t)
                if not clash:
                    break
                x, y = t.get_position()
                t.set_position((x + 0.01 * np.cos(mid), y + 0.01 * np.sin(mid)))

    nudge_cortical([])
    ctx_boxes = [bb for _, _, bb in boxes_of(ctx_texts)]
    lo = min([b[1] for b in ctx_boxes] + [-R_outer])
    hi = max([b[3] for b in ctx_boxes] + [R_outer])

    # Subcortical labels: one column right of the ring, ordered like the
    # segments, spaced by a fixed line pitch, as close as possible to the
    # height of their segment, kept within the vertical extent of the
    # cortical labels when possible; leader lines join label and segment.
    probe = ax.text(0, 0, "Ag", fontsize=fontsize, fontweight=fontweight)
    _, py0, _, py1 = _data_bbox(ax, probe, renderer)
    probe.remove()
    h = py1 - py0
    pitch = h * line_pitch
    order = sorted(range(len(sub_items)), key=lambda k: np.sin(sub_items[k][2]))
    ys = _spread([r_label * np.sin(sub_items[k][2]) for k in order], pitch)
    if ys:
        top_c, bot_c = ys[-1] + h / 2, ys[0] - h / 2
        if top_c - bot_c > hi - lo:
            shift = 0.5 * (hi + lo) - 0.5 * (top_c + bot_c)
        elif top_c > hi:
            shift = hi - top_c
        elif bot_c < lo:
            shift = lo - bot_c
        else:
            shift = 0.0
        ys = [y + shift for y in ys]
    leaders = []
    for k, y in zip(order, ys):
        net, t, mid = sub_items[k]
        t.set_y(y)
        p0 = ((R_outer + 0.012) * np.cos(mid), (R_outer + 0.012) * np.sin(mid))
        p1 = ((R_outer + 0.045) * np.cos(mid), (R_outer + 0.045) * np.sin(mid))
        p2 = (column_x - 0.025, y)
        ax.plot([p0[0], p1[0], p2[0]], [p0[1], p1[1], p2[1]],
                color="#6E6E6E", linewidth=leader_lw, zorder=1.5,
                solid_capstyle="butt")
        leaders.append((net, [p0, p1, p2]))

    nudge_cortical(boxes_of(sub_items))

    # ── Collision checks (labels vs labels, leaders, ring) ──
    boxes = boxes_of(ctx_texts) + boxes_of(sub_items)
    problems = []
    if check:
        for (na, ta, ba), (nb, tb, bb) in itertools.combinations(boxes, 2):
            if _boxes_overlap(ba, bb):
                problems.append(f"label overlap: {NET_SHORT[na]} / {NET_SHORT[nb]}")
        for net, pts in leaders:
            for nb, tb, bb in boxes:
                if nb == net:
                    continue
                if any(_segment_hits_box(pts[j], pts[j + 1], bb)
                       for j in range(len(pts) - 1)):
                    problems.append(f"leader of {NET_SHORT[net]} crosses "
                                    f"label {NET_SHORT[nb]}")
        for (na, pa), (nb, pb) in itertools.combinations(leaders, 2):
            if any(_segments_cross(pa[i], pa[i + 1], pb[j], pb[j + 1])
                   for i in range(len(pa) - 1) for j in range(len(pb) - 1)):
                problems.append(f"leaders cross: {NET_SHORT[na]} / {NET_SHORT[nb]}")
        for nb, tb, bb in boxes:
            if _ring_distance(bb) < R_outer + 0.5 * ring_clearance:
                problems.append(f"label {NET_SHORT[nb]} touches the ring")
        if problems:
            raise RuntimeError("circos label layout problems: " +
                               "; ".join(problems))

    xs = [b[0] for _, _, b in boxes] + [b[2] for _, _, b in boxes] + [-1, 1]
    ys_ = [b[1] for _, _, b in boxes] + [b[3] for _, _, b in boxes] + [-1, 1]
    if title_str:
        ax.set_title(title_str, fontsize=fontsize + 1.5, fontweight="bold")
    return {
        "label_texts": [t for _, t, _ in ctx_texts] + [t for _, t, _ in sub_items],
        "label_boxes": {NET_SHORT[net]: tuple(map(float, bb)) for net, _, bb in boxes},
        "leaders": leaders,
        "node_xy": node_xy,
        "node_segment": node_segment,
        "edges_drawn": edges_drawn,
        "block_bbox": (min(xs), min(ys_), max(xs), max(ys_)),
        "layout_problems": problems,
    }


def audit_texts(fig, texts, min_pt: float = 6.0, pad_px: float = 0.0):
    """Printed-size text audit for a figure authored at its print size.

    Returns a list of problems: text outside the canvas, text smaller than
    min_pt, or overlapping text boxes.
    """
    renderer = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    boxes, problems = [], []
    for t in texts:
        if not t.get_visible() or not t.get_text().strip():
            continue
        bb = t.get_window_extent(renderer)
        if bb.x0 < 0 or bb.y0 < 0 or bb.x1 > W or bb.y1 > H:
            problems.append(f"outside canvas: {t.get_text()!r}")
        if t.get_fontsize() < min_pt - 1e-9:
            problems.append(f"{t.get_fontsize():.2f} pt < {min_pt} pt: "
                            f"{t.get_text()!r}")
        boxes.append((t, (bb.x0, bb.y0, bb.x1, bb.y1)))
    for (a, ba), (b, bb) in itertools.combinations(boxes, 2):
        if _boxes_overlap(ba, bb, pad=pad_px):
            problems.append(f"overlap: {a.get_text()!r} / {b.get_text()!r}")
    return problems


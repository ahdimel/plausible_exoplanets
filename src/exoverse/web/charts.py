"""Server-rendered SVG charts for the web UI.

Follows the dataviz method: thin marks, hairline grid, recessive axes, muted
ink for labels, categorical hues assigned to entities in fixed slot order
(never cycled), hover tooltips via data-tip attributes handled in base.html.
Colors are CSS custom properties defined in base.html so light/dark modes
swap in one place; charts reference roles (var(--series-1)) not raw hex.
"""

from __future__ import annotations

import math
from typing import List, Sequence

# Fixed entity -> categorical slot assignment (never re-ordered by rank)
CLASS_SLOTS = {
    "rocky": 1,
    "sub-neptune": 2,
    "neptunian": 3,
    "giant": 5,
}


def _esc(s: str) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def _ticks_log(lo: float, hi: float) -> List[float]:
    lo_e = math.floor(math.log10(lo))
    hi_e = math.ceil(math.log10(hi))
    return [10.0 ** e for e in range(int(lo_e), int(hi_e) + 1)
            if lo <= 10.0 ** e <= hi]


def _fmt(v: float) -> str:
    if v >= 1000:
        return f"{v:,.0f}"
    if v >= 10:
        return f"{v:.0f}"
    if v >= 1:
        return f"{v:.1f}".rstrip("0").rstrip(".")
    return f"{v:g}"


def scatter_log_log(points: Sequence[dict], width=760, height=440,
                    x_label="Orbital period (days)",
                    y_label="Planet radius (R⊕)") -> str:
    """points: dicts with x, y, cls (composition class), tip, href, ring
    (bool: highlight ring for transiting planets)."""
    if not points:
        return "<p class='muted'>No data.</p>"
    ml, mr, mt, mb = 64, 16, 12, 46
    pw, ph = width - ml - mr, height - mt - mb
    xs = [p["x"] for p in points]
    ys = [p["y"] for p in points]
    x0, x1 = min(xs) * 0.8, max(xs) * 1.25
    y0, y1 = min(ys) * 0.8, max(ys) * 1.25

    def sx(v):
        return ml + (math.log10(v) - math.log10(x0)) / \
            (math.log10(x1) - math.log10(x0)) * pw

    def sy(v):
        return mt + ph - (math.log10(v) - math.log10(y0)) / \
            (math.log10(y1) - math.log10(y0)) * ph

    parts = [f'<svg viewBox="0 0 {width} {height}" role="img" '
             f'aria-label="Planet radius vs orbital period scatter">']
    for tv in _ticks_log(x0, x1):
        x = sx(tv)
        parts.append(f'<line x1="{x:.1f}" y1="{mt}" x2="{x:.1f}" '
                     f'y2="{mt+ph}" class="grid"/>')
        parts.append(f'<text x="{x:.1f}" y="{mt+ph+18}" class="tick" '
                     f'text-anchor="middle">{_fmt(tv)}</text>')
    for tv in _ticks_log(y0, y1):
        y = sy(tv)
        parts.append(f'<line x1="{ml}" y1="{y:.1f}" x2="{ml+pw}" '
                     f'y2="{y:.1f}" class="grid"/>')
        parts.append(f'<text x="{ml-8}" y="{y+4:.1f}" class="tick" '
                     f'text-anchor="end">{_fmt(tv)}</text>')
    parts.append(f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" '
                 f'class="axis"/>')
    for p in points:
        x, y = sx(p["x"]), sy(p["y"])
        slot = CLASS_SLOTS.get(p.get("cls", ""), 1)
        ring = (' stroke="var(--surface-1)" stroke-width="1.5"'
                if not p.get("ring") else
                ' stroke="var(--text-primary)" stroke-width="1.5"')
        dot = (f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.5" '
               f'fill="var(--series-{slot})"{ring} '
               f'data-tip="{_esc(p.get("tip", ""))}"/>')
        if p.get("href"):
            dot = f'<a href="{_esc(p["href"])}">{dot}</a>'
        parts.append(dot)
    parts.append(f'<text x="{ml+pw/2:.0f}" y="{height-8}" class="axis-label" '
                 f'text-anchor="middle">{_esc(x_label)}</text>')
    parts.append(f'<text x="14" y="{mt+ph/2:.0f}" class="axis-label" '
                 f'text-anchor="middle" '
                 f'transform="rotate(-90 14 {mt+ph/2:.0f})">{_esc(y_label)}</text>')
    parts.append("</svg>")
    return "".join(parts)


def hbar(rows: Sequence[tuple], width=760, row_h=34, label_w=210,
         value_suffix="") -> str:
    """rows: (label, value, tip). Single measure -> single hue (series-1)."""
    if not rows:
        return "<p class='muted'>No data.</p>"
    vmax = max(v for _, v, _ in rows) or 1
    height = row_h * len(rows) + 8
    pw = width - label_w - 70
    parts = [f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Bar chart">']
    for i, (label, v, tip) in enumerate(rows):
        y = i * row_h + 6
        bw = max(v / vmax * pw, 1.5)
        parts.append(f'<text x="{label_w-10}" y="{y+row_h/2+1}" class="tick" '
                     f'text-anchor="end" dominant-baseline="middle">'
                     f'{_esc(label)}</text>')
        parts.append(f'<rect x="{label_w}" y="{y+6}" width="{bw:.1f}" '
                     f'height="{row_h-16}" rx="3" fill="var(--series-1)" '
                     f'data-tip="{_esc(tip)}"/>')
        parts.append(f'<text x="{label_w+bw+8:.1f}" y="{y+row_h/2+1}" '
                     f'class="value" dominant-baseline="middle">'
                     f'{_fmt(v)}{value_suffix}</text>')
    parts.append("</svg>")
    return "".join(parts)


def overlay_step_hist(series: Sequence[dict], lo: float, hi: float,
                      n_bins=28, width=760, height=300,
                      x_label="") -> str:
    """Two-series comparison as step-outline histograms (density-normalized).
    series: dicts with name, values, slot (categorical slot number)."""
    if not series or all(not s["values"] for s in series):
        return "<p class='muted'>No data.</p>"
    ml, mr, mt, mb = 56, 16, 34, 44   # top margin reserves the legend band
    pw, ph = width - ml - mr, height - mt - mb
    bw = (hi - lo) / n_bins
    dens = []
    for s in series:
        counts = [0] * n_bins
        for v in s["values"]:
            if lo <= v < hi:
                counts[int((v - lo) / bw)] += 1
        total = max(sum(counts), 1)
        dens.append([c / total / bw for c in counts])
    dmax = max(max(d) for d in dens) or 1.0

    def sx(i):
        return ml + i / n_bins * pw

    def sy(v):
        return mt + ph - v / dmax * ph

    parts = [f'<svg viewBox="0 0 {width} {height}" role="img" '
             f'aria-label="Distribution comparison">']
    parts.append(f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" '
                 f'class="axis"/>')
    n_xticks = 6
    for k in range(n_xticks + 1):
        xv = lo + (hi - lo) * k / n_xticks
        parts.append(f'<text x="{ml + pw*k/n_xticks:.1f}" y="{mt+ph+18}" '
                     f'class="tick" text-anchor="middle">{xv:.1f}</text>')
    for s, d in zip(series, dens):
        pts = []
        for i, v in enumerate(d):
            pts.append(f"{sx(i):.1f},{sy(v):.1f}")
            pts.append(f"{sx(i+1):.1f},{sy(v):.1f}")
        parts.append(f'<polyline points="{" ".join(pts)}" fill="none" '
                     f'stroke="var(--series-{s["slot"]})" stroke-width="2" '
                     f'data-tip="{_esc(s["name"])} (n={len(s["values"])})"/>')
    # Legend (2 series -> always present), in the reserved band above the plot
    lx = ml + 12
    for s in series:
        parts.append(f'<rect x="{lx}" y="{mt-20}" width="14" height="4" rx="2" '
                     f'fill="var(--series-{s["slot"]})"/>')
        parts.append(f'<text x="{lx+20}" y="{mt-13}" class="tick">'
                     f'{_esc(s["name"])}</text>')
        lx += 44 + 7.2 * len(s["name"])
    parts.append(f'<text x="{ml+pw/2:.0f}" y="{height-6}" class="axis-label" '
                 f'text-anchor="middle">{_esc(x_label)}</text>')
    parts.append("</svg>")
    return "".join(parts)


def light_curve_svg(t_hours: Sequence[float], flux: Sequence[float],
                    depth_ppm: float, width=760, height=300) -> str:
    ml, mr, mt, mb = 72, 16, 14, 44
    pw, ph = width - ml - mr, height - mt - mb
    t0, t1 = min(t_hours), max(t_hours)
    dppm = [(f - 1.0) * 1e6 for f in flux]
    y0, y1 = min(dppm) * 1.12, max(max(dppm), 0.0) + abs(min(dppm)) * 0.06

    def sx(v):
        return ml + (v - t0) / (t1 - t0) * pw

    def sy(v):
        return mt + ph - (v - y0) / (y1 - y0) * ph

    parts = [f'<svg viewBox="0 0 {width} {height}" role="img" '
             f'aria-label="Transit light curve, depth {depth_ppm:.0f} ppm">']
    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        yv = y0 + (y1 - y0) * frac
        y = sy(yv)
        parts.append(f'<line x1="{ml}" y1="{y:.1f}" x2="{ml+pw}" y2="{y:.1f}" '
                     f'class="grid"/>')
        parts.append(f'<text x="{ml-8}" y="{y+4:.1f}" class="tick" '
                     f'text-anchor="end">{yv:,.0f}</text>')
    for tv in range(int(math.floor(t0)), int(math.ceil(t1)) + 1):
        if t0 <= tv <= t1:
            parts.append(f'<text x="{sx(tv):.1f}" y="{mt+ph+18}" class="tick" '
                         f'text-anchor="middle">{tv:+d}h</text>')
    pts = " ".join(f"{sx(t):.1f},{sy(d):.1f}" for t, d in zip(t_hours, dppm))
    parts.append(f'<polyline points="{pts}" fill="none" '
                 f'stroke="var(--series-1)" stroke-width="2"/>')
    parts.append(f'<text x="{ml+pw/2:.0f}" y="{height-6}" class="axis-label" '
                 f'text-anchor="middle">Time from mid-transit (hours)</text>')
    parts.append(f'<text x="16" y="{mt+ph/2:.0f}" class="axis-label" '
                 f'text-anchor="middle" transform="rotate(-90 16 {mt+ph/2:.0f})">'
                 f'Flux change (ppm)</text>')
    parts.append("</svg>")
    return "".join(parts)

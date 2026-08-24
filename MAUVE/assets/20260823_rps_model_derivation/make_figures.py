#!/usr/bin/env python3
"""Create dependency-free SVG figures for the 2026-08-23 RPS derivation."""

from __future__ import annotations

import html
import math
from pathlib import Path


OUT = Path(__file__).resolve().parent
# Quick Look's SVG thumbnailer uses a square viewport.  Keep the authored
# content in the upper 900 units of a square canvas, then trim the raster.
W, H = 1500, 1500
COL = {
    "navy": "#16324F",
    "blue": "#2F6690",
    "cyan": "#3A8D9D",
    "green": "#4A7C59",
    "gold": "#D5A021",
    "orange": "#C76D3A",
    "red": "#A23E48",
    "plum": "#6D597A",
    "ink": "#20252A",
    "muted": "#5C646B",
    "paper": "#F7F4EE",
    "grid": "#D8D4CB",
    "pale_blue": "#EAF1F5",
    "pale_red": "#F8E9EA",
    "pale_green": "#EBF3ED",
}


def esc(value: str) -> str:
    return html.escape(value, quote=True)


def svg_start(title: str, subtitle: str = "") -> list[str]:
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<defs>',
        '<marker id="arrow" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto" markerUnits="strokeWidth"><path d="M0,0 L0,6 L9,3 z" fill="#5C646B"/></marker>',
        '<marker id="arrowBlue" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto" markerUnits="strokeWidth"><path d="M0,0 L0,6 L9,3 z" fill="#2F6690"/></marker>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#20252A}.title{font-size:34px;font-weight:700;fill:#16324F}.subtitle{font-size:19px;fill:#5C646B}.boxTitle{font-size:21px;font-weight:700}.body{font-size:17px}.small{font-size:15px;fill:#5C646B}.axis{font-size:17px;fill:#20252A}.label{font-size:18px;font-weight:700}</style>',
        '</defs>',
        f'<text class="title" x="70" y="65">{esc(title)}</text>',
    ]
    if subtitle:
        parts.append(f'<text class="subtitle" x="70" y="98">{esc(subtitle)}</text>')
    return parts


def end(parts: list[str], stem: str) -> None:
    parts.append('</svg>')
    (OUT / f"{stem}.svg").write_text("\n".join(parts) + "\n", encoding="utf-8")


def line(parts: list[str], x1: float, y1: float, x2: float, y2: float, *, colour: str = "#5C646B", width: float = 3, arrow: bool = False, dash: str | None = None) -> None:
    attrs = f'stroke="{colour}" stroke-width="{width}" fill="none"'
    if arrow:
        attrs += ' marker-end="url(#arrow)"'
    if dash:
        attrs += f' stroke-dasharray="{dash}"'
    parts.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" {attrs}/>' )


def box(parts: list[str], x: float, y: float, w: float, h: float, title: str, lines: list[str], colour: str, fill: str = "white") -> None:
    parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="18" fill="{fill}" stroke="{colour}" stroke-width="3"/>')
    parts.append(f'<text class="boxTitle" x="{x + 24}" y="{y + 38}" fill="{colour}">{esc(title)}</text>')
    yy = y + 72
    for value in lines:
        parts.append(f'<text class="body" x="{x + 24}" y="{yy}">{esc(value)}</text>')
        yy += 27


def figure_architecture() -> None:
    p = svg_start(
        "Observation-to-stripping model architecture",
        "Every arrow carries an assumption; branch structure is retained instead of collapsed.",
    )
    specs = [
        (65, 155, 250, 185, "Observed cluster", ["centre, distance", "R_p, Delta v_los", "ICM beta profile", "mass profile"], COL["navy"], COL["pale_blue"]),
        (375, 155, 250, 185, "Radial deprojection", ["choose energy/apocentre", "solve implicit r roots", "infall/outgoing linked", "retain all valid roots"], COL["blue"], "white"),
        (685, 155, 250, 185, "Three-dimensional wind", ["v_gal = epsilon v r_hat", "u_0 = -v_gal", "P_ram = rho |u_0|^2", "optional local rotation"], COL["cyan"], COL["pale_blue"]),
        (995, 155, 250, 185, "Observed galaxy", ["q, q_0, PA", "Sigma_g, sigma_z or h_g", "stellar/bulge/halo mass", "rotation curve"], COL["green"], COL["pale_green"]),
    ]
    for s in specs:
        box(p, *s)
    for x1, x2 in [(315, 375), (625, 685), (935, 995)]:
        line(p, x1 + 8, 248, x2 - 8, 248, arrow=True)
    box(p, 170, 430, 320, 180, "Finite-thickness geometry", ["rho_g(R,z) from HSE", "Sigma_w(b) = integral rho ds", "finite at exact edge-on", "coherence is a closure"], COL["plum"], "white")
    box(p, 590, 430, 320, 180, "Galaxy response", ["face-on long-pulse force", "arbitrary-direction impulse", "or coupled vector ODE", "never freeze L_z if w_phi != 0"], COL["orange"], "white")
    box(p, 1010, 430, 320, 180, "Reported products", ["ray-space susceptibility", "R_strip(phi) only with mapping", "branch-weighted uncertainty", "validity and omission flags"], COL["red"], COL["pale_red"])
    line(p, 810, 340, 350, 430, arrow=True)
    line(p, 810, 340, 750, 430, arrow=True)
    line(p, 1120, 340, 750, 430, arrow=True)
    line(p, 490, 520, 590, 520, arrow=True)
    line(p, 910, 520, 1010, 520, arrow=True)
    p.append(f'<rect x="110" y="690" width="1280" height="125" rx="20" fill="{COL["pale_red"]}" stroke="{COL["red"]}" stroke-width="3"/>')
    p.append(f'<text class="boxTitle" x="145" y="730" fill="{COL["red"]}">Critical validity split</text>')
    p.append('<text class="body" x="145" y="766">Pressure-only loss of equilibrium is a long-pulse, torque-free result. An arbitrary azimuthal wind changes L_z.</text>')
    p.append('<text class="body" x="145" y="797">Use the exact vector impulse threshold for a short pulse, or integrate the vector equations for a finite-duration pressure history.</text>')
    end(p, "figure_01_model_architecture")


def figure_wind_column() -> None:
    p = svg_start(
        "Wind-aligned column geometry",
        "Finite thickness regularizes exact edge-on incidence, but a complete chord is a coherent-column assumption.",
    )
    # Oblique panel.
    p.append(f'<rect x="70" y="140" width="645" height="615" rx="18" fill="{COL["paper"]}" stroke="{COL["grid"]}" stroke-width="2"/>')
    p.append(f'<text class="boxTitle" x="105" y="185" fill="{COL["blue"]}">A. Oblique incidence</text>')
    p.append(f'<ellipse cx="390" cy="485" rx="245" ry="75" fill="{COL["pale_blue"]}" stroke="{COL["navy"]}" stroke-width="4"/>')
    p.append(f'<ellipse cx="390" cy="485" rx="245" ry="25" fill="none" stroke="{COL["blue"]}" stroke-width="2" stroke-dasharray="8 7"/>')
    line(p, 160, 650, 610, 260, colour=COL["orange"], width=7, arrow=True)
    line(p, 390, 485, 390, 250, colour=COL["green"], width=4, arrow=True)
    p.append('<text class="label" x="525" y="315">wind w_hat</text>')
    p.append('<text class="label" x="410" y="275">disk normal n_hat</text>')
    p.append('<path d="M390,355 A130,130 0 0,1 485,397" fill="none" stroke="#6D597A" stroke-width="4"/>')
    p.append('<text class="label" x="462" y="365">theta</text>')
    p.append('<circle cx="390" cy="485" r="9" fill="#A23E48"/>')
    p.append('<text class="body" x="110" y="690">For |w_z| &gt; 0, each ray has a unique midplane crossing.</text>')
    p.append('<text class="small" x="110" y="722">Thin limit: Sigma_w -> Sigma_g / |cos(theta)|.</text>')
    # Edge-on panel.
    p.append(f'<rect x="785" y="140" width="645" height="615" rx="18" fill="{COL["paper"]}" stroke="{COL["grid"]}" stroke-width="2"/>')
    p.append(f'<text class="boxTitle" x="820" y="185" fill="{COL["red"]}">B. Exactly edge-on</text>')
    p.append(f'<rect x="875" y="410" width="450" height="115" rx="55" fill="{COL["pale_green"]}" stroke="{COL["green"]}" stroke-width="4"/>')
    p.append(f'<ellipse cx="1100" cy="467" rx="225" ry="57" fill="none" stroke="{COL["green"]}" stroke-width="2" stroke-dasharray="8 7"/>')
    line(p, 850, 345, 1345, 345, colour=COL["orange"], width=7, arrow=True)
    line(p, 850, 467, 1345, 467, colour=COL["plum"], width=4, arrow=True)
    p.append('<text class="label" x="1190" y="325">w_hat lies in disk</text>')
    p.append('<text class="label" x="1010" y="450">full chord</text>')
    p.append('<text class="body" x="830" y="595">Sigma_parallel(y,z) = integral rho(sqrt(x^2+y^2),z) dx</text>')
    p.append('<text class="body" x="830" y="635">Finite if h_g &gt; 0 and the radial disk is finite.</text>')
    p.append('<text class="small" x="830" y="680">The chord crosses multiple radii and may mix opposite circular velocities.</text>')
    p.append('<text class="small" x="830" y="710">It is a ray-space mass column, not automatically one material parcel.</text>')
    end(p, "figure_02_wind_column_geometry")


def adaptive_simpson(f, a: float, b: float, eps: float = 1e-10, depth: int = 20) -> float:
    def simp(x0: float, x1: float) -> float:
        xm = 0.5 * (x0 + x1)
        return (x1 - x0) * (f(x0) + 4.0 * f(xm) + f(x1)) / 6.0

    whole = simp(a, b)

    def rec(x0: float, x1: float, val: float, tol: float, left: int) -> float:
        xm = 0.5 * (x0 + x1)
        vl, vr = simp(x0, xm), simp(xm, x1)
        if left <= 0 or abs(vl + vr - val) < 15.0 * tol:
            return vl + vr + (vl + vr - val) / 15.0
        return rec(x0, xm, vl, tol / 2.0, left - 1) + rec(xm, x1, vr, tol / 2.0, left - 1)

    return rec(a, b, whole, eps, depth)


def xk1(x: float) -> float:
    if x == 0:
        return 1.0
    # x K_1(x) = integral_0^infinity x cosh(t) exp[-x cosh(t)] dt.
    f = lambda t: x * math.cosh(t) * math.exp(-x * math.cosh(t))
    return adaptive_simpson(f, 0.0, 14.0, eps=2e-9)


def plot_frame(parts: list[str], x0: float, y0: float, w: float, h: float, xmin: float, xmax: float, ymin: float, ymax: float) -> tuple:
    def tx(x: float) -> float:
        return x0 + (x - xmin) * w / (xmax - xmin)

    def ty(y: float) -> float:
        return y0 + h - (y - ymin) * h / (ymax - ymin)

    parts.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="white" stroke="{COL["ink"]}" stroke-width="2"/>')
    return tx, ty


def figure_edge_profile() -> None:
    p = svg_start(
        "Exact edge-on column of a finite double-exponential disk",
        "The normalized chord x K_1(x) remains finite at the centre and exceeds the local face-on exponential away from it.",
    )
    x0, y0, pw, ph = 145, 165, 1100, 565
    tx, ty = plot_frame(p, x0, y0, pw, ph, 0, 5, 0, 1.05)
    for i in range(6):
        xx = float(i)
        line(p, tx(xx), y0, tx(xx), y0 + ph, colour=COL["grid"], width=1)
        p.append(f'<text class="axis" x="{tx(xx)}" y="{y0 + ph + 34}" text-anchor="middle">{i}</text>')
    for j in range(6):
        yy = j / 5
        line(p, x0, ty(yy), x0 + pw, ty(yy), colour=COL["grid"], width=1)
        p.append(f'<text class="axis" x="{x0 - 18}" y="{ty(yy) + 6}" text-anchor="end">{yy:.1f}</text>')
    xs = [5.0 * i / 160 for i in range(161)]
    edge = [xk1(x) for x in xs]
    face = [math.exp(-x) for x in xs]
    pts_edge = " ".join(f"{tx(x):.2f},{ty(y):.2f}" for x, y in zip(xs, edge))
    pts_face = " ".join(f"{tx(x):.2f},{ty(y):.2f}" for x, y in zip(xs, face))
    p.append(f'<polyline points="{pts_edge}" fill="none" stroke="{COL["red"]}" stroke-width="6"/>')
    p.append(f'<polyline points="{pts_face}" fill="none" stroke="{COL["blue"]}" stroke-width="5" stroke-dasharray="12 8"/>')
    p.append(f'<circle cx="{tx(0)}" cy="{ty(1)}" r="8" fill="{COL["red"]}"/>')
    p.append(f'<text class="label" x="{tx(2.8)}" y="{ty(edge[90]) - 25}" fill="{COL["red"]}">edge-on chord: x K_1(x)</text>')
    p.append(f'<text class="label" x="{tx(3.0)}" y="{ty(0.055) - 16}" fill="{COL["blue"]}">face-on local column: exp(-x)</text>')
    p.append(f'<text class="body" x="{tx(0.18)}" y="{ty(0.93)}">lim x K_1(x) = 1</text>')
    p.append(f'<text class="axis" x="{x0 + pw / 2}" y="{y0 + ph + 78}" text-anchor="middle">x = |y| / R_g</text>')
    p.append(f'<text class="axis" x="48" y="{y0 + ph / 2}" text-anchor="middle" transform="rotate(-90 48 {y0 + ph / 2})">column normalized to central edge-on value</text>')
    p.append('<text class="small" x="145" y="845">For rho = [Sigma_0/(2h_g)] exp(-R/R_g) f_z(z), Sigma_parallel = (Sigma_0/h_g)|y|K_1(|y|/R_g) f_z(z).</text>')
    end(p, "figure_03_edge_on_column_profile")


def bisect_root(level: float, lo: float, hi: float) -> float:
    f = lambda x: 1.0 / x - 1.0 / x**3 - level
    flo = f(lo)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def figure_deprojection() -> None:
    p = svg_start(
        "Radial-orbit deprojection is not generally unique",
        "Point-mass example: the LOS constraint can intersect the radial model twice even before observational uncertainties are propagated.",
    )
    x0, y0, pw, ph = 145, 165, 1100, 565
    tx, ty = plot_frame(p, x0, y0, pw, ph, 1, 6, 0, 0.42)
    for i in range(1, 7):
        line(p, tx(i), y0, tx(i), y0 + ph, colour=COL["grid"], width=1)
        p.append(f'<text class="axis" x="{tx(i)}" y="{y0 + ph + 34}" text-anchor="middle">{i}</text>')
    for j in range(8):
        yy = 0.05 * j
        line(p, x0, ty(yy), x0 + pw, ty(yy), colour=COL["grid"], width=1)
        p.append(f'<text class="axis" x="{x0 - 18}" y="{ty(yy) + 6}" text-anchor="end">{yy:.2f}</text>')
    xs = [1 + 5 * i / 400 for i in range(401)]
    ys = [1 / x - 1 / x**3 for x in xs]
    points = " ".join(f"{tx(x):.2f},{ty(y):.2f}" for x, y in zip(xs, ys))
    p.append(f'<polyline points="{points}" fill="none" stroke="{COL["navy"]}" stroke-width="6"/>')
    peak_x = math.sqrt(3)
    peak_y = 1 / peak_x - 1 / peak_x**3
    p.append(f'<circle cx="{tx(peak_x)}" cy="{ty(peak_y)}" r="9" fill="{COL["red"]}"/>')
    p.append(f'<text class="label" x="{tx(peak_x) + 25}" y="{ty(peak_y) - 18}">maximum at r/R_p = sqrt(3)</text>')
    level = 0.20
    line(p, tx(1), ty(level), tx(6), ty(level), colour=COL["orange"], width=4, dash="12 8")
    r1 = bisect_root(level, 1.000001, peak_x)
    r2 = bisect_root(level, peak_x, 6.0)
    for rr in (r1, r2):
        p.append(f'<circle cx="{tx(rr)}" cy="{ty(level)}" r="10" fill="white" stroke="{COL["orange"]}" stroke-width="5"/>')
        line(p, tx(rr), ty(level), tx(rr), ty(0), colour=COL["orange"], width=2, dash="7 7")
    p.append(f'<text class="body" x="{tx(3.4)}" y="{ty(level) - 18}">one measured Delta v_los level</text>')
    p.append(f'<text class="label" x="{tx(r1)}" y="{ty(0.035)}" text-anchor="middle">inner root</text>')
    p.append(f'<text class="label" x="{tx(r2)}" y="{ty(0.035)}" text-anchor="middle">outer root</text>')
    p.append(f'<text class="axis" x="{x0 + pw / 2}" y="{y0 + ph + 78}" text-anchor="middle">x = r / R_p</text>')
    p.append(f'<text class="axis" x="48" y="{y0 + ph / 2}" text-anchor="middle" transform="rotate(-90 48 {y0 + ph / 2})">H(x) / [2 G M f_v^2 / R_p]</text>')
    p.append('<text class="small" x="145" y="845">H(x) = 1/x - 1/x^3 for v = f_v v_esc in a point-mass potential. Above the maximum there is no radial-orbit solution.</text>')
    end(p, "figure_04_radial_deprojection_roots")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    figure_architecture()
    figure_wind_column()
    figure_edge_profile()
    figure_deprojection()


if __name__ == "__main__":
    main()

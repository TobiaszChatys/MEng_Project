"""
POD Presentation Animation  (Slide 4 — Theoretical Bridge)
===========================================================
A ~90-second Manim animation explaining Proper Orthogonal Decomposition
for the MEng viva presentation.  Designed for the presenter to narrate
over, with generous pauses between key visual beats.

Uses real eigenvalue / cumulative-energy data extracted from the L8_G12
POD results (scripts/pod_data_cache.json).

Render
------
    cd scripts
    source .venv/bin/activate
    manim -ql pod_presentation.py PODPresentation     # quick preview  (480p 15fps)
    manim -qm pod_presentation.py PODPresentation     # medium         (720p 30fps)
    manim -qh pod_presentation.py PODPresentation     # final render   (1080p 60fps)

Output lands in  scripts/media/videos/pod_presentation/<quality>/
"""

from __future__ import annotations

import json
from pathlib import Path

from manim import *
import numpy as np

# ── MATLAB-style colour palette (winter colourmap) ───────────────────
# Winter runs  [0,0,1] → [0,1,0.5]  i.e.  blue → teal
BG       = WHITE                   # presentation-ready white slide
TEXT     = "#222222"               # near-black for readability
SIGNAL   = "#0055d4"              # strong blue — organised motion / modes
TEAL     = "#008080"              # mid-winter — detail / secondary
SEAFOAM  = "#00b366"              # end-of-winter — accents
NOISE    = "#cc2222"              # red — turbulent noise / crap vectors
ENERGY   = "#d48a00"              # amber — energy / time coefficients
MUTED    = "#888888"              # soft grey — de-emphasised elements

# Set global text colour so all Text/MathTex inherit it
Tex.set_default(color=TEXT)
MathTex.set_default(color=TEXT)
Text.set_default(color=TEXT)

# ── Load real POD data ───────────────────────────────────────────────
_DATA_PATH = Path(__file__).with_name("pod_data_cache.json")
with open(_DATA_PATH) as _f:
    _POD = json.load(_f)

EIGENVALUES   = np.asarray(_POD["eigenvalues"][:60])
CUM_ENERGY    = np.asarray(_POD["cumulative_energy"][:60])
N_MODES_99    = _POD["n_modes_99pct"]        # 931 modes for 99 %
GRID_ROWS     = _POD["grid_rows"]            # 99
GRID_COLS     = _POD["grid_cols"]            # 50
TOTAL_MODES   = _POD["total_modes"]          # 3749


class PODPresentation(Scene):
    """
    Four scenes answering the three critical questions:

        1. Why are we doing POD?
        2. How does POD filtering work?   (Snapshot Method → Decomposition)
        3. What level of filtering should we use?
    """

    def construct(self):
        self.camera.background_color = BG

        self.why_pod()
        self.wipe()

        self.snapshot_method()
        self.wipe()

        self.how_pod_works()
        self.wipe()

        self.energy_ranking()

    # ── helpers ────────────────────────────────────────────────────────
    def wipe(self):
        """Fade everything out between scenes."""
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=0.7)
        self.wait(0.3)

    def section_title(self, text: str):
        """Write a bold title → pause → shrink to upper-left corner."""
        t = Text(text, font_size=56, weight="BOLD", color=SIGNAL)
        self.play(FadeIn(t, shift=UP * 0.3), run_time=0.9)
        self.wait(2)
        self.play(t.animate.scale(0.45).to_corner(UL), run_time=0.6)
        return t

    # ══════════════════════════════════════════════════════════════════
    # SCENE 1 — Why POD?                                       (~20 s)
    # ══════════════════════════════════════════════════════════════════
    def why_pod(self):
        self.section_title("Why POD?")

        ax = Axes(
            x_range=[0, 4 * PI, PI],
            y_range=[-2.5, 2.5, 1],
            x_length=10, y_length=4,
            axis_config={"color": MUTED, "stroke_width": 1.5},
        ).shift(DOWN * 0.2)

        def signal(x):
            return np.sin(x) + 0.5 * np.sin(2.5 * x)

        def noise(x):
            return (0.30 * np.sin(7.3  * x)
                  + 0.20 * np.sin(11.7 * x)
                  + 0.15 * np.sin(17.1 * x))

        raw_plot = ax.plot(lambda x: signal(x) + noise(x),
                           color=TEXT, stroke_width=2)
        raw_lbl  = Text("Measured PIV velocity", font_size=24
                        ).next_to(ax, UP, buff=0.3)

        self.play(Create(ax), FadeIn(raw_lbl), run_time=0.8)
        self.play(Create(raw_plot), run_time=1.5)
        self.wait(4)

        # Separate signal from noise — clean crossfade
        sig_plot = ax.plot(signal, color=SIGNAL, stroke_width=3)
        noi_plot = ax.plot(noise,  color=NOISE,  stroke_width=1.5)

        sig_lbl = Text("Organised wave motion", font_size=22, color=SIGNAL
                       ).next_to(ax, UP, buff=0.3).shift(LEFT * 2)
        noi_lbl = Text("Turbulent noise", font_size=22, color=NOISE
                       ).next_to(ax, UP, buff=0.3).shift(RIGHT * 2.5)

        self.play(
            FadeOut(raw_plot, run_time=0.5),
            FadeOut(raw_lbl,  run_time=0.5),
        )
        self.play(
            FadeIn(sig_plot), FadeIn(sig_lbl),
            FadeIn(noi_plot), FadeIn(noi_lbl),
            run_time=1,
        )
        self.wait(5)

    # ══════════════════════════════════════════════════════════════════
    # SCENE 2 — The Snapshot Method (Sirovich)                  (~22 s)
    # ══════════════════════════════════════════════════════════════════
    def snapshot_method(self):
        self.section_title("The Snapshot Method")

        # --- PIV arrow grid using real data ---
        piv_file = Path(__file__).with_name("piv_snapshot.npz")
        if piv_file.exists():
            piv = np.load(piv_file)
            U_raw, V_raw = piv["U"], piv["V"]
        else:
            np.random.seed(7)
            U_raw = 0.3 + 0.15 * np.random.randn(6, 6)
            V_raw = 0.1 * np.random.randn(6, 6)

        nr, nc = U_raw.shape
        arrows = VGroup()
        spacing_x, spacing_y = 0.7, 0.6
        ox = -(nc - 1) * spacing_x / 2 - 2.0
        oy = (nr - 1) * spacing_y / 2 - 0.3

        scale_factor = 2.5 / max(np.nanmax(np.abs(U_raw)), np.nanmax(np.abs(V_raw)), 1e-9)

        for i in range(nr):
            for j in range(nc):
                u, v = U_raw[i, j], V_raw[i, j]
                if np.isnan(u) or np.isnan(v) or (u == 0 and v == 0):
                    continue
                origin = np.array([ox + j * spacing_x, oy - i * spacing_y, 0])
                delta  = np.array([u * scale_factor * 0.25,
                                   v * scale_factor * 0.25, 0])
                mag = np.sqrt(u**2 + v**2)
                norm_mag = min(mag / 1.0, 1.0)
                col = interpolate_color(
                    ManimColor(SIGNAL), ManimColor(SEAFOAM), norm_mag)
                arrows.add(
                    Arrow(origin, origin + delta, buff=0, color=col,
                          stroke_width=2, max_tip_length_to_length_ratio=0.3)
                )

        border   = SurroundingRectangle(arrows, color=MUTED, buff=0.25,
                                        corner_radius=0.1)
        snap_lbl = Text("PIV snapshot  (L8 G9)", font_size=20, color=MUTED
                        ).next_to(border, DOWN, buff=0.2)

        grid_info = Text(
            f"{GRID_ROWS} × {GRID_COLS} grid   ·   U, V components",
            font_size=16, color=MUTED,
        ).next_to(snap_lbl, DOWN, buff=0.1)

        self.play(FadeIn(arrows, lag_ratio=0.02),
                  Create(border), FadeIn(snap_lbl), FadeIn(grid_info),
                  run_time=1.2)
        self.wait(3)

        # --- flatten into column vector ---
        flat_arrow = MathTex(r"\longrightarrow", font_size=48
                             ).shift(LEFT * 0.3)
        col_vec = MathTex(
            r"\mathbf{x}_k = \begin{bmatrix}"
            r" V_{x_1} \\ V_{y_1} \\ \vdots \\"
            r" V_{x_n} \\ V_{y_n} \end{bmatrix}",
            font_size=32,
        ).shift(RIGHT * 2.5)

        self.play(FadeIn(flat_arrow, shift=RIGHT * 0.2), run_time=0.4)
        self.play(FadeIn(col_vec, shift=RIGHT * 0.3), run_time=1)
        self.wait(3)

        # --- assemble the snapshot matrix ---
        piv_group = VGroup(arrows, border, snap_lbl, grid_info, flat_arrow, col_vec)
        self.play(FadeOut(piv_group), run_time=0.6)

        mat = MathTex(
            r"X", r"=",
            r"\begin{bmatrix}"
            r" | & | & & | \\"
            r" \mathbf{x}_1 & \mathbf{x}_2 & \cdots"
            r" & \mathbf{x}_{2500} \\"
            r" | & | & & |"
            r"\end{bmatrix}",
            font_size=36,
        )
        mat[0].set_color(ENERGY)

        self.play(FadeIn(mat, shift=UP * 0.3), run_time=1)

        brace     = Brace(mat[2], DOWN, color=MUTED)
        brace_lbl = Text("2,500 snapshots", font_size=22, color=MUTED
                         ).next_to(brace, DOWN, buff=0.1)

        dim_note = Text(
            f"{2 * GRID_ROWS * GRID_COLS:,} rows  (U + V stacked)",
            font_size=16, color=MUTED,
        ).next_to(mat, LEFT, buff=0.6)

        self.play(GrowFromCenter(brace), FadeIn(brace_lbl),
                  FadeIn(dim_note, shift=LEFT * 0.2), run_time=0.8)
        self.wait(4)

    # ══════════════════════════════════════════════════════════════════
    # SCENE 3 — How does POD work?  (Fourier analogy + equation) (~25 s)
    # ══════════════════════════════════════════════════════════════════
    def how_pod_works(self):
        self.section_title("How does POD work?")

        # ---- LEFT: Fourier (greyed-out — predefined) ----
        f_title = Text("Fourier", font_size=30, color=MUTED
                       ).shift(LEFT * 3.5 + UP * 2.5)
        f_sub   = Text("Predefined sine waves", font_size=18, color=MUTED
                       ).next_to(f_title, DOWN, buff=0.15)

        ax_f = Axes(
            x_range=[0, 2 * PI], y_range=[-1.3, 1.3],
            x_length=4, y_length=2.2,
            axis_config={"color": MUTED, "stroke_width": 1},
        ).shift(LEFT * 3.5 + DOWN * 0.3)

        sines = VGroup(
            ax_f.plot(lambda x: np.sin(x),
                      color=MUTED, stroke_width=2),
            ax_f.plot(lambda x: np.sin(2 * x),
                      color=MUTED, stroke_width=2, stroke_opacity=0.6),
            ax_f.plot(lambda x: np.sin(3 * x),
                      color=MUTED, stroke_width=2, stroke_opacity=0.3),
        )

        # ---- RIGHT: POD (highlighted — data-driven) ----
        p_title = Text("POD", font_size=30, weight="BOLD", color=SIGNAL
                       ).shift(RIGHT * 3.5 + UP * 2.5)
        p_sub   = Text("Modes from the data itself", font_size=18, color=TEAL
                       ).next_to(p_title, DOWN, buff=0.15)

        ax_p = Axes(
            x_range=[0, 2 * PI], y_range=[-1.3, 1.3],
            x_length=4, y_length=2.2,
            axis_config={"color": SIGNAL, "stroke_width": 1},
        ).shift(RIGHT * 3.5 + DOWN * 0.3)

        pods = VGroup(
            ax_p.plot(
                lambda x: (np.sin(x) * np.exp(-0.08 * (x - PI) ** 2)
                           + 0.3 * np.cos(1.8 * x)),
                color=SIGNAL, stroke_width=3),
            ax_p.plot(
                lambda x: 0.6 * np.sin(1.4 * x + 0.5)
                          - 0.35 * np.cos(0.9 * x),
                color=TEAL, stroke_width=2.5),
        )

        div = DashedLine(UP * 3, DOWN * 2.8,
                         color=MUTED, stroke_width=1)

        self.play(
            FadeIn(f_title, shift=DOWN * 0.2),
            FadeIn(p_title, shift=DOWN * 0.2),
            Create(div),
            run_time=0.7,
        )
        self.play(
            FadeIn(f_sub), FadeIn(p_sub),
            Create(ax_f), *[Create(s) for s in sines],
            Create(ax_p), *[Create(p) for p in pods],
            run_time=2,
        )
        self.wait(5)

        # ---- transition to the decomposition equation ----
        comp = VGroup(f_title, f_sub, ax_f, *sines,
                      p_title, p_sub, ax_p, *pods, div)
        self.play(FadeOut(comp), run_time=0.7)

        eq = MathTex(
            r"u(\mathbf{x}, t)", r"=",
            r"\sum_{k=1}^{K}",
            r"a_k(t)",          r"\;",
            r"\phi_k(\mathbf{x})",
            font_size=56,
        )
        eq[3].set_color(ENERGY)
        eq[5].set_color(SIGNAL)

        self.play(FadeIn(eq, shift=UP * 0.3), run_time=1.2)
        self.wait(1.5)

        # labelled arrows
        t_lbl = Text("Time coefficients", font_size=22, color=ENERGY)
        t_lbl.next_to(eq[3], DOWN, buff=0.7)
        t_arr = Arrow(t_lbl.get_top(), eq[3].get_bottom(),
                      color=ENERGY, stroke_width=2, buff=0.08)

        m_lbl = Text("Spatial modes", font_size=22, color=SIGNAL)
        m_lbl.next_to(eq[5], DOWN, buff=0.7)
        m_arr = Arrow(m_lbl.get_top(), eq[5].get_bottom(),
                      color=SIGNAL, stroke_width=2, buff=0.08)

        self.play(FadeIn(t_lbl, shift=UP * 0.15), GrowArrow(t_arr),
                  FadeIn(m_lbl, shift=UP * 0.15), GrowArrow(m_arr),
                  run_time=1)
        self.wait(6)

    # ══════════════════════════════════════════════════════════════════
    # SCENE 4 — What level of filtering?  (Energy ranking)      (~25 s)
    #           Uses REAL eigenvalue data from POD_Results_L8_G12.mat
    # ══════════════════════════════════════════════════════════════════
    def energy_ranking(self):
        self.section_title("What level of filtering?")

        N_SHOW = 50
        eig = EIGENVALUES[:N_SHOW]
        eig_norm = eig / eig.max()
        cum = CUM_ENERGY[:N_SHOW]

        ax = Axes(
            x_range=[0, N_SHOW + 2, 10],
            y_range=[0, 1.08, 0.2],
            x_length=9.5, y_length=4.5,
            axis_config={"color": TEXT, "stroke_width": 1.5,
                         "include_numbers": False},
        ).shift(DOWN * 0.3)

        x_lbl = Text("Mode number", font_size=22
                     ).next_to(ax.x_axis, DOWN, buff=0.35)
        y_lbl = Text("Normalised eigenvalue", font_size=22
                     ).rotate(PI / 2).next_to(ax.y_axis, LEFT, buff=0.35)
        x_nums = VGroup(*[
            Text(str(n), font_size=16
                 ).next_to(ax.c2p(n, 0), DOWN, buff=0.12)
            for n in [1, 10, 20, 30, 40, 50]
        ])

        self.play(Create(ax), FadeIn(x_lbl), FadeIn(y_lbl),
                  FadeIn(x_nums), run_time=1)

        # ---- eigenvalue bars — real data with winter colourmap ----
        bars = VGroup()
        bar_w = 9.5 / (N_SHOW + 2) * 0.7
        for i in range(N_SHOW):
            bot = ax.c2p(i + 1, 0)
            top = ax.c2p(i + 1, eig_norm[i])
            h   = abs(top[1] - bot[1])
            bar = Rectangle(
                width=bar_w, height=max(h, 0.01),
                fill_opacity=0.85,
                stroke_width=0.3, stroke_color=TEXT,
            )
            # Winter colour gradient: blue (low mode) → teal → seafoam (high)
            t = i / max(N_SHOW - 1, 1)
            if t < 0.5:
                col = interpolate_color(ManimColor(SIGNAL), ManimColor(TEAL), t * 2)
            else:
                col = interpolate_color(ManimColor(TEAL), ManimColor(SEAFOAM), (t - 0.5) * 2)
            bar.set_fill(col)
            bar.move_to(bot, DOWN)
            bars.add(bar)

        self.play(
            LaggedStart(*[GrowFromEdge(b, DOWN) for b in bars],
                        lag_ratio=0.03),
            run_time=2,
        )
        self.wait(1)

        # annotation: dominant modes vs noise
        brace_top = Brace(
            VGroup(bars[0], bars[4]), UP, color=SIGNAL)
        top_lbl = Text("Dominant wave modes", font_size=18,
                       weight="BOLD", color=SIGNAL
                       ).next_to(brace_top, UP, buff=0.08)

        brace_tail = Brace(
            VGroup(bars[25], bars[N_SHOW - 1]), UP, color=NOISE)
        tail_lbl = Text("Noise / spurious vectors", font_size=18,
                        color=NOISE).next_to(brace_tail, UP, buff=0.08)

        self.play(
            GrowFromCenter(brace_top), FadeIn(top_lbl, shift=DOWN * 0.1),
            GrowFromCenter(brace_tail), FadeIn(tail_lbl, shift=DOWN * 0.1),
            run_time=0.8,
        )
        self.wait(3)

        # ---- cumulative energy overlay (right y-axis) ----
        cum_pts = [ax.c2p(i + 1, cum[i]) for i in range(N_SHOW)]
        cum_line = VMobject(color=ENERGY, stroke_width=3)
        cum_line.set_points_smoothly(cum_pts)

        cum_lbl = Text("Cumulative energy", font_size=17, color=ENERGY
                       ).next_to(cum_pts[-1], UR, buff=0.12)

        pct_labels = VGroup()
        for pct in [0.10, 0.20, 0.30, 0.40]:
            idx = int(np.searchsorted(cum, pct))
            if idx < N_SHOW:
                pct_labels.add(
                    Text(f"{int(pct*100)}%", font_size=14, color=ENERGY
                         ).next_to(ax.c2p(idx + 1, cum[idx]), UR, buff=0.06)
                )

        self.play(Create(cum_line), FadeIn(cum_lbl),
                  FadeIn(pct_labels, lag_ratio=0.1),
                  run_time=1.5)
        self.wait(2)

        # ---- data footnote ----
        footnote = Text(
            f"L8 G12 liquid phase  ·  {TOTAL_MODES:,} total modes  "
            f"·  {N_MODES_99} modes for 99% energy",
            font_size=15, color=MUTED,
        ).to_edge(DOWN, buff=0.3)

        self.play(FadeIn(footnote, shift=UP * 0.1), run_time=0.6)

        # ---- closing statement ----
        conclusion = Text(
            "Reconstruct with dominant modes  →  filter out noise",
            font_size=28, weight="BOLD", color=TEXT,
        ).next_to(footnote, UP, buff=0.3)

        self.play(FadeIn(conclusion, shift=UP * 0.2), run_time=1)
        self.wait(3)

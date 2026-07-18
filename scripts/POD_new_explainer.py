"""
POD Explainer Animation for MMME4086 Presentation (Slide 4)
Render:  manim -qh pod_explainer.py PODExplainer
         (-qh = 1080p, use -qk for 4K, -ql for quick preview)
Output:  media/videos/pod_explainer/1080p60/PODExplainer.mp4
"""

from manim import *
import numpy as np

# ---------- colour palette (professional, presentation-safe) ----------
BG_COLOR = "#1e1e2e"
ACCENT_BLUE = "#89b4fa"
ACCENT_GREEN = "#a6e3a1"
ACCENT_RED = "#f38ba8"
ACCENT_YELLOW = "#f9e2af"
ACCENT_MAUVE = "#cba6f7"
TEXT_COLOR = "#cdd6f4"
DIM_COLOR = "#585b70"


class PODExplainer(Scene):
    def construct(self):
        self.camera.background_color = BG_COLOR

        self.scene_title()
        self.scene_snapshots()
        self.scene_build_matrix()
        self.scene_decompose()
        self.scene_energy_ranking()
        self.scene_filter()
        self.scene_summary()

    # ------------------------------------------------------------------ #
    #  SCENE 1 — Title card                                               #
    # ------------------------------------------------------------------ #
    def scene_title(self):
        title = Text("Proper Orthogonal Decomposition", font_size=44, color=TEXT_COLOR)
        subtitle = Text("(POD)", font_size=36, color=ACCENT_BLUE)
        aka = Text("a.k.a. PCA for flow fields", font_size=24, color=DIM_COLOR)
        group = VGroup(title, subtitle, aka).arrange(DOWN, buff=0.35)

        self.play(Write(title), run_time=1.2)
        self.play(FadeIn(subtitle, shift=UP * 0.2), run_time=0.6)
        self.play(FadeIn(aka, shift=UP * 0.15), run_time=0.5)
        self.wait(1.5)
        self.play(FadeOut(group), run_time=0.6)

    # ------------------------------------------------------------------ #
    #  SCENE 2 — 2 500 noisy snapshots fly in                            #
    # ------------------------------------------------------------------ #
    def scene_snapshots(self):
        label = Text("2,500 PIV snapshots", font_size=30, color=TEXT_COLOR).to_edge(UP, buff=0.5)
        self.play(FadeIn(label))

        np.random.seed(42)

        # Create a grid of mini "velocity field" thumbnails
        def make_thumbnail(noisy=True):
            arrows = VGroup()
            for r in range(4):
                for c in range(5):
                    angle = -0.3 + 0.15 * r
                    if noisy:
                        angle += np.random.normal(0, 0.6)
                    length = 0.18 if not noisy else 0.18 * (0.5 + np.random.random())
                    a = Arrow(
                        start=ORIGIN, end=length * np.array([np.cos(angle), np.sin(angle), 0]),
                        buff=0, stroke_width=1.5,
                        color=interpolate_color(
                            ManimColor(ACCENT_BLUE), ManimColor(ACCENT_RED),
                            (angle + 1) / 2
                        ),
                        max_tip_length_to_length_ratio=0.25,
                    )
                    a.move_to(np.array([c * 0.28 - 0.56, r * 0.28 - 0.42, 0]))
                    arrows.add(a)
            box = SurroundingRectangle(arrows, buff=0.08, color=DIM_COLOR, stroke_width=1)
            return VGroup(box, arrows)

        # Show a few representative frames cascading in
        frames = VGroup()
        positions = [LEFT * 4 + DOWN * 0.3, LEFT * 1.5 + DOWN * 0.3,
                     RIGHT * 1 + DOWN * 0.3, RIGHT * 3.5 + DOWN * 0.3]
        for i, pos in enumerate(positions):
            f = make_thumbnail(noisy=True).scale(0.9).move_to(pos)
            frames.add(f)

        dots = Text("...", font_size=48, color=DIM_COLOR).move_to(RIGHT * 5.5 + DOWN * 0.3)

        self.play(
            LaggedStart(*[FadeIn(f, shift=DOWN * 0.3) for f in frames],
                        lag_ratio=0.25),
            run_time=1.8
        )
        self.play(FadeIn(dots))

        noise_label = Text("noisy vectors everywhere!", font_size=22, color=ACCENT_RED)
        noise_label.next_to(frames, DOWN, buff=0.5)
        self.play(FadeIn(noise_label, shift=UP * 0.15))
        self.wait(1.5)

        self.snapshot_frames = frames
        self.play(FadeOut(dots), FadeOut(noise_label), FadeOut(label))

    # ------------------------------------------------------------------ #
    #  SCENE 3 — Stack snapshots into a matrix                           #
    # ------------------------------------------------------------------ #
    def scene_build_matrix(self):
        heading = Text("Stack into a snapshot matrix", font_size=28, color=TEXT_COLOR).to_edge(UP, buff=0.5)
        self.play(FadeIn(heading))

        # Animate frames collapsing into thin column strips
        cols = VGroup()
        target_x_start = -2.0
        for i, f in enumerate(self.snapshot_frames):
            col = Rectangle(width=0.15, height=2.5, fill_opacity=0.6,
                            fill_color=ACCENT_BLUE, stroke_width=0.8, stroke_color=DIM_COLOR)
            col.move_to(np.array([target_x_start + i * 0.22, -0.2, 0]))
            cols.add(col)

        # Extra columns to suggest 2500
        for j in range(12):
            col = Rectangle(width=0.15, height=2.5, fill_opacity=max(0.6 - j * 0.04, 0.15),
                            fill_color=ACCENT_BLUE, stroke_width=0.8, stroke_color=DIM_COLOR)
            col.move_to(np.array([target_x_start + (4 + j) * 0.22, -0.2, 0]))
            cols.add(col)

        ellipsis_dots = Text("···", font_size=40, color=DIM_COLOR)
        ellipsis_dots.next_to(cols, RIGHT, buff=0.15)

        bracket_l = MathTex(r"\Bigg[", font_size=72, color=TEXT_COLOR).next_to(cols, LEFT, buff=0.15)
        bracket_r = MathTex(r"\Bigg]", font_size=72, color=TEXT_COLOR).next_to(ellipsis_dots, RIGHT, buff=0.15)

        # Labels
        col_label = MathTex(r"\mathbf{u}_1 \;\; \mathbf{u}_2 \;\; \mathbf{u}_3 \;\; \cdots \;\; \mathbf{u}_{2500}",
                            font_size=22, color=DIM_COLOR)
        col_label.next_to(bracket_r, RIGHT, buff=0.3).align_to(cols, DOWN).shift(DOWN*0.4)

        row_label = Text("Vx & Vy\nconcatenated", font_size=18, color=ACCENT_YELLOW)
        row_label.next_to(bracket_l, LEFT, buff=0.25)

        self.play(
            *[ReplacementTransform(f, cols[i]) for i, f in enumerate(self.snapshot_frames)],
            run_time=1.5
        )
        self.play(
            LaggedStart(*[FadeIn(c, shift=RIGHT * 0.1) for c in cols[4:]], lag_ratio=0.06),
            run_time=1.0
        )
        self.play(FadeIn(ellipsis_dots), FadeIn(bracket_l), FadeIn(bracket_r))
        self.play(FadeIn(col_label), FadeIn(row_label))
        self.wait(1.5)

        self.play(
            *[FadeOut(m) for m in [heading, cols, ellipsis_dots, bracket_l,
                                    bracket_r, col_label, row_label]],
            run_time=0.8
        )

    # ------------------------------------------------------------------ #
    #  SCENE 4 — Decomposition equation + PCA analogy                    #
    # ------------------------------------------------------------------ #
    def scene_decompose(self):
        heading = Text("Decompose into modes", font_size=28, color=TEXT_COLOR).to_edge(UP, buff=0.5)

        eq = MathTex(
            r"\mathbf{u}(x,t)", r"=", r"\bar{\mathbf{u}}(x)", r"+",
            r"\sum_{k=1}^{N}", r"a_k(t)", r"\boldsymbol{\phi}_k(x)",
            font_size=38, color=TEXT_COLOR
        )
        eq.move_to(UP * 1.2)

        # Colour-code the parts
        eq[0].set_color(ACCENT_RED)       # original field
        eq[2].set_color(DIM_COLOR)        # mean
        eq[5].set_color(ACCENT_YELLOW)    # temporal coefficients
        eq[6].set_color(ACCENT_GREEN)     # spatial modes

        labels = VGroup(
            Text("original\nsnapshot", font_size=16, color=ACCENT_RED),
            Text("time-\naverage", font_size=16, color=DIM_COLOR),
            Text("temporal\ncoefficients", font_size=16, color=ACCENT_YELLOW),
            Text("spatial\nmodes", font_size=16, color=ACCENT_GREEN),
        )
        labels[0].next_to(eq[0], DOWN, buff=0.5)
        labels[1].next_to(eq[2], DOWN, buff=0.5)
        labels[2].next_to(eq[5], DOWN, buff=0.5)
        labels[3].next_to(eq[6], DOWN, buff=0.5)

        arrows_to_eq = VGroup(*[
            Arrow(l.get_top(), target.get_bottom(), buff=0.08,
                  stroke_width=2, color=l.color, max_tip_length_to_length_ratio=0.2)
            for l, target in zip(labels, [eq[0], eq[2], eq[5], eq[6]])
        ])

        # PCA scatter-plot analogy
        analogy_title = Text("Same idea as PCA:", font_size=22, color=ACCENT_MAUVE).shift(DOWN * 1.8 + LEFT * 3)

        # Simple 2D scatter with rotated axes
        np.random.seed(7)
        pts_raw = np.random.multivariate_normal([0, 0], [[2, 1.5], [1.5, 2]], 60)
        dots = VGroup(*[
            Dot(point=np.array([p[0] * 0.3, p[1] * 0.3, 0]), radius=0.03,
                color=ACCENT_BLUE, fill_opacity=0.5)
            for p in pts_raw
        ]).shift(DOWN * 2.2 + RIGHT * 1.5)

        # Principal axes
        angle = np.pi / 4
        pc1 = Arrow(ORIGIN, 1.8 * np.array([np.cos(angle), np.sin(angle), 0]),
                     buff=0, color=ACCENT_GREEN, stroke_width=3)
        pc2 = Arrow(ORIGIN, 1.0 * np.array([-np.sin(angle), np.cos(angle), 0]),
                     buff=0, color=ACCENT_YELLOW, stroke_width=3)
        axes_group = VGroup(pc1, pc2).move_to(dots.get_center())

        pc1_label = Text("Mode 1\n(most energy)", font_size=14, color=ACCENT_GREEN)
        pc1_label.next_to(pc1.get_end(), UR, buff=0.1)
        pc2_label = Text("Mode 2", font_size=14, color=ACCENT_YELLOW)
        pc2_label.next_to(pc2.get_end(), UL, buff=0.1)

        self.play(FadeIn(heading), Write(eq), run_time=1.5)
        self.play(
            LaggedStart(*[FadeIn(l, shift=UP * 0.1) for l in labels], lag_ratio=0.2),
            LaggedStart(*[GrowArrow(a) for a in arrows_to_eq], lag_ratio=0.2),
            run_time=1.5
        )
        self.wait(1)
        self.play(FadeIn(analogy_title), FadeIn(dots, lag_ratio=0.02), run_time=1)
        self.play(GrowArrow(pc1), GrowArrow(pc2), run_time=1)
        self.play(FadeIn(pc1_label), FadeIn(pc2_label))
        self.wait(2)

        self.play(
            *[FadeOut(m) for m in [heading, eq, labels, arrows_to_eq,
                                    analogy_title, dots, axes_group,
                                    pc1_label, pc2_label]],
            run_time=0.8
        )

    # ------------------------------------------------------------------ #
    #  SCENE 5 — Energy ranking bar chart                                #
    # ------------------------------------------------------------------ #
    def scene_energy_ranking(self):
        heading = Text("Modes ranked by energy", font_size=28, color=TEXT_COLOR).to_edge(UP, buff=0.5)

        # Simulated eigenvalue spectrum (exponential-ish decay)
        n_bars = 15
        energies = [40, 18, 10, 6, 4, 3, 2.2, 1.8, 1.5, 1.2, 1, 0.8, 0.6, 0.5, 0.4]
        max_h = 3.5
        bar_w = 0.35
        bars = VGroup()
        bar_labels = VGroup()
        x_start = -3.5

        for i, e in enumerate(energies):
            h = e / energies[0] * max_h
            bar = Rectangle(
                width=bar_w, height=h,
                fill_opacity=0.8, stroke_width=0.5,
                fill_color=ACCENT_GREEN if i < 5 else (ACCENT_YELLOW if i < 10 else ACCENT_RED),
                stroke_color=DIM_COLOR,
            )
            bar.move_to(np.array([x_start + i * (bar_w + 0.12), -1.8 + h / 2, 0]))
            bars.add(bar)

            lbl = Text(str(i + 1), font_size=12, color=DIM_COLOR)
            lbl.next_to(bar, DOWN, buff=0.08)
            bar_labels.add(lbl)

        x_label = Text("Mode number →", font_size=18, color=DIM_COLOR).next_to(bar_labels, DOWN, buff=0.25)

        brace_phys = Brace(bars[:5], UP, color=ACCENT_GREEN)
        brace_phys_label = Text("real physics", font_size=18, color=ACCENT_GREEN)
        brace_phys_label.next_to(brace_phys, UP, buff=0.1)

        brace_noise = Brace(bars[10:], UP, color=ACCENT_RED)
        brace_noise_label = Text("noise", font_size=18, color=ACCENT_RED)
        brace_noise_label.next_to(brace_noise, UP, buff=0.1)

        self.play(FadeIn(heading))
        self.play(
            LaggedStart(*[GrowFromEdge(b, DOWN) for b in bars], lag_ratio=0.1),
            LaggedStart(*[FadeIn(l) for l in bar_labels], lag_ratio=0.1),
            run_time=2.0
        )
        self.play(FadeIn(x_label))
        self.play(
            GrowFromCenter(brace_phys), FadeIn(brace_phys_label),
            GrowFromCenter(brace_noise), FadeIn(brace_noise_label),
            run_time=1.2
        )
        self.wait(2)

        self.play(
            *[FadeOut(m) for m in [heading, bars, bar_labels, x_label,
                                    brace_phys, brace_phys_label,
                                    brace_noise, brace_noise_label]],
            run_time=0.8
        )

    # ------------------------------------------------------------------ #
    #  SCENE 6 — Filtering: keep signal, discard noise                   #
    # ------------------------------------------------------------------ #
    def scene_filter(self):
        heading = Text("Reconstruct → filter", font_size=28, color=TEXT_COLOR).to_edge(UP, buff=0.5)

        np.random.seed(99)

        def make_field(noisy, label_text, color):
            arrows = VGroup()
            for r in range(6):
                for c in range(8):
                    base_angle = -0.4 + 0.2 * np.sin(c * 0.8) + 0.1 * r
                    if noisy:
                        base_angle += np.random.normal(0, 0.7)
                    length = 0.22 if not noisy else 0.22 * (0.4 + 0.6 * np.random.random())
                    a = Arrow(
                        start=ORIGIN,
                        end=length * np.array([np.cos(base_angle), np.sin(base_angle), 0]),
                        buff=0, stroke_width=1.8,
                        color=interpolate_color(ManimColor(ACCENT_BLUE), ManimColor(color),
                                                abs(base_angle) / 1.5),
                        max_tip_length_to_length_ratio=0.25,
                    )
                    a.move_to(np.array([c * 0.32 - 1.12, r * 0.32 - 0.8, 0]))
                    arrows.add(a)
            box = SurroundingRectangle(arrows, buff=0.12, color=DIM_COLOR, stroke_width=1.5)
            label = Text(label_text, font_size=20, color=color)
            label.next_to(box, DOWN, buff=0.2)
            return VGroup(box, arrows, label)

        raw = make_field(noisy=True, label_text="Raw PIV snapshot", color=ACCENT_RED)
        raw.shift(LEFT * 3.2)

        clean = make_field(noisy=False, label_text="POD-filtered (90% energy)", color=ACCENT_GREEN)
        clean.shift(RIGHT * 3.2)

        arrow_between = Arrow(LEFT * 0.6, RIGHT * 0.6, buff=0, color=ACCENT_YELLOW, stroke_width=4)
        arrow_label = Text("keep top\nmodes only", font_size=16, color=ACCENT_YELLOW)
        arrow_label.next_to(arrow_between, UP, buff=0.15)

        self.play(FadeIn(heading))
        self.play(FadeIn(raw, shift=LEFT * 0.3), run_time=1)
        self.play(GrowArrow(arrow_between), FadeIn(arrow_label))
        self.play(FadeIn(clean, shift=RIGHT * 0.3), run_time=1)
        self.wait(2)

        self.play(
            *[FadeOut(m) for m in [heading, raw, clean, arrow_between, arrow_label]],
            run_time=0.8
        )

    # ------------------------------------------------------------------ #
    #  SCENE 7 — Closing summary                                         #
    # ------------------------------------------------------------------ #
    def scene_summary(self):
        points = VGroup(
            Text("✓  Decomposes flow into energy-ranked modes", font_size=24, color=ACCENT_GREEN),
            Text("✓  High-energy modes = real wave structures", font_size=24, color=ACCENT_BLUE),
            Text("✓  Low-energy modes = noise → discard", font_size=24, color=ACCENT_RED),
            Text("✓  Vx and Vy must be concatenated together", font_size=24, color=ACCENT_YELLOW),
            Text("✓  Water phase only — converges faster", font_size=24, color=ACCENT_MAUVE),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.35).move_to(ORIGIN)

        self.play(
            LaggedStart(*[FadeIn(p, shift=RIGHT * 0.3) for p in points], lag_ratio=0.4),
            run_time=3
        )
        self.wait(2.5)
        self.play(FadeOut(points), run_time=0.8)

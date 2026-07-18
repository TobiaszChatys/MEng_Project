from manim import *
import numpy as np

# --- ROSÉ PINE DAWN PALETTE ---
RP_BASE = "#faf4ed"
RP_TEXT = "#575279"
RP_LOVE = "#b4637a"
RP_GOLD = "#ea9d34"
RP_ROSE = "#d7827e"
RP_PINE = "#286983"
RP_FOAM = "#56949f"
RP_IRIS = "#907aa9"
RP_MUTED = "#9893a5"

Tex.set_default(color=RP_TEXT)
MathTex.set_default(color=RP_TEXT)
Text.set_default(color=RP_TEXT)

class IncrementalPOD(Scene):
    def construct(self):
        self.camera.background_color = RP_BASE

        self.scene_1_normalization()
        self.clear()
        
        self.scene_2_incremental_blocks()
        self.clear()
        
        self.scene_3_energy_threshold()
        self.clear()
        
        self.scene_4_hilbert_pairs()

    # ---------------------------------------------------------
    # SCENE 1: Normalisation (X - mean = Fluctuations)
    # ---------------------------------------------------------
    def scene_1_normalization(self):
        title = Title("1. Normalisation (Fluctuations)", color=RP_TEXT)
        self.add(title)

        # Represent the matrices
        eq = MathTex(
            "X", "-", "\\bar{X}", "=", "X'"
        ).scale(1.5)
        
        # Color specific parts
        eq[0].set_color(RP_PINE)  # Raw Data
        eq[2].set_color(RP_GOLD)  # Mean
        eq[4].set_color(RP_LOVE)  # Fluctuations

        desc_raw = Text("Raw Snapshot Matrix", font_size=20).next_to(eq[0], UP)
        desc_mean = Text("Mean Matrix", font_size=20).next_to(eq[2], DOWN)
        desc_fluc = Text("Fluctuation Matrix", font_size=20).next_to(eq[4], UP)

        self.play(Write(eq[0]), FadeIn(desc_raw))
        self.wait(0.5)
        self.play(Write(eq[1:3]), FadeIn(desc_mean))
        self.wait(0.5)
        self.play(Write(eq[3:]), FadeIn(desc_fluc))
        self.wait(2)

    # ---------------------------------------------------------
    # SCENE 2: The Incremental Algorithm (Nested Loops)
    # ---------------------------------------------------------
    def scene_2_incremental_blocks(self):
        title = Title("2. Incremental Covariance (Block Math)", color=RP_TEXT)
        self.add(title)

        # Show the large fluctuation matrix being split into blocks
        fluc_matrix = Rectangle(width=6, height=3, color=RP_LOVE, fill_opacity=0.2).shift(LEFT * 3.5 + UP * 1)
        fluc_label = MathTex("X'").move_to(fluc_matrix)
        
        self.play(Create(fluc_matrix), Write(fluc_label))
        self.wait(1)

        # Slice it into 3 blocks (representing block_size = 350)
        b1 = Rectangle(width=2, height=3, color=RP_LOVE, fill_opacity=0.5).move_to(fluc_matrix.get_left() + RIGHT * 1)
        b2 = Rectangle(width=2, height=3, color=RP_PINE, fill_opacity=0.5).move_to(b1.get_right() + RIGHT * 1)
        b3 = Rectangle(width=2, height=3, color=RP_GOLD, fill_opacity=0.5).move_to(b2.get_right() + RIGHT * 1)

        b1_label = MathTex("B_1").move_to(b1)
        b2_label = MathTex("B_2").move_to(b2)
        b3_label = MathTex("B_3").move_to(b3)

        self.play(
            FadeOut(fluc_label),
            Transform(fluc_matrix.copy(), VGroup(b1, b2, b3)),
            Write(VGroup(b1_label, b2_label, b3_label))
        )
        self.wait(1)

        # Set up the Covariance Matrix Grid (C)
        cov_matrix = Square(side_length=3, color=RP_TEXT).shift(RIGHT * 3.5 + DOWN * 0.5)
        # Fixed the double superscript error by using {X'}^T
        cov_label = MathTex(r"C = \frac{1}{N} {X'}^T X'").next_to(cov_matrix, UP)
        # Grid lines for the 3x3 block sections in C
        grid = VGroup(
            Line(cov_matrix.get_top() + DOWN*1, cov_matrix.get_bottom() + UP*2, color=RP_MUTED),
            Line(cov_matrix.get_top() + DOWN*2, cov_matrix.get_bottom() + UP*1, color=RP_MUTED),
            Line(cov_matrix.get_left() + RIGHT*1, cov_matrix.get_right() + LEFT*2, color=RP_MUTED),
            Line(cov_matrix.get_left() + RIGHT*2, cov_matrix.get_right() + LEFT*1, color=RP_MUTED),
        )

        self.play(Create(cov_matrix), Create(grid), Write(cov_label))

        # --- Animate the nested loop computation ---
        math_eq = MathTex("C_{1,2} = B_1^T \\times B_2", color=RP_TEXT).shift(DOWN * 3)
        
        # 1. Take Block 1
        target_b1 = b1.copy()
        self.play(target_b1.animate.shift(DOWN * 3).rotate(-PI/2).set_fill(opacity=0.8)) # Transpose it
        
        # 2. Take Block 2
        target_b2 = b2.copy()
        self.play(target_b2.animate.next_to(target_b1, RIGHT, buff=0.5))
        
        self.play(Write(math_eq))
        self.wait(1)

        # 3. Multiply and fill the grid square [1, 2] in C
        computed_block = Rectangle(width=1, height=1, color=RP_PINE, fill_opacity=0.8)
        # Position it in row 1, col 2 of the C matrix
        computed_block.move_to(cov_matrix.get_corner(UL) + RIGHT * 1.5 + DOWN * 0.5)

        self.play(
            ReplacementTransform(VGroup(target_b1, target_b2, math_eq), computed_block)
        )
        
        # Show symmetry (if block ~= upper_block: C[j,i] = computed')
        sym_block = Rectangle(width=1, height=1, color=RP_PINE, fill_opacity=0.4)
        sym_block.move_to(cov_matrix.get_corner(UL) + RIGHT * 0.5 + DOWN * 1.5)
        
        sym_text = Text("Symmetric Transpose", font_size=16, color=RP_MUTED).next_to(cov_matrix, DOWN)
        self.play(TransformFromCopy(computed_block, sym_block), Write(sym_text))
        self.wait(2)

    # ---------------------------------------------------------
    # SCENE 3: Eigen-Decomposition & 99% Threshold
    # ---------------------------------------------------------
    def scene_3_energy_threshold(self):
        title = Title(r"3. Eigenvalues \& Energy Threshold", color=RP_TEXT)
        self.add(title)

        # Plotting the Cumulative Energy Graph
        ax = Axes(
            x_range=[0, 50, 10], y_range=[0, 1.1, 0.2], 
            x_length=7, y_length=5, 
            axis_config={"color": RP_TEXT}
        ).shift(DOWN * 0.5)
        
        x_label = ax.get_x_axis_label("Modes")
        y_label = ax.get_y_axis_label("Cumulative Energy")

        self.play(Create(ax), Write(x_label), Write(y_label))

        # Generate cumulative energy curve (looks like 1 - e^(-x))
        curve = ax.plot(lambda x: 1 - np.exp(-0.15 * x), x_range=[0, 50], color=RP_LOVE, stroke_width=4)
        
        self.play(Create(curve), run_time=2)

        # Draw the 99% threshold line
        threshold_line = ax.get_horizontal_line(ax.c2p(50, 0.99), color=RP_PINE, line_func=DashedLine)
        threshold_label = MathTex("99\\% \\text{ Energy}", font_size=24, color=RP_PINE).next_to(threshold_line, UP, aligned_edge=LEFT)

        self.play(Create(threshold_line), Write(threshold_label))
        
        # Find intersection
        dot = Dot(ax.c2p(30.7, 0.99), color=RP_GOLD)
        v_line = ax.get_vertical_line(ax.c2p(30.7, 0.99), color=RP_GOLD, line_func=DashedLine)
        modes_label = MathTex("K \\approx 31 \\text{ Modes}", font_size=24, color=RP_GOLD).next_to(v_line, DOWN)

        self.play(Create(dot), Create(v_line), Write(modes_label))
        self.wait(2)

    # ---------------------------------------------------------
    # SCENE 4: Hilbert Transform & Mode Pairing
    # ---------------------------------------------------------
    def scene_4_hilbert_pairs(self):
        title = Title("4. Mode Pairing via Hilbert Transform", color=RP_TEXT)
        self.add(title)

        ax = Axes(
            x_range=[0, 4*PI, PI/2], y_range=[-1.5, 1.5, 1], 
            x_length=8, y_length=4, 
            axis_config={"color": RP_TEXT}
        ).shift(DOWN * 1)

        self.play(Create(ax))

        # Plot two out-of-phase modes
        mode_a = ax.plot(lambda x: np.sin(x), color=RP_LOVE)
        mode_b = ax.plot(lambda x: np.sin(x - PI/4), color=RP_GOLD) # Shifted by 45 degrees

        label_a = Text("Mode A (Eigenvector)", font_size=20, color=RP_LOVE).next_to(mode_a, UP).shift(LEFT*2)
        label_b = Text("Mode B (Eigenvector)", font_size=20, color=RP_GOLD).next_to(mode_b, UP).shift(RIGHT*2)

        self.play(Create(mode_a), Write(label_a))
        self.play(Create(mode_b), Write(label_b))
        self.wait(1)
        # We use r"..." (raw strings) and \big( instead of \left( to fix the LaTeX error
        hilbert_eq = MathTex(
            r"\Delta \phi", r"=", r"\text{median}\big(", r"\text{angle}\big(", r"\frac{\mathcal{H}(A)}{\mathcal{H}(B)}", r"\big)\big)"
        ).scale(0.8).next_to(title, DOWN)
        
        hilbert_eq[0].set_color(RP_PINE) # Delta phi
        # Show the math from your script
        self.play(Write(hilbert_eq))

        # Show the phase shift visually
        arrow = DoubleArrow(
            start=ax.c2p(PI/2, 1), 
            end=ax.c2p(PI/2 + PI/4, 1), 
            color=RP_PINE, buff=0
        )
        phase_text = MathTex("45^{\\circ} \\text{ Phase Shift}", font_size=24, color=RP_PINE).next_to(arrow, UP)

        self.play(GrowArrow(arrow), Write(phase_text))
        
        conclusion = Text("Conclusion: Modes A & B form a traveling wave pair.", font_size=24, color=RP_TEXT).next_to(ax, DOWN)
        self.play(Write(conclusion))
        self.wait(3)

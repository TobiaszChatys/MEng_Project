from manim import *
import numpy as np

class PODExplainer(Scene):
    def construct(self):
        # Set a background color standard for 3b1b style
        self.camera.background_color = "#141414"

        self.scene_1_chaos()
        self.clear()
        
        self.scene_2_snapshots()
        self.clear()
        
        self.scene_3_energy_orthogonality()
        self.clear()
        
        self.scene_4_math()
        self.clear()
        
        self.scene_5_reconstruction()

    # ---------------------------------------------------------
    # SCENE 1: The Problem (Chaos in Fluid Dynamics)
    # ---------------------------------------------------------
    def scene_1_chaos(self):
        title = Title("1. The Problem: Drowning in Data")
        self.play(Write(title))

        # Create a chaotic-looking vector field simulating gas-sheared film
        def fluid_func(pos):
            x, y = pos[:2]
            u = np.sin(y) + np.cos(x * 0.5) * np.sin(y * 2)
            v = np.cos(x) - np.sin(x * 0.5)
            return np.array([u, v, 0])
        
        vector_field = ArrowVectorField(
            fluid_func, 
            x_range=[-6, 6, 0.75], 
            y_range=[-3, 3, 0.75],
            colors=[BLUE_E, BLUE_C, TEAL_A] # Fixed!
       )
        
        # Animate the vectors appearing as if flowing
        self.play(Create(vector_field), run_time=3)
        self.wait(2)

        # "The fluid freezes. A grid overlays the fluid."
        grid = NumberPlane(
            x_range=[-7, 7, 1], 
            y_range=[-4, 4, 1],
            background_line_style={"stroke_opacity": 0.3}
        )
        self.play(Create(grid), vector_field.animate.set_opacity(0.5))
        self.wait(2)

    # ---------------------------------------------------------
    # SCENE 2: The Method of Snapshots (Data Matrix)
    # ---------------------------------------------------------
    def scene_2_snapshots(self):
        title = Title("2. The Method of Snapshots")
        self.add(title)

        # Represent a PIV grid flattening into a column vector
        grid_matrix = Matrix([
            ["u_{11}", "u_{12}"],
            ["u_{21}", "u_{22}"]
        ]).shift(LEFT * 4)
        
        arrow = Arrow(LEFT * 2, RIGHT * 2)
        
        col_vector = Matrix([
            ["u_{11}"],
            ["u_{12}"],
            ["u_{21}"],
            ["u_{22}"]
        ]).shift(RIGHT * 4)

        # Flattening animation
        self.play(Write(grid_matrix))
        self.play(GrowArrow(arrow))
        self.play(TransformFromCopy(grid_matrix, col_vector))
        self.wait(2)

        self.play(FadeOut(grid_matrix), FadeOut(arrow), col_vector.animate.shift(LEFT * 4))

        # Build the large Data Matrix X
        dots = MathTex("\\dots").next_to(col_vector, RIGHT)
        
        big_matrix = Matrix([
            ["u_1(t_1)", "u_1(t_2)", "\\dots"],
            ["u_2(t_1)", "u_2(t_2)", "\\dots"],
            ["\\vdots", "\\vdots", "\\ddots"],
            ["u_n(t_1)", "u_n(t_2)", "\\dots"]
        ]).shift(RIGHT * 2)
        
        matrix_name = MathTex("X = ").next_to(big_matrix, LEFT)

        self.play(Write(dots))
        self.play(
            ReplacementTransform(VGroup(col_vector, dots), big_matrix),
            Write(matrix_name)
        )
        self.wait(2)

    # ---------------------------------------------------------
    # SCENE 3: The Core Idea (Energy & Orthogonality)
    # ---------------------------------------------------------
    def scene_3_energy_orthogonality(self):
        title = Title("3. Finding the Best Perspective")
        self.add(title)

        # Create a stretched, rotated cloud of data points (scatter plot)
        axes = Axes(x_range=[-4, 4], y_range=[-3, 3], x_length=8, y_length=6)
        
        # Generate points
        np.random.seed(42)
        points = []
        for _ in range(150):
            x = np.random.normal(0, 2)
            y = np.random.normal(0, 0.4) # Compressed Y to make it an ellipse
            # Rotate points by 30 degrees (PI/6)
            rotated_x = x * np.cos(PI/6) - y * np.sin(PI/6)
            rotated_y = x * np.sin(PI/6) + y * np.cos(PI/6)
            points.append([rotated_x, rotated_y, 0])
            
        dots = VGroup(*[Dot(axes.c2p(p[0], p[1]), color=BLUE_C, radius=0.05) for p in points])

        self.play(Create(axes))
        self.play(FadeIn(dots, lag_ratio=0.05), run_time=2)
        self.wait(1)

        # Draw the Principal Axes (Orthogonal Modes)
        origin = axes.c2p(0, 0)
        mode1_end = axes.c2p(3 * np.cos(PI/6), 3 * np.sin(PI/6))
        mode2_end = axes.c2p(-1 * np.sin(PI/6), 1 * np.cos(PI/6))

        mode1_vector = Arrow(origin, mode1_end, buff=0, color=YELLOW, stroke_width=6, max_tip_length_to_length_ratio=0.1)
        mode2_vector = Arrow(origin, mode2_end, buff=0, color=RED, stroke_width=6, max_tip_length_to_length_ratio=0.2)

        mode1_label = Text("Mode 1 (Max Energy)", color=YELLOW, font_size=24).next_to(mode1_end, RIGHT)
        mode2_label = Text("Mode 2", color=RED, font_size=24).next_to(mode2_end, UP)

        self.play(GrowArrow(mode1_vector), Write(mode1_label))
        self.wait(1)
        self.play(GrowArrow(mode2_vector), Write(mode2_label))
        self.wait(2)

    # ---------------------------------------------------------
    # SCENE 4: The Math (SVD / Eigen-decomposition)
    # ---------------------------------------------------------
    def scene_4_math(self):
        title = Title("4. The Math (Eigen-decomposition)")
        self.add(title)

        # Equation 1: Covariance
        eq1 = MathTex("C", "=", "X^T", "X", font_size=60)
        eq1_desc = Text("Covariance Matrix (Spatial Correlation)", font_size=24).next_to(eq1, DOWN)
        
        self.play(Write(eq1))
        self.play(FadeIn(eq1_desc))
        self.wait(2)
        self.play(eq1.animate.shift(UP * 2), eq1_desc.animate.shift(UP * 2))

        # Equation 2: Eigenvalue problem
        eq2 = MathTex("C", "\\mathbf{v}", "=", "\\lambda", "\\mathbf{v}", font_size=60)
        
        self.play(Write(eq2))
        self.wait(1)

        # Highlight components
        box_v = SurroundingRectangle(eq2[1], color=BLUE)
        box_v2 = SurroundingRectangle(eq2[4], color=BLUE)
        label_v = Text("Spatial Modes", color=BLUE, font_size=24).next_to(box_v, DOWN)

        box_lambda = SurroundingRectangle(eq2[3], color=YELLOW)
        label_lambda = Text("Energy", color=YELLOW, font_size=24).next_to(box_lambda, UP)

        self.play(Create(box_v), Create(box_v2), Write(label_v))
        self.play(Create(box_lambda), Write(label_lambda))
        self.wait(3)

    # ---------------------------------------------------------
    # SCENE 5: Reconstruction (Order from Chaos)
    # ---------------------------------------------------------
    def scene_5_reconstruction(self):
        title = Title("5. Reconstruction: Order from Chaos")
        self.add(title)

        # Split screen setup using graphs to represent the complex wave
        ax_left = Axes(x_range=[0, 10], y_range=[-3, 3], x_length=5, y_length=4).shift(LEFT * 3 + DOWN * 0.5)
        ax_right = Axes(x_range=[0, 10], y_range=[-3, 3], x_length=5, y_length=4).shift(RIGHT * 3 + DOWN * 0.5)

        title_left = Text("Original Data", font_size=24).next_to(ax_left, UP)
        title_right = Text("POD Reconstruction", font_size=24).next_to(ax_right, UP)

        self.play(Create(ax_left), Create(ax_right), Write(title_left), Write(title_right))

        # Generate a noisy/complex signal
        def chaotic_wave(x):
            return np.sin(x) + 0.5 * np.sin(3 * x + 1) + 0.3 * np.sin(7 * x) + np.random.normal(0, 0.1)

        # The pure modes
        mode1_func = lambda x: np.sin(x)
        mode2_func = lambda x: np.sin(x) + 0.5 * np.sin(3 * x + 1)
        mode3_func = lambda x: np.sin(x) + 0.5 * np.sin(3 * x + 1) + 0.3 * np.sin(7 * x)

        original_graph = ax_left.plot(chaotic_wave, color=RED, use_smoothing=False)
        self.play(Create(original_graph), run_time=2)

        # Reconstruct step by step on the right
        recon_text = Text("Mode 1", font_size=20, color=YELLOW).next_to(ax_right, DOWN)
        recon_graph = ax_right.plot(mode1_func, color=YELLOW)
        
        self.play(Create(recon_graph), Write(recon_text))
        self.wait(1)

        recon_text_2 = Text("Mode 1 + Mode 2", font_size=20, color=GREEN).next_to(ax_right, DOWN)
        recon_graph_2 = ax_right.plot(mode2_func, color=GREEN)
        
        self.play(
            Transform(recon_graph, recon_graph_2),
            Transform(recon_text, recon_text_2)
        )
        self.wait(1)

        recon_text_3 = Text("Mode 1 + 2 + 3", font_size=20, color=BLUE).next_to(ax_right, DOWN)
        recon_graph_3 = ax_right.plot(mode3_func, color=BLUE)

        self.play(
            Transform(recon_graph, recon_graph_3),
            Transform(recon_text, recon_text_3)
        )
        self.wait(3)

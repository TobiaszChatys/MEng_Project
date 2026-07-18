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

class PODExplainer(Scene):
    def construct(self):
        self.camera.background_color = RP_BASE

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
    # SCENE 1: Dynamic Fluid Flow
    # ---------------------------------------------------------
    def scene_1_chaos(self):
        title = Title("1. The Problem: Drowning in Data", color=RP_TEXT)
        self.add(title)

        # Use a ValueTracker to make the fluid actually move over time
        t_tracker = ValueTracker(0)

        def fluid_func(pos):
            t = t_tracker.get_value()
            x, y = pos[:2]
            # Add 't' to the trig functions so the waves travel
            u = np.sin(y + t) + np.cos(x * 0.5 - t) * np.sin(y * 2)
            v = np.cos(x + t) - np.sin(x * 0.5 + t)
            return np.array([u, v, 0])

        vector_field = ArrowVectorField(
            fluid_func, 
            x_range=[-6, 6, 0.75], 
            y_range=[-3, 3, 0.75],
            colors=[RP_PINE, RP_FOAM, RP_IRIS]
        )
        
        # Add updater so the field redraws every frame based on t_tracker
        vector_field.add_updater(lambda vf: vf.become(
            ArrowVectorField(fluid_func, x_range=[-6, 6, 0.75], y_range=[-3, 3, 0.75], colors=[RP_PINE, RP_FOAM, RP_IRIS])
        ))

        self.play(FadeIn(vector_field))
        # Animate the flow for 4 seconds
        self.play(t_tracker.animate.set_value(3), run_time=4, rate_func=linear)

        # "The fluid freezes. A grid overlays the fluid."
        vector_field.clear_updaters() # Freeze the fluid
        
        grid = NumberPlane(
            x_range=[-7, 7, 1], y_range=[-4, 4, 1],
            background_line_style={"stroke_color": RP_MUTED, "stroke_opacity": 0.4},
            faded_line_style={"stroke_color": RP_MUTED, "stroke_opacity": 0.1},
            axis_config={"stroke_color": RP_TEXT}
        )
        self.play(Create(grid), vector_field.animate.set_opacity(0.3))
        self.wait(2)

    # ---------------------------------------------------------
    # SCENE 2: Fluid Matrix Assembly
    # ---------------------------------------------------------
    def scene_2_snapshots(self):
        title = Title("2. The Method of Snapshots", color=RP_TEXT)
        self.add(title)

        # Start with the 2D grid
        grid_matrix = Matrix([
            ["u_{11}", "u_{12}"],
            ["u_{21}", "u_{22}"]
        ], element_to_mobject_config={"color": RP_TEXT}).shift(LEFT * 4)
        grid_matrix.get_brackets().set_color(RP_TEXT)
        
        # We flatten it into a column
        col_vector = Matrix([
            ["u_{11}"], ["u_{12}"], ["u_{21}"], ["u_{22}"]
        ], element_to_mobject_config={"color": RP_TEXT})
        col_vector.get_brackets().set_color(RP_TEXT)

        self.play(Write(grid_matrix))
        self.wait(1)
        
        # Smoothly morph the 2D matrix into a 1D column
        self.play(ReplacementTransform(grid_matrix, col_vector), run_time=1.5)
        self.play(col_vector.animate.shift(LEFT * 4))

        # Now, instead of just appearing, we slide the column into the big matrix
        matrix_name = MathTex("X = ").shift(LEFT * 1.5)
        
        big_matrix = Matrix([
            ["u_1(t_1)", "u_1(t_2)", "\\dots", "u_1(t_m)"],
            ["u_2(t_1)", "u_2(t_2)", "\\dots", "u_2(t_m)"],
            ["\\vdots", "\\vdots", "\\ddots", "\\vdots"],
            ["u_n(t_1)", "u_n(t_2)", "\\dots", "u_n(t_m)"]
        ], element_to_mobject_config={"color": RP_TEXT}).next_to(matrix_name, RIGHT)
        big_matrix.get_brackets().set_color(RP_TEXT)

        # Highlight the first column of the big matrix
        first_col_box = SurroundingRectangle(big_matrix.get_columns()[0], color=RP_FOAM)

        self.play(Write(matrix_name))
        self.play(
            ReplacementTransform(col_vector, big_matrix.get_columns()[0]),
            Create(big_matrix.get_brackets()),
            Create(first_col_box)
        )
        # Fade in the rest of the snapshots to represent time
        self.play(
            FadeIn(big_matrix.get_columns()[1:]),
            first_col_box.animate.stretch_to_fit_width(big_matrix.width).move_to(big_matrix),
            run_time=2
        )
        self.play(FadeOut(first_col_box))
        self.wait(1)

    # ---------------------------------------------------------
    # SCENE 3: Core Idea (Unchanged for brevity, it's already solid)
    # ---------------------------------------------------------
    def scene_3_energy_orthogonality(self):
        title = Title("3. Finding the Best Perspective", color=RP_TEXT)
        self.add(title)
        axes = Axes(x_range=[-4, 4], y_range=[-3, 3], x_length=8, y_length=6, axis_config={"color": RP_TEXT})
        
        np.random.seed(42)
        points = []
        for _ in range(150):
            x = np.random.normal(0, 2)
            y = np.random.normal(0, 0.4)
            points.append([x * np.cos(PI/6) - y * np.sin(PI/6), x * np.sin(PI/6) + y * np.cos(PI/6), 0])
            
        dots = VGroup(*[Dot(axes.c2p(p[0], p[1]), color=RP_FOAM, radius=0.06) for p in points])

        self.play(Create(axes), FadeIn(dots, lag_ratio=0.02), run_time=2)
        
        mode1_vector = Arrow(axes.c2p(0, 0), axes.c2p(3 * np.cos(PI/6), 3 * np.sin(PI/6)), buff=0, color=RP_LOVE, stroke_width=6)
        mode2_vector = Arrow(axes.c2p(0, 0), axes.c2p(-1 * np.sin(PI/6), 1 * np.cos(PI/6)), buff=0, color=RP_GOLD, stroke_width=6)

        self.play(GrowArrow(mode1_vector), Write(Text("Mode 1 (Max Energy)", color=RP_LOVE, font_size=24).next_to(mode1_vector, RIGHT)))
        self.play(GrowArrow(mode2_vector), Write(Text("Mode 2", color=RP_GOLD, font_size=24).next_to(mode2_vector, UP)))
        self.wait(1)

    # ---------------------------------------------------------
    # SCENE 4: The Math (Unchanged)
    # ---------------------------------------------------------
    def scene_4_math(self):
        title = Title("4. The Math (Eigen-decomposition)", color=RP_TEXT)
        self.add(title)

        eq1 = MathTex("C", "=", "X^T", "X", font_size=60)
        eq1_desc = Text("Covariance Matrix (Spatial Correlation)", font_size=24, color=RP_MUTED).next_to(eq1, DOWN)
        self.play(Write(eq1), FadeIn(eq1_desc))
        self.play(eq1.animate.shift(UP * 2), eq1_desc.animate.shift(UP * 2))

        eq2 = MathTex("C", "\\mathbf{v}", "=", "\\lambda", "\\mathbf{v}", font_size=60)
        self.play(Write(eq2))

        box_v = SurroundingRectangle(eq2[1], color=RP_PINE, corner_radius=0.1)
        box_lambda = SurroundingRectangle(eq2[3], color=RP_GOLD, corner_radius=0.1)
        self.play(Create(box_v), Write(Text("Spatial Modes", color=RP_PINE, font_size=24).next_to(box_v, DOWN)))
        self.play(Create(box_lambda), Write(Text("Energy", color=RP_GOLD, font_size=24).next_to(box_lambda, UP)))
        self.wait(1)

    # ---------------------------------------------------------
    # SCENE 5: Up to 10 Modes Merging
    # ---------------------------------------------------------
    def scene_5_reconstruction(self):
        title = Title("5. Reconstruction: Order from Chaos", color=RP_TEXT)
        self.add(title)

        ax_config = {"color": RP_TEXT, "stroke_width": 2}
        ax_left = Axes(x_range=[0, 10], y_range=[-4, 4], x_length=5, y_length=4, axis_config=ax_config).shift(LEFT * 3 + DOWN * 0.5)
        ax_right = Axes(x_range=[0, 10], y_range=[-4, 4], x_length=5, y_length=4, axis_config=ax_config).shift(RIGHT * 3 + DOWN * 0.5)

        self.play(
            Create(ax_left), Create(ax_right), 
            Write(Text("Original Data", font_size=24).next_to(ax_left, UP)), 
            Write(Text("POD Reconstruction", font_size=24).next_to(ax_right, UP))
        )

        # Generate 10 specific wave frequencies and amplitudes to simulate POD modes
        np.random.seed(42)
        amplitudes = [1.5, 0.8, 0.6, 0.4, 0.3, 0.2, 0.15, 0.1, 0.08, 0.05]
        frequencies = [1.0, 2.5, 4.1, 5.8, 7.4, 9.1, 11.0, 13.5, 16.2, 19.0]
        phases = np.random.uniform(0, 2 * PI, 10)

        # Function factories to prevent variable scoping issues in the loop
        def get_single_mode_func(idx):
            return lambda x: amplitudes[idx] * np.sin(frequencies[idx] * x + phases[idx])

        def get_sum_modes_func(num_modes):
            return lambda x: sum(amplitudes[j] * np.sin(frequencies[j] * x + phases[j]) for j in range(num_modes))

        # Plot the highly chaotic original signal on the left
        original_graph = ax_left.plot(get_sum_modes_func(10), color=RP_ROSE, use_smoothing=True)
        self.play(Create(original_graph), run_time=2)

        # Start the reconstruction on the right with Mode 1
        current_recon_graph = ax_right.plot(get_sum_modes_func(1), color=RP_PINE)
        mode_text = Text("Mode 1", font_size=20, color=RP_PINE).next_to(ax_right, DOWN)
        
        self.play(Create(current_recon_graph), Write(mode_text))
        self.wait(0.5)

        # Loop to add modes 2 through 10
        for i in range(1, 10):
            # 1. Show the isolated new mode fading in lightly
            new_mode_graph = ax_right.plot(get_single_mode_func(i), color=RP_GOLD).set_opacity(0.6)
            self.play(FadeIn(new_mode_graph), run_time=0.4)

            # 2. Calculate what the graph will look like when added together
            next_recon_graph = ax_right.plot(get_sum_modes_func(i + 1), color=RP_PINE)
            next_mode_text = Text(f"Modes 1 to {i + 1}", font_size=20, color=RP_PINE).next_to(ax_right, DOWN)

            # 3. Visually morph the existing graph AND the new small mode into the combined graph
            # Speeds up as it goes to higher modes to simulate "filling in the details"
            speed = max(0.2, 1.0 - (i * 0.1)) 
            
            self.play(
                ReplacementTransform(VGroup(current_recon_graph, new_mode_graph), next_recon_graph),
                Transform(mode_text, next_mode_text),
                run_time=speed
            )
            
            # Update the reference for the next loop
            current_recon_graph = next_recon_graph

        self.wait(2)

        # Final flourish: Show that the right side now matches the left side perfectly
        match_text = Text("Noise Filtered. Dynamics Captured.", font_size=24, color=RP_LOVE).next_to(mode_text, DOWN)
        self.play(Write(match_text), current_recon_graph.animate.set_color(RP_ROSE))
        self.wait(3)

from manim import *

class PODIntro(Scene):
    def construct(self):
        # 1. Create a title
        title = Tex("Proper Orthogonal Decomposition").to_edge(UP)
        self.play(Write(title))
        
        # 2. Represent PIV Snapshots as column vectors
        snapshot_1 = Matrix([["u_1(t_1)"], ["u_2(t_1)"], ["\\vdots"], ["u_n(t_1)"]])
        snapshot_2 = Matrix([["u_1(t_2)"], ["u_2(t_2)"], ["\\vdots"], ["u_n(t_2)"]])
        
        # Position them on screen
        snapshot_1.shift(LEFT * 3)
        snapshot_2.next_to(snapshot_1, RIGHT, buff=1)
        
        # Animate the vectors appearing
        self.play(FadeIn(snapshot_1), FadeIn(snapshot_2))
        self.wait(1)
        
        # 3. Combine them into the Snapshot Matrix 'X'
        matrix_X = Matrix([
            ["u_1(t_1)", "u_1(t_2)", "\\dots"],
            ["u_2(t_1)", "u_2(t_2)", "\\dots"],
            ["\\vdots", "\\vdots", "\\ddots"],
            ["u_n(t_1)", "u_n(t_2)", "\\dots"]
        ])
        
        matrix_name = MathTex("X = ").next_to(matrix_X, LEFT)
        X_group = VGroup(matrix_name, matrix_X)
        
        # Transform the separate snapshots into the big matrix
        self.play(
            ReplacementTransform(snapshot_1, matrix_X),
            ReplacementTransform(snapshot_2, matrix_X),
            Write(matrix_name)
        )
        self.wait(2)
        
        # 4. Show the Eigenvalue equation
        eq = MathTex("C", "\\mathbf{v}", "=", "\\lambda", "\\mathbf{v}")
        eq.next_to(X_group, DOWN, buff=1)
        
        self.play(Write(eq))
        
        # Highlight Spatial Modes vs Energy
        framebox_v = SurroundingRectangle(eq[1], color=BLUE)
        label_v = Text("Spatial Mode", color=BLUE).scale(0.5).next_to(framebox_v, DOWN)
        
        framebox_lambda = SurroundingRectangle(eq[3], color=YELLOW)
        label_lambda = Text("Energy", color=YELLOW).scale(0.5).next_to(framebox_lambda, UP)
        
        self.play(Create(framebox_v), Write(label_v))
        self.play(Create(framebox_lambda), Write(label_lambda))
        self.wait(3)

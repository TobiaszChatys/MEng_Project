"""
Slide 4 — Static figures for POD explanation
=============================================
Generates three publication-ready PNGs for the presentation:

    1. pca_analogy.png       — scatter plot with rotated principal axes
    2. snapshot_matrix.png   — PIV grid → column vector → matrix X diagram
    3. eigenvalue_spectrum.png — real eigenvalue bars + cumulative energy

Run
---
    cd scripts
    source .venv/bin/activate
    python slide4_figures.py

Output lands in  scripts/slide4_output/
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import numpy as np

# ── Paths ────────────────────────────────────────────────────────────
OUT_DIR = Path(__file__).parent / "slide4_output"
OUT_DIR.mkdir(exist_ok=True)

DATA_PATH = Path(__file__).with_name("pod_data_cache.json")
with open(DATA_PATH) as f:
    POD = json.load(f)

EIGENVALUES = np.asarray(POD["eigenvalues"][:60])
CUM_ENERGY  = np.asarray(POD["cumulative_energy"][:60])
N_MODES_99  = POD["n_modes_99pct"]
TOTAL_MODES = POD["total_modes"]
GRID_ROWS   = POD["grid_rows"]
GRID_COLS   = POD["grid_cols"]

# ── MATLAB winter colourmap colours ──────────────────────────────────
SIGNAL  = "#0055d4"     # blue  — modes / structure
TEAL    = "#008080"     # mid-winter
SEAFOAM = "#00b366"     # end-winter
NOISE   = "#cc2222"     # red   — noise
ENERGY  = "#d48a00"     # amber — energy curve
MUTED   = "#888888"
TEXT    = "#222222"

# ── Global plot style ────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":       "serif",
    "mathtext.fontset":  "cm",
    "axes.labelsize":    18,
    "axes.titlesize":    20,
    "xtick.labelsize":   14,
    "ytick.labelsize":   14,
    "axes.linewidth":    1.2,
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
    "savefig.facecolor": "white",
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
    "savefig.pad_inches": 0.15,
})


# ══════════════════════════════════════════════════════════════════════
# Figure 1 — PCA scatter-plot analogy
# ══════════════════════════════════════════════════════════════════════
def fig_pca_analogy():
    fig, ax = plt.subplots(figsize=(7, 5.5))

    np.random.seed(42)
    n = 200
    angle = np.radians(30)
    R = np.array([[np.cos(angle), -np.sin(angle)],
                  [np.sin(angle),  np.cos(angle)]])
    raw = np.column_stack([
        np.random.normal(0, 2.0, n),
        np.random.normal(0, 0.5, n),
    ])
    pts = (R @ raw.T).T

    ax.scatter(pts[:, 0], pts[:, 1], s=12, alpha=0.5, color=TEAL,
               edgecolors="none", zorder=2)

    # Original axes (grey, thin)
    ax.axhline(0, color=MUTED, lw=0.8, ls="--", zorder=1)
    ax.axvline(0, color=MUTED, lw=0.8, ls="--", zorder=1)

    # Principal component arrows
    mode1_dir = R @ np.array([3.5, 0])
    mode2_dir = R @ np.array([0, 1.5])

    ax.annotate("", xy=mode1_dir, xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color=SIGNAL,
                                lw=2.5, mutation_scale=18),
                zorder=5)
    ax.annotate("", xy=mode2_dir, xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color=NOISE,
                                lw=2.5, mutation_scale=18),
                zorder=5)

    ax.text(mode1_dir[0] + 0.15, mode1_dir[1] + 0.25,
            r"$\phi_1$ — max energy", fontsize=15, color=SIGNAL,
            fontweight="bold", ha="left", va="bottom")
    ax.text(mode2_dir[0] - 0.15, mode2_dir[1] + 0.15,
            r"$\phi_2$", fontsize=15, color=NOISE,
            fontweight="bold", ha="right", va="bottom")

    ax.set_xlim(-5, 5)
    ax.set_ylim(-3.5, 3.5)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$V_x$")
    ax.set_ylabel(r"$V_y$")
    ax.set_title("POD as PCA: rotating axes to align with maximum energy",
                 fontsize=16, pad=12)

    # Equations in lower-right
    eq_box = (
        r"$u(\mathbf{x}, t) = \overline{u}(\mathbf{x})"
        r" + \sum_{k} a_k(t)\,\phi_k(\mathbf{x})$"
    )
    ax.text(0.97, 0.05, eq_box, transform=ax.transAxes,
            fontsize=14, ha="right", va="bottom", color=TEXT,
            bbox=dict(boxstyle="round,pad=0.4", fc="white",
                      ec=MUTED, alpha=0.9))

    ax.tick_params(direction="in", top=False, right=False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    out = OUT_DIR / "pca_analogy.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  ✓ {out}")


# ══════════════════════════════════════════════════════════════════════
# Figure 2 — Snapshot matrix diagram
# ══════════════════════════════════════════════════════════════════════
def fig_snapshot_matrix():
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.set_xlim(-0.5, 10.5)
    ax.set_ylim(-1.8, 3.8)
    ax.set_aspect("equal")
    ax.axis("off")

    # --- left: small PIV grid ---
    gx, gy = 4, 4
    piv_ox, piv_oy = 0.3, 0.5
    sp = 0.6
    np.random.seed(7)

    for i in range(gy):
        for j in range(gx):
            x0 = piv_ox + j * sp
            y0 = piv_oy + (gy - 1 - i) * sp
            dx = 0.25 + 0.1 * np.random.randn()
            dy = 0.08 * np.random.randn()
            ax.annotate("", xy=(x0 + dx, y0 + dy), xytext=(x0, y0),
                        arrowprops=dict(arrowstyle="-|>", color=SIGNAL,
                                        lw=1.2, mutation_scale=10))

    grid_rect = mpatches.FancyBboxPatch(
        (piv_ox - 0.2, piv_oy - 0.25),
        gx * sp + 0.1, gy * sp - 0.05,
        boxstyle="round,pad=0.08", linewidth=1.2,
        edgecolor=MUTED, facecolor="none")
    ax.add_patch(grid_rect)
    ax.text(piv_ox + gx * sp / 2 - 0.1, piv_oy - 0.6,
            "PIV snapshot", fontsize=11, ha="center", color=MUTED)

    # --- arrow →  column vector ---
    ax.annotate("", xy=(3.6, 1.5), xytext=(3.0, 1.5),
                arrowprops=dict(arrowstyle="-|>", color=TEXT, lw=1.5))

    # --- column vector ---
    col_x = 4.1
    col_entries = [r"$V_{x_1}$", r"$V_{y_1}$", r"$\vdots$",
                   r"$V_{x_n}$", r"$V_{y_n}$"]
    col_h = len(col_entries) * 0.45
    col_bot = 1.5 - col_h / 2

    rect_col = mpatches.FancyBboxPatch(
        (col_x - 0.1, col_bot - 0.05), 0.75, col_h + 0.1,
        boxstyle="round,pad=0.05", linewidth=1,
        edgecolor=TEXT, facecolor="#f0f7ff")
    ax.add_patch(rect_col)

    for idx, txt in enumerate(col_entries):
        ax.text(col_x + 0.28, col_bot + col_h - (idx + 0.5) * 0.45,
                txt, fontsize=11, ha="center", va="center", color=TEXT)

    ax.text(col_x + 0.28, col_bot - 0.35,
            r"$\mathbf{x}_k$", fontsize=13, ha="center", color=SIGNAL,
            fontweight="bold")

    # --- arrow →  matrix X ---
    ax.annotate("", xy=(5.6, 1.5), xytext=(5.15, 1.5),
                arrowprops=dict(arrowstyle="-|>", color=TEXT, lw=1.5))

    # --- matrix X ---
    mat_x = 6.0
    mat_w = 4.0
    mat_h = 2.6
    mat_bot = 1.5 - mat_h / 2

    rect_mat = mpatches.FancyBboxPatch(
        (mat_x - 0.1, mat_bot), mat_w, mat_h,
        boxstyle="round,pad=0.08", linewidth=1.2,
        edgecolor=TEXT, facecolor="#f0f7ff")
    ax.add_patch(rect_mat)

    col_positions = [0.4, 1.2, 2.4, 3.4]
    col_labels = [r"$\mathbf{x}_1$", r"$\mathbf{x}_2$",
                  r"$\cdots$", r"$\mathbf{x}_{2500}$"]

    for cx, lbl in zip(col_positions, col_labels):
        xx = mat_x + cx
        if lbl != r"$\cdots$":
            ax.plot([xx, xx], [mat_bot + 0.3, mat_bot + mat_h - 0.3],
                    color=SIGNAL, lw=1, alpha=0.4)
        ax.text(xx, 1.5, lbl, fontsize=12, ha="center", va="center",
                color=TEXT)

    ax.text(mat_x - 0.35, 1.5, r"$X =$", fontsize=16,
            ha="right", va="center", color=ENERGY, fontweight="bold")

    # braces
    brace_y = mat_bot - 0.15
    ax.annotate("", xy=(mat_x + mat_w - 0.1, brace_y - 0.1),
                xytext=(mat_x - 0.1, brace_y - 0.1),
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0))
    ax.text(mat_x + mat_w / 2 - 0.05, brace_y - 0.35,
            "2,500 snapshots", fontsize=11, ha="center", color=MUTED)

    dim_text = f"{2 * GRID_ROWS * GRID_COLS:,} rows"
    ax.text(mat_x + mat_w + 0.25, 1.5, dim_text,
            fontsize=10, ha="left", va="center", color=MUTED,
            rotation=90)

    out = OUT_DIR / "snapshot_matrix.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  ✓ {out}")


# ══════════════════════════════════════════════════════════════════════
# Figure 3 — Eigenvalue spectrum + cumulative energy (real data)
# ══════════════════════════════════════════════════════════════════════
def fig_eigenvalue_spectrum():
    N = 50
    eig = EIGENVALUES[:N]
    eig_norm = eig / eig.max()
    cum = CUM_ENERGY[:N]

    fig, ax1 = plt.subplots(figsize=(9, 5))

    # Winter colourmap gradient for bars
    colours = []
    for i in range(N):
        t = i / (N - 1)
        r = 0
        g = t
        b = 1.0 - 0.5 * t
        colours.append((r, g, b))

    ax1.bar(np.arange(1, N + 1), eig_norm, color=colours,
            edgecolor="white", linewidth=0.3, width=0.8, zorder=3)

    ax1.set_xlabel("Mode number")
    ax1.set_ylabel(r"Normalised eigenvalue $(\lambda_k / \lambda_1)$")
    ax1.set_xlim(0, N + 1)
    ax1.set_ylim(0, 1.12)
    ax1.tick_params(direction="in")

    for spine in ("top", "right"):
        ax1.spines[spine].set_visible(False)

    # Cumulative energy on secondary axis
    ax2 = ax1.twinx()
    ax2.plot(np.arange(1, N + 1), cum * 100, color=ENERGY,
             lw=2.5, marker="o", markersize=3, zorder=4)
    ax2.set_ylabel(r"Cumulative energy (%)", color=ENERGY)
    ax2.tick_params(axis="y", colors=ENERGY, direction="in")
    ax2.set_ylim(0, 65)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_color(ENERGY)

    # Percentage annotations
    for pct in [10, 20, 30, 40, 50]:
        frac = pct / 100
        idx = int(np.searchsorted(cum, frac))
        if idx < N:
            ax2.annotate(
                f"{pct}%",
                xy=(idx + 1, cum[idx] * 100),
                xytext=(idx + 4, cum[idx] * 100 + 3),
                fontsize=10, color=ENERGY,
                arrowprops=dict(arrowstyle="-|>", color=ENERGY,
                                lw=0.8, mutation_scale=8),
                ha="left",
            )

    # Equation
    eq = r"$\lambda_1 > \lambda_2 > \cdots > \lambda_n$"
    ax1.text(0.97, 0.92, eq, transform=ax1.transAxes,
             fontsize=16, ha="right", va="top", color=TEXT,
             bbox=dict(boxstyle="round,pad=0.4", fc="white",
                       ec=MUTED, alpha=0.9))

    # Footnote
    footnote = (
        f"Liquid phase only  |  "
        f"{TOTAL_MODES:,} total modes  |  "
        f"{N_MODES_99} modes for 99% energy"
    )
    ax1.text(0.5, -0.14, footnote, transform=ax1.transAxes,
             fontsize=10, ha="center", color=MUTED)

    out = OUT_DIR / "eigenvalue_spectrum.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  ✓ {out}")


# ══════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("Generating Slide 4 figures...")
    fig_pca_analogy()
    fig_snapshot_matrix()
    fig_eigenvalue_spectrum()
    print(f"\nAll figures saved to {OUT_DIR}/")

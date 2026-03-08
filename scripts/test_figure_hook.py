"""Smoke test for the Claude Code figure quality hook.

Run this script once to verify the PostToolUse hook fires correctly:

    python scripts/test_figure_hook.py

Expected behaviour:
  1. Script saves a simple 2-panel figure to
     projects/visium/ggo_visium/figures/exploratory/hook_smoke_test.png
  2. Claude Code detects the save (keyword: savefig) and fires the hook.
  3. Hook finds the new PNG, builds the evaluation prompt, prints it to stdout.
  4. Claude reads the figure and writes the analytics file:
     projects/visium/ggo_visium/figures/exploratory/hook_smoke_test_claude_analytics_and_legend.txt

If the analytics file appears within a few seconds of the script completing,
the hook is working end-to-end.

To clean up afterwards:
    rm projects/visium/ggo_visium/figures/exploratory/hook_smoke_test*
"""

import os
import sys

try:
    import matplotlib
    matplotlib.use("Agg")  # non-interactive backend
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError:
    print("matplotlib / numpy not installed in this environment.", file=sys.stderr)
    print("Activate an environment that has them, e.g.:", file=sys.stderr)
    print("  conda activate sc_tools && python scripts/test_figure_hook.py", file=sys.stderr)
    sys.exit(1)

# Output path — must be inside a recognised figures/ subdir (not scratch/)
OUT_DIR = "projects/visium/ggo_visium/figures/exploratory"
OUT_PATH = os.path.join(OUT_DIR, "hook_smoke_test.png")

os.makedirs(OUT_DIR, exist_ok=True)

rng = np.random.default_rng(42)

fig, axes = plt.subplots(1, 2, figsize=(8, 4))

# Panel A — scatter plot
x = rng.normal(0, 1, 80)
y = 0.7 * x + rng.normal(0, 0.5, 80)
axes[0].scatter(x, y, color="#E69F00", alpha=0.7, edgecolors="none", s=40)
axes[0].set_xlabel("Gene A expression (log1p)")
axes[0].set_ylabel("Gene B expression (log1p)")
axes[0].set_title("(a) Gene correlation")

# Panel B — bar plot with error bars
groups = ["Control", "Treatment"]
means = [1.2, 2.8]
sems = [0.15, 0.22]
axes[1].bar(groups, means, yerr=sems, color=["#56B4E9", "#D55E00"],
            capsize=5, edgecolor="black", linewidth=0.8)
axes[1].set_ylabel("Normalised expression")
axes[1].set_title("(b) Group comparison")
axes[1].set_ylim(0, 3.5)

fig.suptitle("Hook smoke test — exploratory figure", fontsize=10)
fig.tight_layout()
fig.savefig(OUT_PATH, dpi=150, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {OUT_PATH}")
print("Hook should now fire and write:")
print(f"  {OUT_PATH.replace('.png', '_claude_analytics_and_legend.txt')}")

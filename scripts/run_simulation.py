import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8.5,
    "mathtext.fontset": "cm",
    "pdf.fonttype": 42,
})

class GalacticEcosystem:
    EMPTY = 0
    INTROVERT = 1
    EXTROVERT = 2
    RUIN = 3

    def __init__(self, grid_size=60, p_birth=0.0001, p_mutate=0.0005,
                 p_expand=0.05, s_travel_cost=0.3, s_purge=0.99,
                 decay_introvert=0.0001, decay_extrovert=0.05, decay_ruin=0.005):
        self.L = grid_size
        self.grid = np.zeros((self.L, self.L), dtype=np.int8)
        self.p_birth = p_birth
        self.p_mutate = p_mutate
        self.p_expand = p_expand
        self.s_travel_cost = s_travel_cost
        self.s_purge = s_purge
        self.decay_introvert = decay_introvert
        self.decay_extrovert = decay_extrovert
        self.decay_ruin = decay_ruin
        self.history = {"empty": [], "introvert": [], "extrovert": [], "ruin": []}

    def _get_neighbors(self, y, x):
        neighbors = []
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if dy == 0 and dx == 0:
                    continue
                neighbors.append(((y + dy) % self.L, (x + dx) % self.L))
        return neighbors

    def step(self):
        new_grid = self.grid.copy()
        for y in range(self.L):
            for x in range(self.L):
                state = self.grid[y, x]
                if state == self.EMPTY:
                    if np.random.rand() < self.p_birth:
                        new_grid[y, x] = self.INTROVERT
                elif state == self.INTROVERT:
                    if np.random.rand() < self.p_mutate:
                        new_grid[y, x] = self.EXTROVERT
                    elif np.random.rand() < self.decay_introvert:
                        new_grid[y, x] = self.EMPTY
                elif state == self.EXTROVERT:
                    if np.random.rand() < self.p_expand:
                        neighbors = self._get_neighbors(y, x)
                        ty, tx = neighbors[np.random.randint(len(neighbors))]
                        t_state = self.grid[ty, tx]
                        if np.random.rand() > self.s_travel_cost:
                            if t_state == self.EMPTY:
                                new_grid[ty, tx] = self.EXTROVERT
                            elif t_state in (self.INTROVERT, self.EXTROVERT, self.RUIN):
                                if np.random.rand() < self.s_purge:
                                    new_grid[y, x] = self.RUIN
                                    if t_state == self.EXTROVERT:
                                        new_grid[ty, tx] = self.RUIN
                    if np.random.rand() < self.decay_extrovert:
                        new_grid[y, x] = self.RUIN
                elif state == self.RUIN:
                    if np.random.rand() < self.decay_ruin:
                        new_grid[y, x] = self.EMPTY

        self.grid = new_grid
        total = self.L * self.L
        self.history["empty"].append(np.sum(self.grid == self.EMPTY) / total)
        self.history["introvert"].append(np.sum(self.grid == self.INTROVERT) / total)
        self.history["extrovert"].append(np.sum(self.grid == self.EXTROVERT) / total)
        self.history["ruin"].append(np.sum(self.grid == self.RUIN) / total)

    def run(self, steps=10000):
        for _ in range(steps):
            self.step()

if __name__ == "__main__":
    print("Running 10,000-step galactic evolution...")
    sim = GalacticEcosystem()
    sim.run(steps=10000)

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.8))
    cmap = ListedColormap(["#0d1117", "#58a6ff", "#f85149", "#8b949e"])
    axes[0].imshow(sim.grid, cmap=cmap, vmin=0, vmax=3)
    axes[0].set_title("Galactic Snapshot (Step 10000)")
    axes[0].axis("off")

    steps_range = range(len(sim.history["introvert"]))
    axes[1].plot(steps_range, sim.history["introvert"], label="Introvert (A - Silent)", color="#1f77b4", lw=1.8)
    axes[1].plot(steps_range, sim.history["extrovert"], label="Extrovert (a - Visible)", color="#d62728", lw=1.5)
    axes[1].plot(steps_range, sim.history["ruin"], label="Ruins / Scars", color="#7f7f7f", linestyle="--")
    axes[1].set_xlabel("Time Steps")
    axes[1].set_ylabel("Occupancy Ratio")
    axes[1].set_title("Civilization Dynamics Over Time")
    axes[1].legend(loc="center right")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("figures/fig1_dynamics.pdf", bbox_inches="tight")
    print("Saved: figures/fig1_dynamics.pdf")

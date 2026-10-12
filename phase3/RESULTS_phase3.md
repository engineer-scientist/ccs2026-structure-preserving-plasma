
## Phase 3 (Oct 12, 2026): locate imbalance diagnosed and fixed, then weak scaling to 16 nodes
Approved 12:40 AM ET Oct 12. Running log: `phase3/NOTES.md`. Scripts: `phase3/scripts/` (copies in `~/grokbot-sp/p3/` on both machines).
Figure: `figs/weak_scaling_v2.png`. Data: `weak_scaling_v2.csv` (combined), `aurora_p3diag_8922161_reduced.csv`, `aurora_p3ws_8922217_reduced.csv`, `polaris_p3ws_7738293_reduced.csv`.

### Step A: why one rank was 18–41× slower in DMPlexLocatePoints
* **Cause (confirmed three ways).** `DMPlexCreateGridHash()` (PETSc 3.26 `plexgeometry.c`) sizes each rank's point-location hash box from *all* of that rank's local coordinates.
  On a periodic 1D mesh, the rank that owns the wrap-around cell (the last rank) also holds the vertex at x = 0, so its hash box spans the
  **whole domain [0, L]** while it gets only ~0.8 × 256 ≈ 204 bins. Its 256 cells then fall into one or two bins, and every particle that left its
  previous cell costs a linear scan over ~256 cells instead of ~1.
  1. `-dm_plex_print_locate` on the box (4 ranks): ranks 0–2 report hash boxes equal to their own subdomain; rank 3 reports `(0 -- 12.566)`.
  2. Per-rank `-log_view :…:ascii_info_detail` on Aurora (job 8922161, 96 ranks): DMPlexLocatePoints time **131 s on rank 95 vs a median of 7.4 s** (18.6×);
     flops on rank 95 are 19× those on any other rank. Every other rank waits for it in the next collective, so the whole step runs at rank 95's speed.
  3. The imbalance follows the formula: with 256·P/204 cells per bin on the wrap rank, it grows with P until a single bin holds all 256 cells,
     which matches the saturation at **≈ 41×** seen at 384–1,536 ranks.
  It is not the initial particle distribution (particles are created per local cell and stay balanced), not periodic coordinate localization itself,
  and not migration volume (message sizes are small and balanced). Brute-force location (`-dm_plex_hash_location 0`) makes it worse: 5.3 s per step at 96 ranks.
* **Fix.** A small local patch to ex2 (`phase3/patch/ex2_gb_wrap_hash.patch`, 27 lines, not pushed): new option `-gb_wrap_hash_factor f`. After the mesh
  is created, any rank whose local bounding box is more than 4× wider than its owned cells (only the wrap rank) sets `-dm_plex_hash_box_faces = f·width/h`
  in its own options database, so its bins are about one cell wide. All other ranks keep the default hash. We used f = 2.
  The upstream fix would be for `DMPlexCreateGridHash()` to size the box from the localized coordinates of owned cells only. That is worth reporting to PETSc.
* **Rejected first fix.** The global option `-dm_plex_hash_box_faces 2·256·P` works at 96 ranks (8.2× faster), but the hash build costs grow faster
  than linearly with the number of bins (box test: 4k bins 0.3 s, 65k 8 s, 262k 145 s). At 1,536 ranks (786k bins on *every* rank) setup never finished,
  and job 8922177 was cancelled after 10 minutes (2.7 node-hours lost). The patch gives fine bins only to the one rank that needs them.
* Effect at 96 ranks (one Aurora node): per-step time **2.29 s → 0.28–0.30 s (7.7–8.2× faster)**. Locate work imbalance (max/min flops) **19× → 1.07–1.17×**.

### Step B: grow-the-box weak scaling with the fix
Setup: 1D-1V Landau damping, k = 0.5, 256 x-cells × 200 v-cells = **51,200 particles per rank**, Δt = 0.05, 100 steps ("before": 20 steps).
**Domain L = π·P** (4π·P/4, a width of π per rank). This differs from the proposal's L = 4π·P/48: that rule gives a subdomain width of 0.26,
below |v|max·Δt = 0.3, so particles would skip a rank and be lost. With width π, every P that is a multiple of 4 holds a whole number of k = 0.5 wavelengths.
Efficiency below is T(1 rank)/T(P) from the TSStep time per step (slowest rank). The 1-node-relative value is in brackets.

| Aurora (96 ranks/node) | 1 rank | 96 (1 node) | 192 (2) | 384 (4) | 768 (8) | 1,536 (16) |
|---|---|---|---|---|---|---|
| s/step after fix | 0.147 | 0.299 | 0.330 | 0.342 | 0.363 | 0.400 |
| efficiency after fix [vs 1 node] | 1.00 | 0.49 [1.00] | 0.45 [0.91] | 0.43 [0.87] | 0.41 [0.82] | 0.37 [0.75] |
| s/step before fix | 0.151 | 2.29 | 4.32 | 4.93 | 4.85 | 4.88 |
| efficiency before fix | 1.00 | 0.066 | 0.035 | 0.031 | 0.031 | 0.031 |
| Locate max/min across ranks, before → after | – | 18.2 → 1.2 | 35.6 → 1.2 | 41.5 → 1.2 | 41.5 → 1.4 | 41.4 → 1.3 |
| particles (all conserved exactly) | 51,200 | 4.9 M | 9.8 M | 19.7 M | 39.3 M | **78.6 M** |

| Polaris (job 7738293, debug queue) | 1 rank | 32 (1 node) | 64 (2 nodes) |
|---|---|---|---|
| s/step after fix, 32 CPU ranks/node | 0.190 | 0.269 | 0.278 |
| efficiency after fix [vs 1 node] | 1.00 | 0.71 [1.00] | 0.68 [0.97] |
| s/step before fix | – | 0.932 | 1.646 |
| Locate max/min before → after | – | 6.8 → 1.2 | 12.8 → 1.1 |
| 4 ranks/node, CPU: s/step (4 / 8 ranks) | | 0.278 | 0.221 |
| 4 ranks/node, CUDA vec/mat: s/step (4 / 8 ranks) | | 0.354 | 0.327 |

* Polaris 4- and 8-node data come from job 7738212 (8 nodes, debug-scaling). It waited in the queue for ~1 h and started at ~1:57 AM ET; it runs unattended
  and writes to `~/grokbot-sp/p3/runs_polaris_p3ws_7738212`. **Those numbers are not in this table or the figure yet.**
* GPU (CUDA vectors/matrices, 1 rank per A100) is still **1.3–1.5× slower** than the same number of CPU ranks, because only the field solve is offloaded (as in Phase 2).

* **Speed-up from the fix at 16 Aurora nodes: 12.2× per step.** Across nodes the fixed code holds 75% efficiency from 1 to 16 nodes on Aurora and 97% from 1 to 2 nodes on Polaris.
  Before the fix, efficiency was 3% and was set entirely by the one wrap rank.
* **What is left.** (i) A drop within one node (1 rank → full node: 0.49 on Aurora). The point-location time per rank roughly doubles (14 s → 25 s per 100 steps)
  even though the work is balanced, which points to shared memory bandwidth and cache, not imbalance. (ii) A slow rise from 1 to 16 nodes, mostly the Poisson solve
  (SNESSolve 0.024 → 0.073 s/step, GAMG on a 393k-cell 1D mesh) and mass-matrix assembly. Location is now flat with P.
* **Invariants (the VERIFY{} contract) at every scale.** Global particle count and total weight (charge) are identical at the first and last step in every run,
  1 to 1,536 ranks, with no particle loss. The relative total-energy drift is the same at every P: **1.14e-5** after 100 steps (Δt = 0.05), and 1.12e-5 for the 20-step runs.
  So the parallel decomposition does not change the physics. (The 1-rank run, L = π, has a different box and gives 1.1e-7.)

### ALCF usage, Phase 3 (walltime × nodes from `qstat -xf`)
* Aurora: 8922147 (1 node, 0:03, failed: bad symlink) ≈ 0.00; 8922161 diagnostic (1 node, 8:27) = 0.14; **8922177 (16 nodes, 10:12, stalled in hash setup, cancelled) = 2.72 (wasted)**;
  8922217 sweep (16 nodes, 19:33) = 5.21 → **Aurora 8.08 node-hours**
* Polaris: 7738204 (cancelled before it started) = 0; 7738288 (2 nodes, 0:02, failed: wrong MPI library at link time) ≈ 0.00; 7738293 (2 nodes, 6:54) = 0.23 → **Polaris 0.23**
* **Phase 3 total so far: 8.31 node-hours** of the 25 cap. Pending: 7738212 (8 nodes, ≤ 20 min) adds ≤ 2.67, so the worst case is ≈ 11.0.

### Caveats
* The fix is a local ex2 patch, not a PETSc change. Other DMSwarm codes on periodic Plex meshes in 2D and 3D likely have the same wrap-rank hash problem,
  and the patch as written only handles 1D (`cdim == 1`).
* "Before" points in panel (a) come from two setups: Phase 2 (fixed box L = 4π, Δt = 0.01, mesh refined with P) and Phase 3 (grow box, Δt = 0.05, 20 steps).
  Both show the same pathology.
* Runs are short (100 steps, t = 5). Times exclude setup. At 1,536 ranks the wrap rank spends roughly 20 s once, at startup, building its fine hash.
* The ex2 "GPU" runs offload only the field solve. Particles stay on the CPU, as in Phase 2.

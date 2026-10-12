# Phase 3 NOTES (resume log)
Approved 12:40 AM ET Mon Oct 12 2026. Cap 25 node-h. Step A diag (Aurora 1 node <=0.5 nh), Step B weak scaling.
## Log
- 00:45 ET started; read RESULTS.md. Note: L=4pi*P/48 gives subdomain width 0.262 < |v|max*dt=0.3 -> unsafe; will pick L=4pi*P/16 (width 0.785) unless Step A says otherwise.
- 00:50 ET: ALCF SSH ControlMaster sockets for BOTH aurora and polaris are gone (~/.ssh/cm empty, dir mtime 00:35). BatchMode ssh -> Permission denied. Need Sarthak to re-login (ssh aurora / ssh polaris in WSL, keep windows open). No jobs submitted.
- Diagnosis from source (PETSc 3.26.0 plexgeometry.c DMPlexCreateGridHash): hash box = all local vertex coords + localized cell coords. Rank owning periodic wrap cell (last rank) has vertex x=0 and cell coords up to L -> its grid hash spans [0,L] with only n=floor(0.8*256)=204 bins -> its 256 cells fall into ~204/P bins -> linear scan of ~256*P/204 cells per lookup. At P=192 all 256 cells in ~1 bin => matches max/min 36x. Fix candidate (options only): -dm_plex_hash_box_faces N with N ~ 0.8*global_cells (wrap rank gets dx-size bins).
- Building MPI PETSc on box (arch-box-mpi, mpich) to verify fix locally.
- 00:45 ET: parent confirms sessions closed ~00:35; user asked to re-login. Plan: box prep, poll ssh -O check every few min, resume after.
- 00:55 ET box (mpich PETSc 3.26 arch-box-mpi) CONFIRMED: -dm_plex_print_locate shows ranks 0..P-2 hash box = own subdomain, last rank hash box = (0 -- L) whole domain. Per-rank -log_view detail at 8 ranks: rank 7 LocatePoints flops 2.34x others (time 3.65 vs ~1.5 s). Fix -dm_plex_hash_box_faces 2048 (=global cells): flop ratio 1.34; 4x global: 1.08. Chosen fix: -dm_plex_hash_box_faces 2*256*P.
- Domain choice for Step B: L = pi*P (4pi*P/4; width pi per rank, dx=pi/256) so |v|max*dt=0.3 << width, and L multiple of 4pi (k=0.5 periodic) for all P multiple of 4. Box test dt=0.05 at P=1,8: particle count exactly conserved.
- Scripts: phase3/scripts/{common_opts_v3.txt,run_common.sh,diag_aurora.pbs,ws_aurora.pbs,ws_polaris.pbs,perrank.py}
- 00:50 ET masters back. Uploaded scripts to ~/grokbot-sp/p3 on both (ex2/ex2_326 symlinks). Queue limits: Aurora debug max 2 nodes, debug-scaling max 256 nodes, max_run 1/user; Polaris debug-scaling max 10 nodes, max_queued 1/user -> submit sequentially.
- SUBMITTED Aurora 8922147 diag_aurora.pbs (debug, 1 node, 20 min). Output dir ~/grokbot-sp/p3/runs_aurora_p3diag_8922147
- SUBMITTED Polaris 7738204 ws_polaris.pbs select=1 (debug-scaling). Output runs_polaris_p3ws_1node_7738204
- Status file push command: phase3/status_cmd.txt (regenerate after editing status.txt)
- 00:58 ET: reduce3.py + poll.sh uploaded to both p3 dirs (poll: ssh host bash grokbot-sp/p3/poll.sh). 8922147 RUNNING; 7738204 queued.
- 01:03 ET: 8922147 FAILED instantly (ex2_326 symlink pointed to a directory; ~1 min of 1 node). Fixed links: Aurora p3/ex2_326 -> ../ex2_326/ex2, Polaris p3/ex2 -> ../ex2/ex2. RESUBMITTED Aurora diag as 8922161. Polaris 7738204 still queued (link fixed before start).
- 01:08 ET (user approved unattended sweep): both debug-scaling queues allow only 1 queued job/user, so each machine's sweep is ONE job using --hostfile subsets.
  * Polaris 7738204 qdel'd (never ran). SUBMITTED Polaris 7738212 ws_polaris_all.pbs (8 nodes, 20 min): 32ppn fix at 8/4/2/1 nodes + np1, GPU(4ppn cuda) & cpu4 at 8/4/2/1, unfixed 20-step baselines. Output ~/grokbot-sp/p3/runs_polaris_p3ws_7738212 (summary.txt, reduced.csv)
  * Aurora 8922172 qdel'd (never ran). SUBMITTED Aurora 8922177 ws_aurora_all.pbs (16 nodes, 20 min): 96ppn fix at 16/8/4/2/1 nodes + np1, unfixed baselines. Output runs_aurora_p3ws_8922177
  * Aurora diag 8922161 RUNNING (debug). Output runs_aurora_p3diag_8922161
  Worst-case charge: 0.33 + 16*0.33 + 8*0.33 = 8.3 nh (+~0.02 failed 8922147).
- Box has no ssh client; still using WSL.
- 01:10 ET Step A DONE: 8922161 (1 node, 04:59:21-05:07:48 UTC = 8.5 min = 0.14 nh). Results (96 ranks):
  D1 fixed-box base: TSStep 69.3 s/30 steps (2.31 s/step); LocatePoints max/min 18.6, hot = rank 95 (last, wrap), flop ratio 19.2
  D2 fixed-box fix: 8.40 s/30 (0.28 s/step) 8.2x faster; Locate flop max/min 1.15
  D3 grow base: 68.7/30 (2.29 s/step), rank 95 hot 18.3x; D4 grow fix: 27.9/100 (0.279 s/step), flop ratio 1.17
  D5 np1 fix 0.137 s/step; D6 np1 base 0.151 s/step -> eff(96) 0.06 -> 0.49
  D7 no hash (brute force), 10 steps: 5.26 s/step (worse). Np & total weight conserved in all; W drift ~1.1e-5 (dt .05).
  Residual: per-rank locate 25 s vs 14 s at np1 (node sharing / memory bandwidth, not imbalance).
- 01:15 ET 8922177 (16 nodes) STUCK in setup at 1536 ranks with global -dm_plex_hash_box_faces 786432: qdel'd after ~10.5 min (05:05:27-~05:16 UTC => ~2.8 nh WASTED). Box test: hash build cost superlinear in #bins (8 ranks: 4k bins 0.3 s, 65k 8 s, 262k 145 s).
- NEW FIX = local ex2 patch (phase3/patch/ex2_gb_wrap_hash.patch): option -gb_wrap_hash_factor f; in CreateMesh, only ranks whose local bbox width > 4x owned-cell width (the wrap rank) set -dm_plex_hash_box_faces = f*width/h in their own options DB. Box 8 ranks f=2: Locate flop max/min 1.12, TSStep 4.9 s/20; wrap rank with 786k bins costs ~20 s one-time setup only (outside TSStep).
- Built on ALCF: ~/grokbot-sp/ex2gb/ex2 (Aurora arch-aurora-cpu-opt, Polaris arch-polaris-cuda-opt). p3/ex2_326 (Aurora) and p3/ex2 (Polaris) now -> ../ex2gb/ex2. run_common.sh FIX=1 now = -gb_wrap_hash_factor 2 (old saved as run_common.sh.v1).
- Polaris 7738212 still queued -> will use patched binary + new run_common.
- RESUBMITTED Aurora sweep 8922217 (ws_aurora_all.pbs, 16 nodes, order 1->16 nodes, 'timeout 300' per run). Output runs_aurora_p3ws_8922217.
  Node-hours so far (Phase 3): 8922147 ~0.02, 8922161 0.14, 8922177 ~2.8 => ~3.0 nh. Remaining worst case: 16*0.33+8*0.33 = 8 nh => <= ~11 nh total.
- 01:41 ET Aurora 8922217: all fix runs (96..1536 + np1) done, Np conserved; base 96,192,384,768 done; base1536 likely cut by walltime. Usage: 8922147 0:03 (1n), 8922161 8:27 (1n), 8922177 10:12 (16n)=2.72nh, 8922217 ~20 min (16n) ~5.3nh
- 01:45 ET Aurora sweep 8922217 DONE (19:33 x16 = 5.21 nh). CSV saved phase3/data/aurora_p3ws_8922217_reduced.csv. TSStep s/step fix: np1 .147, 96 .299, 192 .330, 384 .342, 768 .363, 1536 .400 (eff vs 1 node: 1,.905,.874,.823,.748). Base: 96 2.29, 192 4.32, 384 4.93, 768 4.85, 1536 4.88; Locate max/min 18->41 base, 1.2-1.4 fix. Np & weight exactly conserved everywhere; W drift 1.14e-5 all scales.
  Aurora Phase-3 total: ~8.07 nh. Polaris 7738212 still queued.
- 01:55 ET figure script phase3/scripts/make_fig.py (venv python) -> figs/weak_scaling_v2.png + data/weak_scaling_v2.csv; Polaris added automatically from data/polaris_p3ws_*_reduced.csv
- RESULTS_phase3.md drafted (placeholders for Polaris + usage)
- ~01:58 ET Polaris 7738212 still Q ('Not enough free nodes', 13 queued in debug-scaling). SUBMITTED fallback Polaris 7738288 ws_polaris_debug.pbs (debug queue, 2 nodes, 20 min; 1-2 node part). Output runs_polaris_p3ws_7738288
- ~01:55 ET Polaris fallback 7738288 FAILED instantly (ex2gb built in non-login ssh shell linked libmpi_gnu_123.so.12, absent at runtime). Rebuilt ~/grokbot-sp/ex2gb/ex2 under bash -l with PrgEnv-gnu+cuda/13.0 (now libmpi_gnu.so.12, same as Phase-2 binary). Failed outputs moved to failed_runs_polaris_p3ws_7738288. RESUBMITTED fallback as 7738293 (debug, 2 nodes). 7738212 (8 nodes) still Q; will use the rebuilt binary.
- 02:00 ET Polaris fallback 7738293 DONE (6:54 x2 = 0.23 nh); CSV data/polaris_p3ws_7738293_reduced.csv. 7738212 (8 nodes) STARTED ~05:57 UTC, unattended, -> runs_polaris_p3ws_7738212.
- Mon Oct 12 02:04 EDT: RESULTS.md updated on Drive; status files pushed to aurora+polaris; total 8.31 nh; 7738212 running unattended (<=2.67 nh more). Final report sent.
- Mon Oct 12 02:19 EDT: 7738212 done (8 nodes, 17:09 = 2.29 nh, exit 0). Fig/CSV regenerated (7738293 csv moved to data/superseded). RESULTS.md + CSV updated on Drive in place; fig uploaded as weak_scaling_v3.png (1Fak8RRpGHN0IOCWSzkCUv2i-PnaOCBOL); logs 1yOO1QGVRZfPmvEWsuf0YM77YsUNzPLTP; csv 1jlIaI8E2Ub-MboMb-B9IaSATMD8eRAn3. FINAL total 10.60 nh. Phase 3 complete.

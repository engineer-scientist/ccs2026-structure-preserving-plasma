# Phase 4 NOTES (strong scaling, long-run conservation, 32/64-node weak scaling)
Approved by Sarthak Mon Oct 12 2026 (~03:30 ET): "Can all 3 be done?" Hard cap 40 nh. Access via laptop WSL ControlMaster (aurora, polaris).
Scripts: phase4/scripts/ (uploaded to ~/grokbot-sp/p4 on both machines; binaries = symlinks to ../ex2gb/ex2 patched, -gb_wrap_hash_factor 2). helper: phase4/w.sh <script> -> wsl one-liner (base64).
## Settings
- Strong (ss_*.pbs): fixed global L=pi*192, NX=49152 cells x 200 v = 9.83M particles, dt .05, 100 steps, v3 opts, fix=1. Aurora 1,2,4,8,16 nodes x96ppn (subdomain width at 1536 ranks 0.39 > |v|max dt 0.3); Polaris 1,2,4,8 x32ppn. timeout 240/300 s per run. Ref = 1 node (np1 would take ~45 min).
- Long-run (lr_*.pbs, debug 1 node, 45 min, timeout 2500 s each, 3 concurrent 1-rank runs, phase-2 physics opts v2, x in [0,4pi]):
  x10_dt05 (256x200=51,200 particles, 6000 steps, t=300), x10_dt01 (6000 steps, t=60), x100_dt05_lowres (64x50=3,200 particles, 60000 steps, t=3000; full res would be ~3.6 h at 0.22 s/step).
- Weak 32/64 (ws_aurora_32/64.pbs): same per-rank load as phase 3 (256 cells/rank, L=pi*P), 3072 / 6144 ranks, 100 steps; ws32 also re-runs a 1-node np96 reference. timeout 600 / 780 s.
- Pusher long runs: NOT on ALCF (serial C code; built on box only). Run on box: scripts/box_pusher_long.sh, 1e8 steps (100x phase 2) uniform+solovev x boris/rk4/dg -> phase4/box_pusher/ (dg uses -ts_max_snes_failures -1; first 1e7 attempt w/o it aborted DIVERGED_NONLINEAR_SOLVE).
## Log
- 03:35 ET Polaris: SUBMITTED 7738512 ss_polaris.pbs (debug-scaling, 8 nodes, 15 min) Q 'Not enough free nodes'; 7738513 lr_polaris.pbs (debug, 1 node, 45 min) started 03:37:51 ET.
- 03:37 ET Aurora: SUBMITTED 8922402 ss_aurora.pbs (debug-scaling 16n 15 min) started 03:39:24 ET; 8922403 lr_aurora.pbs (debug 1n 45 min) started 03:38:56 ET.
  prod queue REJECTED 32/64-node jobs: Aurora prod resources_min.nodect=256. Moved to debug-scaling (max 256 nodes, max_run 1/user, 1 queued job/user).
  8922409 (ws32) queued, ws64 rejected by 1-queued limit (also with -W depend). qdel 8922409 (never ran). RESUBMITTED ws32 as 8922414; it qsubs ws_aurora_64.pbs from its first lines once running (status line 'chained ws64 submit: <id>' in ~/grokbot-status.txt).
- 03:41 ET health check: all running jobs producing E: lines; rates ok (Aurora x10 ~22 min, x100 lowres ~14 min; Polaris x10 ~31 min). "PETSc installed without X Windows" ERROR line in outputs is the -draw_type null warning, runs continue.
## Node-hour worst case (walltime x nodes)
ssA 16x0.25=4.0, ssP 8x0.25=2.0, lrA 0.75, lrP 0.75, ws32 32x0.25=8.0, ws64 64x0.30=19.2 -> 34.7 nh (<= 40). Expected actual ~20-25.
## Collect later
On each machine: bash ~/grokbot-sp/p4/collect.sh -> ~/grokbot-sp/p4/p4_results_<host>.tgz + prints summaries/status. Quick check cmd: phase4/cmd_check.txt (run on laptop machineId f08c...).

ALCF job scripts used by Grok Bot (copies; originals live in ~/grokbot-sp/ on each machine).
common_opts_v1.txt : ex2 options for 600-step Landau runs (GMRES+GAMG, ksp atol 1e-14)  -- jobs 7737974, 8921878
common_opts_v2.txt : CG+GAMG, atol 1e-12 (jobs 7738002, 8921903)
pilot4 (final weak-scaling) = v2 + "-em_ksp_type gmres -em_ksp_rtol 1e-8 -em_ksp_atol 1e-9 -ts_time_step 0.01 -ts_max_time 1.0", 100 steps
Builds: Polaris ~/grokbot-sp/build_petsc.sh (PrgEnv-gnu, cuda/13.0, craype-accel-nvidia80, darshan/xalt unloaded, CUDA, opt)
        Aurora  ~/grokbot-sp/build_petsc_aurora.sh (below)

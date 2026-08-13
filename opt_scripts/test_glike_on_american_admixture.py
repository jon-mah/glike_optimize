"""
Benchmark: original glike vs glike variants on the American Admixture model.

Reproduces the American Admixture experiment from the glike paper,
comparing the original optimizer against glike variants (reparameterized
optimizer + state pruning + tree-level multiprocessing).

Usage:
  # Step 1: Simulate trees (SLURM array 0-19, one per replicate)
  python benchmark_american_admixture.py simulate <idx>

  # Step 2: Estimate parameters (SLURM array 0-399, 20 reps x 20 init sets;
  #         idx = rep * N_INIT_SETS + set_idx)
  python benchmark_american_admixture.py estimate <idx>

  # Step 3: Plot results (run once after all estimation jobs finish)
  python benchmark_american_admixture.py plot
"""

import logging
import os
import glob
import math
import sys
import json
import time
import re
import ast
import pdb
import glike
import tskit
import msprime
import numpy as np
import importlib.util


# ===================== Helpers =====================

def _config_suffix():
    """Return '' when defaults, else '_reps{N_REPLICATES}-trees{N_TREES}'.

    Also append '_msprime{version}' if msprime version != 1.3.4, and
    '_haps{N_SAMPLES}' if N_SAMPLES != 1000.
    """
    # if N_REPLICATES != 20 or N_TREES != 20:
    suffix = f"_reps{N_REPLICATES}-trees{N_TREES}"
    if msprime.__version__ != "1.3.4":
        suffix += f"_msprime{msprime.__version__}"
    if N_SAMPLES != 1000:
        suffix += f"_haps{N_SAMPLES}"
    return suffix
    #return ""


def _results_suffix():
    """Suffix for the results dir and plot filenames.

    Extends _config_suffix() with the init-sets count and the init mode.
    Results and plots depend on these, but the simulated trees do not —
    so this is kept separate from _config_suffix().
    """
    if INIT_MODE == 'fixed-init-seeded':
        return _config_suffix() + f"_{N_INIT_SETS}-seeds-fixed-init"
    return _config_suffix() + f"_{N_INIT_SETS}-sets-of-init-vals"


# ===================== Random initial values =====================

def _sample_time_chain(rng, time_params, bounds, param_names):
    """Sample sorted time values respecting per-param numeric bounds and chain order.

    For each time param, the effective lower bound is the prefix-max of numeric
    lowers (over earlier times in the chain) and the effective upper is the
    suffix-min of numeric uppers. Values are then sampled sequentially with a
    tiny multiplicative buffer to enforce strict ordering t_1 < t_2 < ... < t_N.
    """
    if not time_params:
        return {}

    n = len(time_params)
    numeric_bounds = []
    for p in time_params:
        idx = param_names.index(p)
        lo, hi = bounds[idx]
        numeric_bounds.append((
            lo if isinstance(lo, (int, float)) else None,
            hi if isinstance(hi, (int, float)) else None,
        ))

    eff_lo = [None] * n
    eff_hi = [None] * n
    cur = None
    for i in range(n):
        nlo = numeric_bounds[i][0]
        if nlo is not None:
            cur = nlo if cur is None else max(cur, nlo)
        eff_lo[i] = cur
    cur = None
    for i in range(n - 1, -1, -1):
        nhi = numeric_bounds[i][1]
        if nhi is not None:
            cur = nhi if cur is None else min(cur, nhi)
        eff_hi[i] = cur

    eff_lo = [v if v is not None else 1.0 for v in eff_lo]
    eff_hi = [v if v is not None else 1e5 for v in eff_hi]

    out = {}
    prev = None
    for i, p in enumerate(time_params):
        lo = eff_lo[i] if prev is None else max(eff_lo[i], prev * (1.0 + 1e-9))
        hi = eff_hi[i]
        v = float(rng.uniform(lo, hi)) if lo < hi else float(lo)
        out[p] = v
        prev = v
    return out


def _sample_proportions(rng, prop_params, bounds, param_names):
    """Sample proportion params. If their bounds cross-reference each other
    (e.g. r1 upper = "1-r2"), sample from the simplex via Dirichlet over n+1
    dims; otherwise sample uniformly per-param within numeric bounds.
    """
    if not prop_params:
        return {}

    numeric_bounds = []
    cross_ref = False
    for p in prop_params:
        idx = param_names.index(p)
        lo, hi = bounds[idx]
        for endpoint in (lo, hi):
            if isinstance(endpoint, str):
                ref = endpoint[2:] if endpoint.startswith("1-") else endpoint
                if ref in prop_params:
                    cross_ref = True
        numeric_bounds.append((
            lo if isinstance(lo, (int, float)) else 0.0,
            hi if isinstance(hi, (int, float)) else 1.0,
        ))

    out = {}
    if cross_ref:
        sample = rng.dirichlet(np.ones(len(prop_params) + 1))
        for p, v, (lo_n, hi_n) in zip(prop_params, sample[:-1], numeric_bounds):
            out[p] = float(min(max(float(v), lo_n), hi_n))
    else:
        for p, (lo_n, hi_n) in zip(prop_params, numeric_bounds):
            out[p] = float(rng.uniform(lo_n, hi_n))
    return out


def _sample_other(rng, params, bounds, param_names):
    """Sample remaining params uniformly within numeric bounds."""
    out = {}
    for p in params:
        idx = param_names.index(p)
        lo, hi = bounds[idx]
        if not isinstance(lo, (int, float)) or not isinstance(hi, (int, float)):
            raise ValueError(f"non-numeric bounds for {p}: ({lo}, {hi})")
        out[p] = float(rng.uniform(lo, hi))
    return out


def _generate_random_init(seed, bounds, param_names, param_types):
    """Generate one random initial-values dict, deterministic from `seed`."""
    rng = np.random.default_rng(seed)
    time_params = sorted(
        [p for p in param_names if param_types.get(p, "").startswith("time:")],
        key=lambda p: int(param_types[p].split(":")[1]),
    )
    prop_params = [p for p in param_names if param_types.get(p) == "proportion"]
    other_params = [p for p in param_names if p not in time_params and p not in prop_params]

    out = {}
    out.update(_sample_time_chain(rng, time_params, bounds, param_names))
    out.update(_sample_proportions(rng, prop_params, bounds, param_names))
    out.update(_sample_other(rng, other_params, bounds, param_names))
    return {p: out[p] for p in param_names}


def _get_or_create_init_values(rep, set_idx, results_dir, bounds, param_names,
                                param_types, test, fixed_x0=None):
    """Return (x0, seed) for (rep, set_idx); generate-and-persist on first use.

    File: {results_dir}/init_values/init_{test}_rep{rep}_set{set_idx}.json.
    The TEST tag is included because different TESTs estimate different
    parameter sets (e.g. fix-r1-r2 drops r1/r2), so init values are not
    interchangeable across TESTs.

    If fixed_x0 is provided (fixed-init-seeded mode), all sets share that x0
    and only the seed varies (seed = rep * 10000 + set_idx). JSON body is
    {"x0": ..., "seed": ...} and the returned seed is non-None — the caller
    should pass it to np.random.seed() before maximize().

    If fixed_x0 is None (random mode), x0 is randomly generated within bounds
    and the JSON body is the x0 dict directly. The returned seed is None.

    Delete the file to force regeneration on the next run.
    """
    init_dir = os.path.join(results_dir, "init_values")
    os.makedirs(init_dir, exist_ok=True)
    path = os.path.join(init_dir, f"init_{test}_rep{rep}_set{set_idx}.json")

    if fixed_x0 is not None:
        if os.path.exists(path):
            with open(path) as f:
                data = json.load(f)
            return data["x0"], data["seed"]
        seed = rep * 10000 + set_idx
        x0 = dict(fixed_x0)
        with open(path, "w") as f:
            json.dump({"x0": x0, "seed": seed}, f, indent=2)
        return x0, seed

    if os.path.exists(path):
        with open(path) as f:
            return json.load(f), None
    seed = rep * 10000 + set_idx
    x0 = _generate_random_init(seed, bounds, param_names, param_types)
    with open(path, "w") as f:
        json.dump(x0, f, indent=2)
    return x0, None


def _find_best_set_per_rep(prefix, test, results_dir, n_replicates):
    """For each replicate, return the set with the highest estimated logp.

    Returns dict: rep -> (set_idx, est_dict, est_logp, elapsed_seconds).
    Replicates with no parseable result files are absent from the dict.
    """
    out = {}
    for rep in range(n_replicates):
        pattern = os.path.join(results_dir, f"{prefix}_{test}_rep{rep}_set*.txt")
        best = None
        for fpath in sorted(glob.glob(pattern)):
            m = re.search(r"_set(\d+)\.txt$", fpath)
            if m is None:
                continue
            set_idx = int(m.group(1))
            try:
                with open(fpath) as f:
                    lines = f.readlines()
                parts = lines[3].split("\t")
                est_dict = ast.literal_eval(parts[0])
                est_logp = float(parts[1].strip())
                elapsed = float(lines[5].strip())
            except (IndexError, ValueError, SyntaxError):
                continue
            if best is None or est_logp > best[2]:
                best = (set_idx, est_dict, est_logp, elapsed)
        if best is not None:
            out[rep] = best
    return out


# ===================== Logging =====================

def get_logger(name, path=None):
    logger = logging.getLogger(name)
    if not logger.handlers:
        logger.propagate = 0
        console = logging.StreamHandler()
        logger.addHandler(console)
        log_format = "[%(asctime)s - %(levelname)s] %(message)s"
        date_format = "%Y-%m-%d %H:%M:%S"
        formatter = logging.Formatter(fmt=log_format, datefmt=date_format)
        console.setFormatter(formatter)
        if path is not None:
            if not path.endswith('.log'):
                path += '.log'
            disk_log_stream = open(path, "w")
            disk_handler = logging.StreamHandler(disk_log_stream)
            logger.addHandler(disk_handler)
            disk_handler.setFormatter(formatter)
    return logger


# ===================== Simulate =====================

def simulate(idx):
    """Simulate ARG and save trees for one replicate.

    Args:
        idx: Replicate index (0-indexed, from SLURM_ARRAY_TASK_ID).
    """
    trees_dir = os.path.join(OUTPUT_DIR, f"trees{_config_suffix()}")
    os.makedirs(trees_dir, exist_ok=True)
    rep = idx
    seed = rep + 1

    log = get_logger(f"sim_{rep}", os.path.join(trees_dir, f"simulate_rep{rep}.log"))
    log.setLevel(logging.INFO)

    log.info(f"Simulating replicate {rep} (seed={seed}) ...")

    demography = glike.american_admixture_demography(**TRUE_PARAMS)
    arg = msprime.sim_ancestry(
        {"admix": N_SAMPLES},
        sequence_length=SEQUENCE_LENGTH,
        recombination_rate=RECOMBINATION_RATE,
        demography=demography,
        ploidy=1,
        random_seed=seed,
    )

    # Save ARG
    arg_path = os.path.join(trees_dir, f"true_{rep}.trees")
    arg.dump(arg_path)

    # Compute and save tree positions
    step = int(SEQUENCE_LENGTH) // (N_TREES + 1)
    positions = list(range(step, int(SEQUENCE_LENGTH), step))[:N_TREES]
    positions_path = os.path.join(trees_dir, f"positions_{rep}.json")
    with open(positions_path, "w") as f:
        json.dump(positions, f)

    log.info(f"Saved ARG ({arg.num_trees} trees) to {arg_path}")
    log.info(f"Extracted {len(positions)} tree positions to {positions_path}")
    log.info(f"Simulation for replicate {rep} complete.")


# ===================== Estimate =====================

def _parse_trajectory_for_resume(trajectory_path):
    """Parse a trajectory file and decide whether resume is applicable.

    Returns (last_epoch_num, last_x, last_lrs) on success.
    Returns None if the file says the previous run stopped early (resume not
    applicable) or if the file is malformed (parse failure).
    """
    with open(trajectory_path) as f:
        lines = f.readlines()

    # Find the last status line; if it says "stopped early", resume is not
    # applicable.
    status_line = None
    for line in reversed(lines):
        if line.startswith("# status:"):
            status_line = line.strip()
            break
    if status_line is None:
        print(f"trajectory file has no status line: {trajectory_path}", flush=True)
        return None
    if "stopped early" in status_line:
        print(f"resume not applicable: {status_line}", flush=True)
        return None

    # Parse the most recent "# epoch N end: {x} y" and "# epoch N lrs: {lrs}".
    end_re = re.compile(r"^# epoch (\d+) end: (\{.*\}) (\S+)\s*$")
    lrs_re = re.compile(r"^# epoch (\d+) lrs: (\{.*\})\s*$")
    last_epoch_num = None
    last_x = None
    last_lrs = None
    for line in reversed(lines):
        line = line.rstrip("\n")
        if last_lrs is None:
            m = lrs_re.match(line)
            if m is not None:
                last_lrs = ast.literal_eval(m.group(2))
                continue
        if last_x is None:
            m = end_re.match(line)
            if m is not None:
                last_epoch_num = int(m.group(1))
                last_x = ast.literal_eval(m.group(2))
                continue
        if last_x is not None and last_lrs is not None:
            break

    if last_x is None or last_lrs is None or last_epoch_num is None:
        print(f"failed to parse last epoch from {trajectory_path}", flush=True)
        return None
    return last_epoch_num, last_x, last_lrs


def estimate(idx, resume=False):
    """Run estimation for one (replicate, method) pair.

    Job layout: 20 replicates x 2 methods = 40 jobs.
        idx 0  -> (rep=0,  method='original')
        idx 1  -> (rep=0,  method='glike_ji_v1')
        idx 2  -> (rep=1,  method='original')
        ...
        idx 39 -> (rep=19, method='glike_ji_v1')

    Args:
        idx: Job index (0-indexed, from SLURM_ARRAY_TASK_ID).
        resume: If True, continue from the last epoch in the trajectory file.
            The trajectory file's last x and learning rates are used as the
            new initial state, and epoch numbering continues without a gap.
            If the previous run already stopped early, this returns without
            running anything (resume not applicable).
            If False (default), fail-fast when the trajectory file already
            exists, so previous results are not silently overwritten.
    """
    results_dir = os.path.join(OUTPUT_DIR, f"results{_results_suffix()}")
    if TREE_PRUNE != 0:
        results_dir += f"_prune{TREE_PRUNE}"
    os.makedirs(results_dir, exist_ok=True)

    # methods = ['original', 'glike_ji_v1']
    # methods = ['glike_ji_v1']
    # methods = ['glike_ji_v2']
    # methods = ['glike_ji_v3']
    # methods = ['glike_ji_v4']
    # idx layout: rep = idx // N_INIT_SETS, set_idx = idx % N_INIT_SETS
    rep = idx // N_INIT_SETS
    set_idx = idx % N_INIT_SETS
    # pdb.set_trace()
    method = 'glike_ji_only_tree_parallel'

    if 'glike_ji' in method:
        _glike_ji_module = 'glike_ji_v1' if method == 'glike_ji_only_tree_parallel' else method
        _glike_ji_path = f'/project2/jitang_1167/projects/{_glike_ji_module}/glike'
        _spec = importlib.util.spec_from_file_location(
            _glike_ji_module,
            os.path.join(_glike_ji_path, "__init__.py"),
            submodule_search_locations=[_glike_ji_path])
        glike_ji = importlib.util.module_from_spec(_spec)
        sys.modules[_glike_ji_module] = glike_ji
        _spec.loader.exec_module(glike_ji)

    outfile = os.path.join(results_dir, f"{method}_{TEST}_rep{rep}_set{set_idx}.txt")
    log = get_logger(f"est_{method}_{TEST}_{rep}_set{set_idx}",
                     os.path.join(results_dir, f"{method}_{TEST}_rep{rep}_set{set_idx}.log"))
    log.setLevel(logging.INFO)
    log.info(f"Replicate {rep}, set {set_idx}, method={method}, idx={idx}, resume={resume}")

    # ---- Trajectory file path & resume handling ----
    trajectory_path = os.path.join(
        results_dir, f"{method}_{TEST}_rep{rep}_set{set_idx}_logp-trajectory.txt")
    start_epoch = 0
    lrs0 = None
    fixed_x0 = X0 if INIT_MODE == 'fixed-init-seeded' else None
    x0_to_use, seed_to_use = _get_or_create_init_values(
        rep, set_idx, results_dir, BOUNDS, PARAM_NAMES, PARAM_TYPES, TEST,
        fixed_x0=fixed_x0)
    log.info(f"INIT_MODE={INIT_MODE}; initial values for (rep={rep}, set={set_idx}): {x0_to_use}")
    if seed_to_use is not None:
        log.info(f"Setting np.random.seed({seed_to_use}) before maximize()")
        np.random.seed(seed_to_use)
    if resume:
        if not os.path.exists(trajectory_path):
            log.error(f"resume requested but trajectory file not found: {trajectory_path}")
            return
        parsed = _parse_trajectory_for_resume(trajectory_path)
        if parsed is None:
            log.info("resume not applicable; exiting without running.")
            return
        last_epoch_num, last_x, last_lrs = parsed
        start_epoch = last_epoch_num + 1
        lrs0 = last_lrs
        x0_to_use = last_x
        log.info(f"resuming from epoch {start_epoch} (previous last epoch = {last_epoch_num})")
    else:
        if os.path.exists(trajectory_path):
            log.error(f"trajectory file already exists: {trajectory_path}")
            log.error("either re-run with resume=True (4th CLI arg = 'resume') "
                      "or delete the file to start fresh.")
            return

    # ---- Load trees ----
    trees_dir = os.path.join(OUTPUT_DIR, f"trees{_config_suffix()}")
    arg_path = os.path.join(trees_dir, f"true_{rep}.trees")
    positions_path = os.path.join(trees_dir, f"positions_{rep}.json")
    if not os.path.exists(arg_path):
        log.error(f"ARG file not found: {arg_path}. Run 'simulate' first.")
        return
    arg = tskit.load(arg_path)
    with open(positions_path) as f:
        positions = json.load(f)
    trees = [arg.at(pos).copy() for pos in positions]
    log.info(f"Loaded {len(trees)} trees from {arg_path}")

    # ---- Sample-to-population mapping ----
    samples = {i: "admix" for i in range(N_SAMPLES)}

    # ---- Define objective function and run optimizer ----
    # pdb.set_trace()
    if method == 'original': # the original glike
        def fun(**kwargs):
            try:
                demo = glike.american_admixture_demo(**kwargs)
                return glike.glike_trees(trees, demo, samples=samples, kappa=KAPPA, spread=spread, prune=TREE_PRUNE)
            except Exception:
                return -math.inf
        
        # ---- Log-likelihood at true parameters ----
        demo_true = glike.american_admixture_demo(**TRUE_PARAMS)
        logp_true = glike.glike_trees(trees, demo_true, samples=samples, kappa=KAPPA, spread=spread, prune=TREE_PRUNE)
        log.info(f"Log-likelihood at true parameters: {logp_true}")

        log.info("Running glike.maximize() (original) ...")
        t_start = time.time()
        x, logp = glike.maximize(fun, x0_to_use, bounds=BOUNDS, epochs=EPOCHS, verbose=True)
        elapsed = time.time() - t_start

    elif method == 'glike_ji_only_tree_parallel':
        def fun(**kwargs):
            try:
                demo = glike_ji.american_admixture_demo(**kwargs)
                return glike_ji.glike_trees(trees, demo, samples=samples, kappa=KAPPA, spread=spread, prune=TREE_PRUNE,
                                         state_prune=math.inf, n_workers=n_workers)
            except Exception:
                return -math.inf

        n_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
        if N_WORKERS_V1 == -1:
            log.info(f"Using all available CPUs: {n_workers}")
        else:
            assert N_WORKERS_V1 <= n_workers, f"Requested {N_WORKERS_V1} workers, but only {n_workers} available"
            n_workers = N_WORKERS_V1
            log.info(f"Using requested CPUs: {n_workers}")

        # ---- Log-likelihood at true parameters ----
        demo_true = glike_ji.american_admixture_demo(**TRUE_PARAMS)
        # pdb.set_trace()
        logp_true = glike_ji.glike_trees(trees, demo_true, samples=samples, kappa=KAPPA, spread=spread, prune=TREE_PRUNE,
                                            state_prune=math.inf, n_workers=n_workers)
        log.info(f"Log-likelihood at true parameters: {logp_true}")

        log.info(f"Running glike_ji.maximize() (only tree-level multiprocessing, "
                 f"state_prune=inf, n_workers={n_workers}) ...")

        t_start = time.time()
        x, logp = glike_ji.maximize(fun, x0_to_use, bounds=BOUNDS, epochs=EPOCHS, verbose=True,
                                    log_file=trajectory_path,
                                    start_epoch=start_epoch, lrs0=lrs0)
        elapsed = time.time() - t_start

    elif 'glike_ji' in method:
        def fun(**kwargs):
            try:
                demo = glike_ji.american_admixture_demo(**kwargs)
                return glike_ji.glike_trees(trees, demo, samples=samples, kappa=KAPPA, spread=spread, prune=TREE_PRUNE,
                                         state_prune=STATE_PRUNE, n_workers=n_workers)
            except Exception:
                return -math.inf

        n_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
        if N_WORKERS_V1 == -1:
            log.info(f"Using all available CPUs: {n_workers}")
        else:
            assert N_WORKERS_V1 <= n_workers, f"Requested {N_WORKERS_V1} workers, but only {n_workers} available"
            n_workers = N_WORKERS_V1
            log.info(f"Using requested CPUs: {n_workers}")
        
        # ---- Log-likelihood at true parameters ----
        demo_true = glike_ji.american_admixture_demo(**TRUE_PARAMS)
        logp_true = glike_ji.glike_trees(trees, demo_true, samples=samples, kappa=KAPPA, spread=spread, prune=TREE_PRUNE,
                                            state_prune=STATE_PRUNE, n_workers=n_workers)
        log.info(f"Log-likelihood at true parameters: {logp_true}")

        log.info(f"Running {method}, "
                 f"state_prune={STATE_PRUNE}, n_workers={n_workers} ...")
            
        t_start = time.time()
        if method in ['glike_ji_v4', 'glike_ji_v6']:
            freeze_graph=False
            if method == 'glike_ji_v6':
                freeze_graph=True
            x, logp = glike_ji.maximize_lbfgs_fd(
                                            model_fn=glike_ji.american_admixture_demo,
                                            trees=trees,
                                            x0=x0_to_use,
                                            param_types=PARAM_TYPES,
                                            samples=samples,
                                            kappa=KAPPA,
                                            spread=spread,
                                            state_prune=STATE_PRUNE,     
                                            n_workers=n_workers, 
                                            maxiter=50,
                                            verbose=True,
                                            freeze_graph=freeze_graph
                                        )
        elif method == 'glike_ji_v5':
            x, logp = glike_ji.maximize_lbfgs(
                                            model_fn=glike_ji.american_admixture_demo,
                                            trees=trees,
                                            x0=x0_to_use,
                                            param_types=PARAM_TYPES,
                                            samples=samples,
                                            kappa=KAPPA,
                                            spread=spread,
                                            state_prune=STATE_PRUNE,     
                                            n_workers=n_workers, 
                                            maxiter=50,
                                            verbose=True,
                                        )
        else:
            x, logp = glike_ji.maximize_reparam(fun, x0_to_use, PARAM_TYPES, epochs=EPOCHS, verbose=True)
        elapsed = time.time() - t_start

    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)
    log.info(f"Finished in {hours}h {minutes}m {elapsed % 60:.1f}s")
    log.info(f"Estimated: {x}")
    log.info(f"logp = {logp}  (true = {logp_true})")

    # ---- Write output (same format as run_sglike_nh_sim.py) ----
    with open(outfile, "w") as f:
        f.write("# True parameters:\n")
        f.write(str(TRUE_PARAMS) + "\t" + str(logp_true) + "\n")
        f.write("# Estimated parameters:\n")
        f.write(str(x) + "\t" + str(logp) + "\n")
        f.write("# Elapsed time (seconds):\n")
        f.write(str(elapsed) + "\n")

    log.info(f"Results written to {outfile}")


# ===================== Plot =====================

def plot_results():
    """Plot relative-error boxplots for both methods, side by side and individually."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib.transforms as mtransforms
    import pandas as pd

    # METHODS = ['glike', 'glike_ji_v1']  # 'glike' and/or 'glike_ji_v1'
    # METHODS = ['glike_ji_v1']
    # METHODS = ['glike_ji_v2']
    # METHODS = ['glike_ji_v3']
    # METHODS = ['glike']
    METHODS = ['glike_ji_only_tree_parallel']

    suffix = _results_suffix()
    if TREE_PRUNE != 0:
        suffix += f"_prune{TREE_PRUNE}"
    results_dir = os.path.join(OUTPUT_DIR, f"results{suffix}")
    
    # ---- X-axis labels and colors, keyed by parameter name so they stay in
    # sync with PARAM_NAMES under configs like 'fix-r1-r2' that drop params. ----
    _GROUP_COLORS = {
        "t":  (205/255, 133/255,  63/255),  # orange
        "r":  ( 16/255,  78/255, 139/255),  # blue
        "N":  (  0/255, 139/255,   0/255),  # green
        "gr": (204/255, 121/255, 167/255),  # pink
    }

    def _label_for(p):
        if p.startswith("gr_"):
            return "$\\mathregular{gr_{" + p[3:] + "}}$"
        if p.startswith("N_"):
            return "$\\mathregular{N_{" + p[2:] + "}}$"
        if p.startswith("t") or p.startswith("r"):
            return "$\\mathregular{" + p[0] + "_{" + p[1:] + "}}$"
        return p

    def _color_for(p):
        if p.startswith("gr_"):
            return _GROUP_COLORS["gr"]
        if p.startswith("N_"):
            return _GROUP_COLORS["N"]
        if p.startswith("r"):
            return _GROUP_COLORS["r"]
        return _GROUP_COLORS["t"]

    names = [_label_for(p) for p in PARAM_NAMES]
    colors = [_color_for(p) for p in PARAM_NAMES]

    n_params = len(PARAM_NAMES)
    truth = [TRUE_PARAMS[p] for p in PARAM_NAMES]
    truth_series = pd.Series(truth, index=PARAM_NAMES)

    # ---- Method name mappings ----
    file_prefix = {'glike': 'original', 'glike_ji_v1': 'glike_ji_v1', 'glike_ji_v2': 'glike_ji_v2', 'glike_ji_v3': 'glike_ji_v3', 'glike_ji_only_tree_parallel': 'glike_ji_only_tree_parallel'}
    display_label = {'glike': 'glike (original)', 'glike_ji_v1': 'glike-ji-v1', 'glike_ji_v2': 'glike-ji-v2', 'glike_ji_v3': 'glike-ji-v3', 'glike_ji_only_tree_parallel': 'glike-ji (only tree parallel)'}

    # ---- Read results for each method (best set per replicate) ----
    method_data = {}
    for method in METHODS:
        prefix = file_prefix[method]
        best_per_rep = _find_best_set_per_rep(prefix, TEST, results_dir, N_REPLICATES)

        all_estimated = []
        all_times = []
        n_skipped_nonfinite = 0
        for rep in range(N_REPLICATES):
            if rep not in best_per_rep:
                print(f"[{method}] rep {rep}: no result files found")
                continue
            set_idx, est_dict, est_logp, elapsed = best_per_rep[rep]
            if not math.isfinite(est_logp):
                print(f"[{method}] rep {rep}: best logp = {est_logp}, skipping (all sets failed)")
                n_skipped_nonfinite += 1
                continue
            all_estimated.append([est_dict[p] for p in PARAM_NAMES])
            all_times.append(elapsed)
            print(f"[{method}] rep {rep}: best set={set_idx}, logp={est_logp:.2f}")

        method_data[method] = {
            'estimated': all_estimated,
            'times': all_times,
        }
        n_loaded = len(all_estimated)
        if n_loaded > 0:
            extra = f" ({n_skipped_nonfinite} skipped: best logp non-finite)" if n_skipped_nonfinite else ""
            print(f"[{method}] Loaded best-of-{N_INIT_SETS} for {n_loaded}/{N_REPLICATES} replicates{extra}. "
                  f"Mean time (best set): {np.mean(all_times)/60:.1f} min")
        else:
            print(f"[{method}] No result files found.")

    # ---- Check we have data ----
    if all(not method_data[m]['estimated'] for m in METHODS):
        print("No result files found. Nothing to plot.")
        return

    # ---- Compute relative errors ----
    def compute_errors(all_estimated):
        data = pd.DataFrame(all_estimated, columns=PARAM_NAMES)
        zero_params = truth_series[truth_series == 0].index.tolist()
        nonzero_params = truth_series[truth_series != 0].index.tolist()
        errors = pd.DataFrame(index=data.index, columns=PARAM_NAMES, dtype=float)
        if nonzero_params:
            errors[nonzero_params] = data[nonzero_params].divide(truth_series[nonzero_params]) - 1
        if zero_params:
            errors[zero_params] = data[zero_params]
        return errors

    # ---- Helper: draw one boxplot panel ----
    def draw_panel(ax, errors, method_label, mean_time_min, n_loaded):
        errors = errors.copy()
        errors[errors > 3] = 3

        n = errors.shape[0]

        # Reference lines
        ax.axhline(y=0, linestyle="solid", color="black", zorder=1)
        ax.axhline(y=-0.3, linestyle="dotted", color="black", zorder=1)
        ax.axhline(y=0.3, linestyle="dotted", color="black", zorder=1)

        # Scatter individual points with jitter
        for i in range(n_params):
            jitter = [i + 1 - 0.2 + 0.4 * np.random.uniform() for _ in range(n)]
            ax.scatter(jitter, errors.iloc[:, i], s=5, facecolors='none',
                       edgecolor=colors[i], linewidths=0.5, zorder=2)

        # Box plot
        def adjust_color(color, alpha):
            if len(color) == 3:
                color = (color[0], color[1], color[2], 1.0)
            return (color[0], color[1], color[2], color[3] * alpha)

        bp = ax.boxplot(errors.values, showfliers=False, showcaps=False,
                        zorder=2, widths=0.75, patch_artist=True)
        for i, color in enumerate(colors):
            bp["boxes"][i].set_edgecolor(color)
            bp["boxes"][i].set_facecolor(adjust_color(color, alpha=0.2))
            bp["boxes"][i].set_linewidth(1.5)
            bp["whiskers"][i * 2].set_color(color)
            bp["whiskers"][i * 2].set_linewidth(1.5)
            bp["whiskers"][i * 2 + 1].set_color(color)
            bp["whiskers"][i * 2 + 1].set_linewidth(1.5)
            bp["medians"][i].set_color(color)
            xm, ym = bp["medians"][i].get_data()
            xn = (xm - (xm.sum() / 2.)) * 1.0 + (xm.sum() / 2.)
            ax.hlines(y=ym, xmin=xn[0], xmax=xn[1], color=color,
                      linewidth=2.5, zorder=4)

        ax.set_xticks(list(range(1, n_params + 1)))
        ax.set_xticklabels(names, rotation=-30, ha="left", fontsize=8)
        ax.set_yticks([-1, 0, 1, 2, 3])
        ax.set_ylim(-1.2, 3.2)
        ax.set_ylabel('Relative error', fontsize=11)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Title with timing info
        if mean_time_min < 60:
            time_str = f"{mean_time_min:.1f} min"
        else:
            time_str = f"{mean_time_min/60:.1f} h"
        ax.set_title(f"{method_label}  (n={n_loaded}, mean time: {time_str})", fontsize=11)

    # ---- Side-by-side comparison (only when both methods are present) ----
    if len(METHODS) == 2:
        fig, axes = plt.subplots(1, 2, figsize=(24, 6), sharey=True)

        for i, method in enumerate(METHODS):
            label = display_label[method]
            ax = axes[i]
            md = method_data[method]
            if not md['estimated']:
                ax.text(0.5, 0.5, 'No results yet', transform=ax.transAxes,
                        ha='center', va='center', fontsize=14, color='gray')
                ax.set_title(f"{label}  (no data)", fontsize=11)
                continue
            errors = compute_errors(md['estimated'])
            mean_time = np.mean(md['times']) / 60
            draw_panel(ax, errors, label, mean_time, len(md['estimated']))

        dx = -7 / 72.
        offset = mtransforms.ScaledTranslation(dx, 0, fig.dpi_scale_trans)
        for ax in axes:
            for lbl in ax.xaxis.get_majorticklabels():
                lbl.set_transform(lbl.get_transform() + offset)

        plt.tight_layout()
        outfig = os.path.join(OUTPUT_DIR, f"benchmark_comparison{suffix}.png")
        fig.savefig(outfig, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {outfig}")
        plt.close(fig)

    # ---- Individual plots per method ----
    for method in METHODS:
        label = display_label[method]
        md = method_data[method]
        if not md['estimated']:
            continue
        errors = compute_errors(md['estimated'])
        mean_time = np.mean(md['times']) / 60

        fig_single, ax_single = plt.subplots(figsize=(14, 5))
        draw_panel(ax_single, errors, label, mean_time, len(md['estimated']))
        dx = -7 / 72.
        offset_s = mtransforms.ScaledTranslation(dx, 0, fig_single.dpi_scale_trans)
        for lbl in ax_single.xaxis.get_majorticklabels():
            lbl.set_transform(lbl.get_transform() + offset_s)
        plt.tight_layout()
        outfig_single = os.path.join(OUTPUT_DIR, f"benchmark_{method}{suffix}_{TEST}.png")
        fig_single.savefig(outfig_single, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {outfig_single}")
        plt.close(fig_single)


def plot_logp_trajectory():
    """Plot logp trajectory across epochs, one curve per replicate (best set).

    For each replicate, picks the init set with the highest final logp (same
    selection as plot_results), then reads the corresponding trajectory file
    {method}_{TEST}_rep{rep}_set{set_idx}_logp-trajectory.txt and parses lines
    of the form `# epoch N end: {x} <logp>`. Each curve plots end-of-epoch
    logp vs. epoch number. Curves may have different lengths because
    maximize() can early-stop before reaching EPOCHS. Replicates whose best
    logp is non-finite (e.g. all init sets failed) are skipped.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    suffix = _results_suffix()
    if TREE_PRUNE != 0:
        suffix += f"_prune{TREE_PRUNE}"
    results_dir = os.path.join(OUTPUT_DIR, f"results{suffix}")
    end_re = re.compile(r"^# epoch (\d+) end: \{.*\} (\S+)\s*$")

    # Logp values below this floor are numerical artifacts (the optimizer
    # occasionally lands in a catastrophically bad spot for a single epoch
    # before recovering). Drop them so a single -1e279 doesn't squash the
    # y-axis for the other curves.
    LOGP_FLOOR = -1e10

    method = 'glike_ji_only_tree_parallel'
    best_per_rep = _find_best_set_per_rep(method, TEST, results_dir, N_REPLICATES)

    fig, ax = plt.subplots(figsize=(10, 6))
    cmap = plt.get_cmap("tab20")
    n_plotted = 0
    n_stopped_early = 0
    n_skipped_nonfinite = 0
    n_outliers_total = 0
    for rep in range(N_REPLICATES):
        if rep not in best_per_rep:
            print(f"Warning: no result files for rep {rep}, skipping")
            continue
        set_idx, _, best_logp, _ = best_per_rep[rep]
        if not math.isfinite(best_logp):
            print(f"Warning: rep {rep} best logp = {best_logp}, skipping (all sets failed)")
            n_skipped_nonfinite += 1
            continue
        path = os.path.join(
            results_dir, f"{method}_{TEST}_rep{rep}_set{set_idx}_logp-trajectory.txt")
        if not os.path.exists(path):
            print(f"Warning: {path} not found, skipping replicate {rep}")
            continue
        epochs = []
        logps = []
        with open(path) as f:
            for line in f:
                m = end_re.match(line.rstrip("\n"))
                if m is None:
                    continue
                try:
                    epochs.append(int(m.group(1)))
                    logps.append(float(m.group(2)))
                except ValueError:
                    pass
        if not logps:
            print(f"Warning: no '# epoch N end' lines parsed from {path}")
            continue
        n_epochs = len(logps)
        if n_epochs < EPOCHS:
            n_stopped_early += 1

        # Drop epochs whose end-of-epoch logp is below the sanity floor.
        keep = [(e, lp) for e, lp in zip(epochs, logps) if lp >= LOGP_FLOOR]
        n_dropped = len(logps) - len(keep)
        if n_dropped > 0:
            n_outliers_total += n_dropped
            print(f"  rep {rep} set {set_idx}: dropped {n_dropped} epoch(s) with "
                  f"logp < {LOGP_FLOOR:.0e}")
        if not keep:
            print(f"Warning: all epochs filtered out for rep {rep}, skipping curve")
            continue
        epochs, logps = (list(x) for x in zip(*keep))

        ax.plot(epochs, logps,
                marker="o", markersize=3, linewidth=1,
                color=cmap(rep % 20), label=f"rep {rep} set {set_idx} ({n_epochs})")
        n_plotted += 1

    ax.set_xlabel("epoch")
    ax.set_ylabel("logp")
    ax.set_title(f"logp trajectory across epochs ({n_plotted} replicates)")
    note = f"{n_stopped_early}/{n_plotted} stopped earlier (< {EPOCHS} epochs)"
    if n_skipped_nonfinite:
        note += f"\n{n_skipped_nonfinite} replicate(s) skipped: best logp non-finite"
    if n_outliers_total:
        note += f"\n{n_outliers_total} epoch point(s) dropped: logp < {LOGP_FLOOR:.0e}"
    ax.text(0.02, 0.98, note,
            transform=ax.transAxes, va="top", ha="left",
            fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.7", alpha=0.9))
    ax.legend(loc="lower right", ncol=2, fontsize=8, title="rep set (n_epochs)")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    outpath = os.path.join(results_dir, "logp_trajectory.png")
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    print(f"Plot saved to {outpath}")
    plt.close(fig)


# ===================== CLI =====================

if __name__ == "__main__":
    # ===================== Configuration =====================
    # OUTPUT_DIR = '/project2/jitang_1167/projects/glike_sglike/glike_ji_v1_on_american_admixture'
    OUTPUT_DIR = '/project2/jitang_1167/projects/glike_sglike_data/glike_on_american_admixture'
    N_WORKERS_V1 = -1  # Using all available CPUs. Request 20 cpus so that the jobs can routed across all the nodes on Endeavour and Discover
    N_REPLICATES, N_TREES = 20, 20
    # N_INIT_SETS = 20  # number of random initial-value sets per replicate; SLURM array spans 0..(N_REPLICATES * N_INIT_SETS - 1)
    N_INIT_SETS = 1
    # Mode for the N_INIT_SETS sets per rep:
    #   'random'             - 20 random x0 sets per rep, each within bounds (np.random.seed not touched).
    #   'fixed-init-seeded'  - all sets share the same x0 (the script's X0 for the current TEST), but each
    #                          gets a distinct np.random seed before maximize() so the stochastic likelihood
    #                          (immigrate_stochastic when num_links > kappa) drives 20 different trajectories.
    INIT_MODE = 'fixed-init-seeded'
    # N_REPLICATES, N_TREES  = 50, 60
    # N_REPLICATES, N_TREES  = 50, 100 # 5 batchs, each with 20 cpus. Or 2 batchs, each with 50 cpus (if the jobs can still be parallelized)
    N_SAMPLES = 1000  # haplotypes from "admix"
    # N_SAMPLES = 500
    SEQUENCE_LENGTH = 3e7
    RECOMBINATION_RATE = 1e-8
    KAPPA = 10000
    TREE_PRUNE = 0
    if N_TREES > 80:
        TREE_PRUNE = 0.5
    STATE_PRUNE = 20
    # STATE_PRUNE = math.inf  # no pruning for main comparison
    EPOCHS = 20
    spread = 1.0

    # TESTs = ['fix-r1-r2', 'kappa50000', 'kappa100000', 'old-init-values']
    TEST = sys.argv[1]
    command = sys.argv[2]

    # init_values_name = 'github-init-values'
    init_values_name = 'old-init-values'

    # True parameters (from paper/american_admixture/simulate.sh)
    TRUE_PARAMS = {
        "t1": 12, "t2": 920, "t3": 2040, "t4": 5920,
        "r1": 0.17, "r2": 0.33,
        "N_afr": 14474,
        "N_eur": round(1032 * math.exp(0.0038 * 920), 1),
        "N_asia": round(554 * math.exp(0.0048 * 920), 1),
        "N_admix": round(30000 * math.exp(0.05 * 12), 1),
        "N_ooa": 1861, "N_anc": 7310,
        "gr_eur": 0.0038, "gr_asia": 0.0048, "gr_admix": 0.05,
    }

    if init_values_name == 'github-init-values':
        # Initial values (from paper/american_admixture/estimate.sh)
        X0 = {
            "t1": 5, "t2": 200, "t3": 500, "t4": 2000,
            "r1": 0.33, "r2": 0.33,
            "N_afr": 10000, "N_eur": 10000, "N_asia": 10000,
            "N_admix": 10000, "N_ooa": 10000, "N_anc": 10000,
            "gr_eur": 0.01, "gr_asia": 0.01, "gr_admix": 0.01,
        }
        
        
        # Bounds (from paper/american_admixture/estimate.sh)
        BOUNDS = [
            (1, "t2"), ("t1", "t3"), ("t2", "t4"), ("t3", 100000),
            (0.0, "1-r2"), (0.0, "1-r1"),
            (100, 1e6), (100, 1e6), (100, 1e6), (100, 1e6), (100, 1e6), (100, 1e6),
            (0.0001, 1), (0.0001, 1), (0.0001, 1),
        ]
    
    if init_values_name == 'old-init-values':
        X0 = {
            "t1": 10, "t2": 500, "t3": 1000, "t4": 3000,
            "r1": 0.33, "r2": 0.33,
            "N_afr": 10000, "N_eur": 50000, "N_asia": 50000,
            "N_admix": 50000, "N_ooa": 1000, "N_anc": 1000,
            "gr_eur": 0.01, "gr_asia": 0.01, "gr_admix": 0.1,
        }
        BOUNDS = [
        (1, 100), (100, "t3"), ("t2", "t4"), ("t3", 20000),
        (0.01, 0.49), (0.01, 0.49),
        (100, 1e6), (100, 1e6), (100, 1e6), (100, 1e6), (100, 1e6), (100, 1e6),
        (0.0001, 1), (0.0001, 1), (0.0001, 1),
        ]
    if TEST == 'kappa50000':
        KAPPA = 50000
    if TEST == 'kappa100000':
        KAPPA = 100000
    # Helper for "fix-*" TEST modes: drop the named params from X0/BOUNDS and
    # resolve any cross-references in the remaining bounds (e.g. ("t2", "t4")
    # or "1-r2") to numeric values from TRUE_PARAMS. Returns (new_X0, new_BOUNDS).
    def _apply_fix(fixed_params):
        def resolve(endpoint):
            if isinstance(endpoint, str):
                ref = endpoint[2:] if endpoint.startswith("1-") else endpoint
                if ref in fixed_params:
                    val = TRUE_PARAMS[ref]
                    return (1 - val) if endpoint.startswith("1-") else val
            return endpoint
        names = list(X0.keys())
        keep = [(n, b) for n, b in zip(names, BOUNDS) if n not in fixed_params]
        new_x0 = {n: X0[n] for n, _ in keep}
        new_bounds = [(resolve(lo), resolve(hi)) for _, (lo, hi) in keep]
        return new_x0, new_bounds

    # Parse any "fix-<p1>[-<p2>...]" TEST: extract the params to fix and
    # canonicalize TEST so that e.g. 'fix-t4-t3' and 'fix-t3-t4' resolve to
    # the same string (using the order in TRUE_PARAMS). Then apply the fix.
    if TEST.startswith('fix-'):
        _tokens = TEST[len('fix-'):].split('-')
        _all_params = list(TRUE_PARAMS.keys())
        for _tok in _tokens:
            if _tok not in _all_params:
                raise ValueError(f"unknown fix param '{_tok}' in TEST='{TEST}'")
        _fixed_params = sorted(set(_tokens), key=_all_params.index)
        TEST = 'fix-' + '-'.join(_fixed_params)
        X0, BOUNDS = _apply_fix(_fixed_params)
    # Parameter types for reparameterized optimizer (glike-ji-v1)
    PARAM_TYPES = {
        "t1": "time:0", "t2": "time:1", "t3": "time:2", "t4": "time:3",
        "r1": "proportion", "r2": "proportion",
        "N_afr": "size", "N_eur": "size", "N_asia": "size",
        "N_admix": "size", "N_ooa": "size", "N_anc": "size",
        "gr_eur": "positive", "gr_asia": "positive", "gr_admix": "positive",
    }

    # Parameter display order
    PARAM_NAMES = list(X0.keys())

    if len(sys.argv) < 2:
        print("Usage:")
        print("  python benchmark_american_admixture.py <TEST> simulate <idx>")
        print("  python benchmark_american_admixture.py <TEST> estimate <idx> [resume]")
        print("  python benchmark_american_admixture.py <TEST> plot")
        sys.exit(1)

    if command == "simulate":
        idx = int(sys.argv[3])
        simulate(idx)
    elif command == "estimate":
        # method = sys.argv[2]
        idx = int(sys.argv[3])
        resume = (len(sys.argv) > 4 and sys.argv[4].lower() == "resume")
        estimate(idx, resume=resume)
    elif command == "plot":
        plot_results()
        plot_logp_trajectory()
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)




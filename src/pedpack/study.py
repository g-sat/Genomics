"""End-to-end comparison on the block simulator.

The oracle table is computed so regret is measurable. Budgeted methods are
only allowed to read a subset of that table; the rest is an analysis device,
not part of the selector.
"""

from __future__ import annotations

import time
from collections import defaultdict

import numpy as np

from pedpack.candidates import build_pool
from pedpack.cost import storage_weight
from pedpack.genotype import build_cohort, mendelian_error_count
from pedpack.layout import prepare_layout
from pedpack.models import train_embeddings
from pedpack.pedigree import simulate_pedigree, trios
from pedpack.selector import (
    analytical_objective,
    measured_objective,
    pilot_count,
    predict_query_seconds,
    select_analytical,
    select_neural,
    select_random,
    stratified_pilot,
)
from pedpack.simulator import measure, measure_many
from pedpack.types import Cluster, Cohort, PreparedLayout, Query
from pedpack.workloads import make_workloads

WORKLOADS = ("fit", "scan", "selective", "heldout", "window")


def make_cluster(raw: dict | None) -> Cluster:
    if not raw:
        return Cluster()
    values = dict(raw)
    if "failed_nodes" in values:
        values["failed_nodes"] = tuple(int(x) for x in values["failed_nodes"])
    return Cluster().with_updates(**values)


def run_study(config: dict, verbose: bool = True) -> dict:
    started = time.perf_counter()
    seeds = [int(s) for s in config["seeds"]]
    rows = [run_seed(int(seed), config, verbose) for seed in seeds]
    report = {
        "config": config,
        "seeds": rows,
        "elapsed_seconds": time.perf_counter() - started,
    }
    return report


def run_seed(seed: int, config: dict, verbose: bool = True) -> dict:
    started = time.perf_counter()
    rng = np.random.default_rng(seed)
    cluster = make_cluster(config.get("cluster"))
    pedigree = simulate_pedigree(
        n_families=int(config["n_families"]),
        generations=int(config["generations"]),
        mean_children=float(config["mean_children"]),
        cross_rate=float(config["cross_rate"]),
        rng=rng,
        min_people=int(config.get("min_people", 24)),
    )
    cohort = build_cohort(
        pedigree,
        n_variants=int(config["n_variants"]),
        rng=rng,
        bp_spacing=int(config.get("bp_spacing", 100)),
    )
    errors = mendelian_error_count(cohort.genotypes, cohort.father, cohort.mother)
    if errors:
        raise RuntimeError(f"simulated genotypes have {errors} Mendelian errors")
    embeddings = train_embeddings(
        cohort,
        steps=int(config.get("encoder_steps", 120)),
        hidden=int(config.get("hidden", 32)),
        seed=seed,
    )
    layouts, heuristic_index = _layouts(cohort, embeddings, cluster, config)
    workloads = make_workloads(
        cohort,
        rng,
        n_fit_scan=int(config.get("n_fit_scan", 4)),
        n_fit_selective=int(config.get("n_fit_selective", 10)),
        n_eval_scan=int(config.get("n_eval_scan", 4)),
        n_eval_selective=int(config.get("n_eval_selective", 6)),
        n_heldout=int(config.get("n_heldout", 6)),
        n_windows=int(config.get("n_windows", 4)),
    )
    catalog, spans = _concat(workloads)
    if verbose:
        print(
            f"seed {seed}: {cohort.n_samples} samples, {cohort.n_variants} variants, "
            f"{len(layouts)} layouts, {len(catalog)} queries",
            flush=True,
        )
    cache = _measure_all(layouts, catalog, cluster, verbose)
    weight = _storage_weight(layouts, cohort, cluster)
    fit_queries = workloads["fit"]
    fit_weights = _fit_weights(fit_queries, config)

    def measurements(index: int, workload: str) -> list:
        start, end = spans[workload]
        return cache[index][start:end]

    def summarize(index: int, workload: str) -> dict:
        rows = measurements(index, workload)
        weights = fit_weights if workload == "fit" else None
        return _summarize(layouts[index], rows, weight, weights)

    def best(indices: list[int], workload: str) -> int:
        weights = fit_weights if workload == "fit" else None
        winner = indices[0]
        best_cost = float("inf")
        for index in indices:
            cost = measured_objective(layouts[index], measurements(index, workload), weight, weights)
            if cost < best_cost:
                best_cost = cost
                winner = index
        return winner

    by_family: dict[str, list[int]] = defaultdict(list)
    for index, layout in enumerate(layouts):
        by_family[layout.stats.family].append(index)
    pedigree_ids = by_family["family"] + by_family["pedigree"] + by_family["heuristic"]
    # Variant layouts stay in the pedigree-free pool so that selector can still
    # fall back to a full scan layout. They are also scored on their own.
    no_pedigree_ids = by_family["variant"] + by_family["generic"]

    methods: dict[str, dict] = {}
    methods["variant_oracle"] = _oracle("variant_oracle", best(by_family["variant"], "fit"), len(by_family["variant"]), summarize)
    methods["family_oracle"] = _oracle("family_oracle", best(by_family["family"], "fit"), len(by_family["family"]), summarize)
    methods["generic_oracle"] = _oracle("generic_oracle", best(by_family["generic"], "fit"), len(by_family["generic"]), summarize)
    methods["pedigree_oracle"] = _oracle("pedigree_oracle", best(pedigree_ids, "fit"), len(pedigree_ids), summarize)
    methods["full_oracle"] = _oracle("full_oracle", best(list(range(len(layouts))), "fit"), len(layouts), summarize)
    methods["full_oracle"]["notes"] = "exhaustive over the whole pool, tuned on the fit mix"
    methods["heuristic_fixed"] = _oracle("heuristic_fixed", heuristic_index, 1, summarize)
    methods["heuristic_fixed"]["oracle"] = False
    methods["heuristic_fixed"]["notes"] = "fixed nuclear families, fixed width and pack"

    budget = int(config.get("budget", 36))
    fraction = float(config.get("pilot_fraction", 0.5))
    hidden = int(config.get("hidden", 32))
    epochs = int(config.get("epochs", 80))

    def cached(workload: str):
        def _measure(index: int):
            return measurements(index, workload)

        return _measure

    analytical = select_analytical(
        layouts, fit_queries, cluster, weight, budget, cached("fit"), query_weights=fit_weights
    )
    methods["analytical_budget"] = _from_selection(analytical, summarize, oracle=False)

    random_sel = select_random(
        layouts,
        fit_queries,
        cluster,
        weight,
        budget,
        np.random.default_rng(seed + 3),
        cached("fit"),
        query_weights=fit_weights,
    )
    methods["random_budget"] = _from_selection(random_sel, summarize, oracle=False)

    n_pilot = pilot_count(len(layouts), budget, fraction)
    pilot = stratified_pilot(layouts, n_pilot, np.random.default_rng(seed + 5))
    neural = select_neural(
        layouts,
        fit_queries,
        cohort,
        cluster,
        embeddings,
        weight,
        budget,
        fraction,
        np.random.default_rng(seed + 5),
        cached("fit"),
        hide_pedigree=False,
        hidden=hidden,
        epochs=epochs,
        seed=seed,
        pilot_ids=pilot,
        query_weights=fit_weights,
    )
    methods["pedpack"] = _from_selection(neural, summarize, oracle=False)

    blinded = select_neural(
        layouts,
        fit_queries,
        cohort,
        cluster,
        embeddings,
        weight,
        budget,
        fraction,
        np.random.default_rng(seed + 5),
        cached("fit"),
        hide_pedigree=True,
        hidden=hidden,
        epochs=epochs,
        seed=seed + 1,
        pilot_ids=pilot,
        query_weights=fit_weights,
    )
    methods["pedpack_no_pedigree_features"] = _from_selection(blinded, summarize, oracle=False)

    sub_layouts = [layouts[i] for i in no_pedigree_ids]

    def measure_sub(local: int):
        return measurements(no_pedigree_ids[local], "fit")

    generic_sel = select_neural(
        sub_layouts,
        fit_queries,
        cohort,
        cluster,
        embeddings,
        weight,
        budget,
        fraction,
        np.random.default_rng(seed + 9),
        measure_sub,
        hide_pedigree=False,
        hidden=hidden,
        epochs=epochs,
        seed=seed + 2,
        query_weights=fit_weights,
    )
    generic_sel.index = no_pedigree_ids[generic_sel.index]
    methods["pedpack_generic_pool"] = _from_selection(generic_sel, summarize, oracle=False)
    methods["pedpack_generic_pool"]["notes"] = "same selector, sample groups ignore the pedigree"

    scan_queries = workloads["scan"]
    scan_sel = select_neural(
        layouts,
        scan_queries,
        cohort,
        cluster,
        embeddings,
        weight,
        budget,
        fraction,
        np.random.default_rng(seed + 11),
        cached("scan"),
        hidden=hidden,
        epochs=epochs,
        seed=seed + 3,
    )
    methods["pedpack_tuned_on_scan"] = _from_selection(scan_sel, summarize, oracle=False)
    methods["pedpack_tuned_on_scan"]["notes"] = "selector objective was the scan workload"

    selective_queries = workloads["selective"]
    selective_sel = select_neural(
        layouts,
        selective_queries,
        cohort,
        cluster,
        embeddings,
        weight,
        budget,
        fraction,
        np.random.default_rng(seed + 13),
        cached("selective"),
        hidden=hidden,
        epochs=epochs,
        seed=seed + 4,
    )
    methods["pedpack_tuned_on_selective"] = _from_selection(selective_sel, summarize, oracle=False)
    methods["pedpack_tuned_on_selective"]["notes"] = "selector objective was the selective workload"

    model = _model_report(
        neural,
        blinded,
        layouts,
        fit_queries,
        cohort,
        cluster,
        embeddings,
        weight,
        pilot,
        measurements,
        fit_weights,
    )
    sensitivity = _sensitivity(methods, layouts, fit_queries, cluster)
    analysis = {
        "best_on_heldout": _oracle(
            "best_on_heldout",
            best(list(range(len(layouts))), "heldout"),
            len(layouts),
            summarize,
        ),
        "best_on_scan": _oracle(
            "best_on_scan",
            best(list(range(len(layouts))), "scan"),
            len(layouts),
            summarize,
        ),
        "best_on_window": _oracle(
            "best_on_window",
            best(list(range(len(layouts))), "window"),
            len(layouts),
            summarize,
        ),
    }
    for bucket in analysis.values():
        bucket["analysis_only"] = True
        bucket["notes"] = "upper bound that was allowed to see this workload; not a method"

    return {
        "seed": seed,
        "n_samples": cohort.n_samples,
        "n_variants": cohort.n_variants,
        "n_trios": len(trios(cohort.father, cohort.mother)),
        "n_layouts": len(layouts),
        "n_fit_queries": len(fit_queries),
        "storage_seconds_per_byte": weight,
        "mendelian_errors": errors,
        "methods": methods,
        "model": model,
        "sensitivity": sensitivity,
        "analysis": analysis,
        "elapsed_seconds": time.perf_counter() - started,
    }


def _layouts(cohort: Cohort, embeddings: np.ndarray, cluster: Cluster, config: dict) -> tuple[list[PreparedLayout], int]:
    block = cluster.block_size
    packs = []
    if config.get("include_naive_pack", True):
        packs.append("naive")
    for count in config.get("pack_counts", [16]):
        packs.append(f"count:{int(count)}")
    fill = int(config.get("block_fill_blocks", 4))
    packs.append(f"bytes:{fill * block}")
    widths = [int(w) for w in config["widths"]]
    widths.append(cohort.n_variants)
    candidates, heuristic_index = build_pool(
        cohort,
        embeddings,
        widths=widths,
        generic_sizes=[int(s) for s in config["generic_sizes"]],
        agg_sizes=[int(s) for s in config["agg_sizes"]],
        min_kinship=[float(x) for x in config["min_kinship"]],
        lams=[float(x) for x in config["lams"]],
        orders=list(config.get("orders", ["group_major", "locus_major"])),
        packs=packs,
        heuristic_width=int(config.get("heuristic_width", 400)),
        heuristic_pack=str(config.get("heuristic_pack", "count:16")),
        heuristic_order=str(config.get("heuristic_order", "group_major")),
        include_three_gen=bool(config.get("include_three_gen", True)),
    )
    prepared = [prepare_layout(candidate, cohort, cluster) for candidate in candidates]
    for layout in prepared:
        if layout.stats.n_files < 1 or layout.stats.primary_bytes <= 0:
            raise RuntimeError(f"empty layout {layout.stats.label}")
    return prepared, heuristic_index


def _concat(workloads: dict[str, list[Query]]) -> tuple[list[Query], dict[str, tuple[int, int]]]:
    catalog: list[Query] = []
    spans: dict[str, tuple[int, int]] = {}
    for name in WORKLOADS:
        start = len(catalog)
        catalog.extend(workloads[name])
        spans[name] = (start, len(catalog))
        if start == len(catalog):
            raise RuntimeError(f"workload {name} is empty")
    return catalog, spans


def _measure_all(layouts, catalog, cluster, verbose: bool) -> list[list]:
    cache = []
    for index, layout in enumerate(layouts):
        cache.append(measure_many(layout, catalog, cluster))
        if verbose and (index + 1) % 40 == 0:
            print(f"  measured {index + 1}/{len(layouts)}", flush=True)
    uncovered = sum(1 for rows in cache if any(item.uncovered for item in rows))
    if uncovered:
        raise RuntimeError(f"{uncovered} layouts do not cover every query sample")
    return cache


def _storage_weight(layouts: list[PreparedLayout], cohort: Cohort, cluster: Cluster) -> float:
    full = [
        layout
        for layout in layouts
        if layout.stats.family == "variant" and layout.stats.width == cohort.n_variants
    ]
    if not full:
        raise RuntimeError("the pool has no full-width variant layout to price storage against")
    reference = min(full, key=lambda layout: (layout.stats.n_files, layout.stats.primary_bytes))
    query = Query(
        name="reference-scan",
        workload="scan",
        samples=tuple(range(cohort.n_samples)),
        v0=0,
        v1=cohort.n_variants,
        trios=tuple(trios(cohort.father, cohort.mother)),
    )
    seconds = measure(reference, query, cluster).seconds
    return storage_weight(seconds, reference.stats.storage_bytes, cluster.storage_scan_equivalents)


def _summarize(
    layout: PreparedLayout,
    rows: list,
    weight: float,
    query_weights: np.ndarray | None = None,
) -> dict:
    seconds = np.asarray([row.seconds for row in rows], dtype=np.float64)
    read = np.asarray([row.bytes_read for row in rows], dtype=np.float64)
    remote = np.asarray([row.remote_bytes for row in rows], dtype=np.float64)
    shuffle = np.asarray([row.shuffle_bytes for row in rows], dtype=np.float64)
    files = np.asarray([row.files_touched for row in rows], dtype=np.float64)
    total_read = float(read.sum())
    return {
        "objective": measured_objective(layout, rows, weight, query_weights),
        "mean_seconds": float(seconds.mean()),
        "seconds_std": float(seconds.std()),
        "mean_bytes": float(read.mean()),
        "remote_fraction": float(remote.sum() / total_read) if total_read else 0.0,
        "mean_shuffle": float(shuffle.mean()),
        "mean_files_touched": float(files.mean()),
        "files": int(layout.stats.n_files),
        "blocks": int(layout.stats.n_blocks),
        "storage_bytes": float(layout.stats.storage_bytes),
        "primary_bytes": float(layout.stats.primary_bytes),
        "build_seconds": float(layout.stats.build_seconds),
        "duplication": float(layout.stats.duplication),
        "trio_colocation": float(layout.stats.trio_colocation),
        "mean_within_kinship": float(layout.stats.mean_within_kinship),
        "mean_group_size": float(layout.stats.mean_group_size),
        "n_groups": int(layout.stats.n_groups),
        "width": int(layout.stats.width),
        "pack": layout.stats.pack,
        "order": layout.stats.order,
        "family": layout.stats.family,
        "label": layout.stats.label,
    }


def _choice(summary: dict) -> dict:
    keys = (
        "label",
        "family",
        "width",
        "pack",
        "order",
        "n_groups",
        "files",
        "blocks",
        "trio_colocation",
        "duplication",
        "mean_group_size",
        "storage_bytes",
        "build_seconds",
    )
    return {key: summary[key] for key in keys}


def _oracle(name: str, index: int, measured: int, summarize) -> dict:
    workloads = {workload: summarize(index, workload) for workload in WORKLOADS}
    return {
        "name": name,
        "index": int(index),
        "oracle": True,
        "layouts_measured": int(measured),
        "notes": "exhaustive on this family, tuned on the fit mix",
        "choice": _choice(workloads["fit"]),
        "workloads": workloads,
    }


def _from_selection(selection, summarize, oracle: bool) -> dict:
    workloads = {workload: summarize(selection.index, workload) for workload in WORKLOADS}
    return {
        "name": selection.notes,
        "index": int(selection.index),
        "oracle": oracle,
        "layouts_measured": int(selection.measured),
        "pilot": int(selection.pilot),
        "shortlist": int(selection.shortlist),
        "notes": selection.notes,
        "choice": _choice(workloads["fit"]),
        "workloads": workloads,
    }


def _model_report(
    neural,
    blinded,
    layouts,
    fit_queries,
    cohort,
    cluster,
    embeddings,
    weight,
    pilot,
    measurements,
    query_weights,
) -> dict:
    pilot_set = set(pilot)
    complement = [i for i in range(len(layouts)) if i not in pilot_set]
    true = np.asarray(
        [measured_objective(layouts[i], measurements(i, "fit"), weight, query_weights) for i in range(len(layouts))],
        dtype=np.float64,
    )
    analytical = np.asarray(
        [analytical_objective(layouts[i], fit_queries, cluster, weight, query_weights) for i in range(len(layouts))],
        dtype=np.float64,
    )
    neural_pred = neural.predicted_objective
    report = {
        "pilot": len(pilot),
        "shortlist": neural.shortlist,
        "layouts_measured": neural.measured,
        "val_mse": None if neural.surrogate is None else neural.surrogate.history.get("val_mse"),
        "blinded_choice_same": bool(neural.index == blinded.index),
        "spearman_neural": _spearman(neural_pred[complement], true[complement]) if complement else float("nan"),
        "spearman_analytical": _spearman(analytical[complement], true[complement]) if complement else float("nan"),
    }
    if complement and neural.surrogate is not None:
        predicted = predict_query_seconds(
            neural.surrogate,
            [layouts[i] for i in complement],
            fit_queries,
            cohort,
            cluster,
            embeddings,
            hide_pedigree=False,
        )
        actual = np.asarray(
            [[row.seconds for row in measurements(i, "fit")] for i in complement],
            dtype=np.float64,
        )
        analytical_seconds = np.asarray(
            [
                [analytical_objective_query(layouts[i], query, cluster) for query in fit_queries]
                for i in complement
            ],
            dtype=np.float64,
        )
        report["seconds_mape_neural"] = _mape(predicted, actual)
        report["seconds_mape_analytical"] = _mape(analytical_seconds, actual)
    else:
        report["seconds_mape_neural"] = float("nan")
        report["seconds_mape_analytical"] = float("nan")
    return report


def _fit_weights(queries: list[Query], config: dict) -> np.ndarray:
    """Per-query multipliers so the fit objective is a stated mixture.

    A raw sum lets one full scan outvote many family-scoped queries, because a
    scan reads the whole matrix. ``horizon_queries`` is how many future queries
    the build is amortized over. ``fit_mixture`` splits that horizon between
    scan-like queries and the rest of the fit list.
    """

    horizon = float(config.get("horizon_queries", 32))
    mixture = config.get("fit_mixture", {"scan": 0.25, "selective": 0.75})
    scan_share = float(mixture.get("scan", 0.25))
    selective_share = float(mixture.get("selective", 0.75))
    scan_at = [i for i, query in enumerate(queries) if query.name.startswith("fit-scan")]
    other_at = [i for i, query in enumerate(queries) if i not in set(scan_at)]
    weights = np.zeros(len(queries), dtype=np.float64)
    if scan_at:
        weights[scan_at] = horizon * scan_share / len(scan_at)
    if other_at:
        weights[other_at] = horizon * selective_share / len(other_at)
    if not scan_at and other_at:
        weights[other_at] = horizon / len(other_at)
    return weights


def analytical_objective_query(layout, query, cluster) -> float:
    from pedpack.cost import analytical_query

    return float(analytical_query(layout, query, cluster).seconds)


def _sensitivity(methods, layouts, fit_queries, cluster) -> dict:
    tracked = (
        "pedpack",
        "generic_oracle",
        "variant_oracle",
        "family_oracle",
        "heuristic_fixed",
        "analytical_budget",
        "full_oracle",
    )
    indexes = {name: methods[name]["index"] for name in tracked}
    failure = cluster.with_updates(failed_nodes=(0,))
    node_settings = (4, 8, 16) if cluster.n_nodes >= 8 else (max(2, cluster.n_nodes // 2), cluster.n_nodes)
    block_settings = (max(4096, cluster.block_size // 4), cluster.block_size, cluster.block_size * 2)
    rep_settings = tuple(dict.fromkeys((1, cluster.replication)))

    def run(name: str, alt: Cluster) -> dict:
        rows = measure_many(layouts[indexes[name]], fit_queries, alt)
        read = sum(row.bytes_read for row in rows)
        remote = sum(row.remote_bytes for row in rows)
        return {
            "mean_seconds": float(np.mean([row.seconds for row in rows])),
            "remote_fraction": float(remote / read) if read else 0.0,
            "unavailable_blocks": int(sum(row.unavailable_blocks for row in rows)),
            "mean_bytes": float(np.mean([row.bytes_read for row in rows])),
        }

    out: dict[str, dict] = {"failure_node0": {}, "nodes": {}, "block_size": {}, "replication": {}}
    for name in tracked:
        out["failure_node0"][name] = run(name, failure)
        out["nodes"][name] = {str(n): run(name, cluster.with_updates(n_nodes=n, failed_nodes=())) for n in node_settings}
        out["block_size"][name] = {
            str(size): run(name, cluster.with_updates(block_size=size, failed_nodes=())) for size in block_settings
        }
        out["replication"][name] = {
            str(rep): run(name, cluster.with_updates(replication=rep, failed_nodes=())) for rep in rep_settings
        }
    return out


def _spearman(left, right) -> float:
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    if left.size < 3 or float(np.std(left)) < 1e-12 or float(np.std(right)) < 1e-12:
        return float("nan")
    return float(np.corrcoef(_ranks(left), _ranks(right))[0, 1])


def _ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=np.float64)
    ranks[order] = np.arange(values.size, dtype=np.float64)
    return ranks


def _mape(pred, true) -> float:
    pred = np.asarray(pred, dtype=np.float64)
    true = np.asarray(true, dtype=np.float64)
    return float(np.mean(np.abs(pred - true) / np.maximum(np.abs(true), 1e-9)))

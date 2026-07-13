#!/usr/bin/env python3
import argparse
import json
import math
from pathlib import Path

import numpy as np

from bgo_features_v2 import load_bgo2


BGO_FEATURES = [
    "maxE", "sumE", "nHit", "nCl4", "nCl8", "cl1Sum4", "cl2Sum4",
    "cl1Size4", "cl1MaxFrac4", "cl1RmsDeg4", "allRmsDeg", "cl2OverCl1",
    "isolatedEFrac4", "seed15NCl4", "seed15Cl1Sum4", "seed15Cl1Size4",
    "e2OverE1", "nHit2T", "sum2T", "local3x3Frac", "meanHitE",
    "thetaSpan", "phiSpan",
]

HODO_FEATURES = [
    "thTotalE_MeV", "thMaxE_MeV", "thNHit",
    "tlcExpectedTotal", "tlcExpectedMax", "tlcNExpectedSeg",
    "matchedExpectedPhotons", "matchedExpectedFraction",
    "matchedTlcMaxExpected", "cherenkovDt_ns", "cherenkovDtValid",
]
TH_FEATURES = ["thTotalE_MeV", "thMaxE_MeV", "thNHit"]
FINAL_TH_FEATURES = [
    "thMatchedDEdx_MeV_per_mm", "absDeltaZ_mm", "absDeltaZValid",
    "absDeltaZLt90",
]
TLC_EXPECTED_FEATURES = [
    "tlcExpectedTotal", "tlcExpectedMax", "tlcNExpectedSeg",
]

# Raw segment IDs and absolute times are exported for diagnostics/reproduction,
# but excluded from the default isotropic classifier to avoid phi artifacts.
DIAGNOSTIC_HODO_FEATURES = [
    "thLeadingSeg", "thLeadingTime_ns", "tlcLeadingSeg",
    "tlcLeadingTime_ns", "matchedTlcLeadingSeg",
]

LOG_FEATURES = {
    "maxE", "sumE", "cl1Sum4", "cl2Sum4", "seed15Cl1Sum4", "sum2T",
    "meanHitE", "thTotalE_MeV", "thMaxE_MeV", "tlcExpectedTotal",
    "tlcExpectedMax", "matchedExpectedPhotons", "matchedTlcMaxExpected",
}
for _tag in ("eff01", "eff02", "eff05", "eff10"):
    LOG_FEATURES.update({
        f"tlcNpeTotal_{_tag}", f"tlcNpeMax_{_tag}",
        f"matchedNpeTotal_{_tag}", f"matchedNpeMax_{_tag}",
    })


def poisson_tlc_features(tag, npe_threshold):
    return [
        f"tlcNpeTotal_{tag}", f"tlcNpeMax_{tag}",
        f"tlcNpeNSegGe{npe_threshold}_{tag}",
    ]


def poisson_hodo_features(tag, npe_threshold):
    return TH_FEATURES + poisson_tlc_features(tag, npe_threshold) + [
        f"matchedNpeTotal_{tag}", f"matchedNpeFraction_{tag}",
        f"matchedNpeMax_{tag}",
        f"matchedNpeNSegGe{npe_threshold}_{tag}",
    ]


def fixed_tlc_binary_feature(tag, npe_threshold, matched=False):
    prefix = "matchedAnySegHit" if matched else "tlcAnySegHit"
    return f"{prefix}Ge{npe_threshold}_{tag}"


def ensure_fixed_binary_features(electron, pion, names):
    names = list(names)
    e_extra, p_extra = [], []
    for tag in ("eff01", "eff02", "eff05", "eff10"):
        for threshold in (1, 2, 3):
            for matched in (False, True):
                output = fixed_tlc_binary_feature(tag, threshold, matched)
                if output in names:
                    continue
                prefix = "matchedNpeNSeg" if matched else "tlcNpeNSeg"
                source = f"{prefix}Ge{threshold}_{tag}"
                if source not in names:
                    continue
                idx = names.index(source)
                e_extra.append((electron[:, idx] > 0).astype(np.float32))
                p_extra.append((pion[:, idx] > 0).astype(np.float32))
                names.append(output)
    if e_extra:
        electron = np.column_stack([electron] + e_extra)
        pion = np.column_stack([pion] + p_extra)
    return electron, pion, names


def fixed_response_metrics(electron, pion, e_hit, p_hit):
    metrics = {
        "e_keep": float(np.mean(e_hit)),
        "p_reject": float(np.mean(~p_hit)),
    }
    return add_intervals(metrics, len(electron), len(pion))


def fixed_tlc_response_table(electron, e_train, e_test,
                             pion, p_train, p_test, names):
    col = {name: i for i, name in enumerate(names)}
    table = []
    for tag in ("eff01", "eff02", "eff05", "eff10"):
        for npe_threshold in (1, 2, 3):
            nseg = f"tlcNpeNSegGe{npe_threshold}_{tag}"
            matched_nseg = f"matchedNpeNSegGe{npe_threshold}_{tag}"
            total = f"tlcNpeTotal_{tag}"
            if not all(name in col for name in (nseg, matched_nseg, total)):
                continue
            scenario = {
                "efficiency": float(tag[-2:]) / 100.0,
                "npe_threshold": npe_threshold,
                "primary_definition":
                    "any detector segment has Npe >= fixed threshold",
                "primary_any_segment": {
                    "all_100k": fixed_response_metrics(
                        electron, pion, electron[:, col[nseg]] > 0,
                        pion[:, col[nseg]] > 0),
                    "train": fixed_response_metrics(
                        e_train, p_train, e_train[:, col[nseg]] > 0,
                        p_train[:, col[nseg]] > 0),
                    "test": fixed_response_metrics(
                        e_test, p_test, e_test[:, col[nseg]] > 0,
                        p_test[:, col[nseg]] > 0),
                },
                "matched_pm1_ablation": {
                    "all_100k": fixed_response_metrics(
                        electron, pion,
                        electron[:, col[matched_nseg]] > 0,
                        pion[:, col[matched_nseg]] > 0),
                    "definition":
                        "any TLC segment in the BGO/TH azimuth +/-1 window has Npe >= fixed threshold",
                    "train": fixed_response_metrics(
                        e_train, p_train,
                        e_train[:, col[matched_nseg]] > 0,
                        p_train[:, col[matched_nseg]] > 0),
                    "test": fixed_response_metrics(
                        e_test, p_test, e_test[:, col[matched_nseg]] > 0,
                        p_test[:, col[matched_nseg]] > 0),
                },
                "event_total_npe_comparison": {
                    "all_100k": fixed_response_metrics(
                        electron, pion,
                        electron[:, col[total]] >= npe_threshold,
                        pion[:, col[total]] >= npe_threshold),
                    "definition":
                        "sum over all detector segments Npe >= fixed threshold",
                    "train": fixed_response_metrics(
                        e_train, p_train,
                        e_train[:, col[total]] >= npe_threshold,
                        p_train[:, col[total]] >= npe_threshold),
                    "test": fixed_response_metrics(
                        e_test, p_test,
                        e_test[:, col[total]] >= npe_threshold,
                        p_test[:, col[total]] >= npe_threshold),
                },
            }
            table.append(scenario)
    return table


def abs_delta_z_baseline(electron, e_train, e_test,
                         pion, p_train, p_test, names):
    required = {"absDeltaZValid", "absDeltaZLt90"}
    if not required.issubset(names):
        return None
    valid_idx = names.index("absDeltaZValid")
    cut_idx = names.index("absDeltaZLt90")

    def evaluate(electron, pion):
        e_valid = electron[:, valid_idx] > 0.5
        p_valid = pion[:, valid_idx] > 0.5
        e_accept = e_valid & (electron[:, cut_idx] > 0.5)
        p_accept = p_valid & (pion[:, cut_idx] > 0.5)
        out = add_intervals({
            "e_keep": float(np.mean(e_accept)),
            "p_reject": float(np.mean(~p_accept)),
        }, len(electron), len(pion))
        out["e_valid_fraction"] = float(np.mean(e_valid))
        out["p_valid_fraction"] = float(np.mean(p_valid))
        out["invalid_policy"] = "reject"
        if np.any(e_valid):
            out["e_keep_given_valid"] = float(np.mean(
                electron[e_valid, cut_idx] > 0.5))
        if np.any(p_valid):
            out["p_reject_given_valid"] = float(np.mean(
                pion[p_valid, cut_idx] <= 0.5))
        return out

    return {
        "definition": "accept abs(zReco-zPred) < 90 mm; invalid rejects",
        "all_100k": evaluate(electron, pion),
        "train": evaluate(e_train, p_train),
        "test": evaluate(e_test, p_test),
    }


def fixed_timing_response(electron, pion, names, pi0=None):
    points = [
        (0.0, "absDeltaZLt90"),
        (0.2, "absDeltaZLt90_sigT0p2ns"),
        (0.5, "absDeltaZLt90_sigT0p5ns"),
        (1.0, "absDeltaZLt90_sigT1p0ns"),
    ]
    if "absDeltaZValid" not in names:
        return None
    valid_idx = names.index("absDeltaZValid")
    result = []
    for sigma, cut_name in points:
        if cut_name not in names:
            continue
        cut_idx = names.index(cut_name)
        e_valid = electron[:, valid_idx] > 0.5
        p_valid = pion[:, valid_idx] > 0.5
        e_accept = e_valid & (electron[:, cut_idx] > 0.5)
        p_accept = p_valid & (pion[:, cut_idx] > 0.5)
        e_pass_given_valid = float(np.mean(
            electron[e_valid, cut_idx] > 0.5)) if np.any(e_valid) else None
        p_reject_given_valid = float(np.mean(
            pion[p_valid, cut_idx] <= 0.5)) if np.any(p_valid) else None
        e_accept_invalid_kept = (~e_valid) | (electron[:, cut_idx] > 0.5)
        p_reject_invalid_kept = p_valid & (pion[:, cut_idx] <= 0.5)
        item = {
            "sigma_t_per_end_ns": sigma,
            "cut": "abs(zReco-zPred) < 90 mm",
            "primary_result": "conditional_on_valid",
            "e_keep_given_valid": e_pass_given_valid,
            "p_reject_given_valid": p_reject_given_valid,
            "full_selection_policy": "require valid and pass cut",
            "e_keep": float(np.mean(e_accept)),
            "e_keep_95ci": wilson(float(np.mean(e_accept)), len(electron)),
            "p_reject": float(np.mean(~p_accept)),
            "p_reject_95ci": wilson(float(np.mean(~p_accept)), len(pion)),
            "invalid_kept_e_keep": float(np.mean(e_accept_invalid_kept)),
            "invalid_kept_p_reject": float(np.mean(p_reject_invalid_kept)),
            "e_valid_fraction": float(np.mean(e_valid)),
            "p_valid_fraction": float(np.mean(p_valid)),
        }
        if pi0 is not None:
            pi0_valid = pi0[:, valid_idx] > 0.5
            pi0_accept = pi0_valid & (pi0[:, cut_idx] > 0.5)
            pi0_keep = float(np.mean(pi0_accept))
            pi0_reject_given_valid = float(np.mean(
                pi0[pi0_valid, cut_idx] <= 0.5)) if np.any(pi0_valid) else None
            pi0_reject_invalid_kept = pi0_valid & (
                pi0[:, cut_idx] <= 0.5)
            item.update({
                "pi0_reject_given_valid": pi0_reject_given_valid,
                "pi0_keep": pi0_keep,
                "pi0_keep_95ci": wilson(pi0_keep, len(pi0)),
                "pi0_reject": 1.0 - pi0_keep,
                "pi0_reject_95ci": wilson(1.0 - pi0_keep, len(pi0)),
                "invalid_kept_pi0_reject": float(np.mean(
                    pi0_reject_invalid_kept)),
                "pi0_valid_fraction": float(np.mean(pi0_valid)),
            })
        result.append(item)
    return {
        "statistics": "all events (fixed detector operating point)",
        "response_model": (
            "independent Gaussian end-time smearing; deterministic event/segment hash; "
            "v_eff from runmeta; provisional because calibration is not supplied"),
        "points": result,
    }


def fixed_th_dedx_response(electron, pion, names, pi0=None):
    required = {"thGeomPath_mm", "thMatchedDEdx_MeV_per_mm"}
    if not required.issubset(names):
        return None
    path_idx = names.index("thGeomPath_mm")
    dedx_idx = names.index("thMatchedDEdx_MeV_per_mm")
    e_valid = electron[:, path_idx] > 0.0
    p_valid = pion[:, path_idx] > 0.0
    points = []
    # Fixed grid declared before inspecting strict timing data; no threshold is
    # selected here. Invalid geometry is kept, as a veto observable is absent.
    for threshold in (0.2, 0.3, 0.4, 0.5, 0.75, 1.0):
        e_accept = (~e_valid) | (electron[:, dedx_idx] < threshold)
        p_accept = (~p_valid) | (pion[:, dedx_idx] < threshold)
        item = {
            "threshold_MeV_per_mm": threshold,
            "rule": "accept dE/dx < threshold; invalid geometry keeps",
            "e_keep": float(np.mean(e_accept)),
            "e_keep_95ci": wilson(float(np.mean(e_accept)), len(electron)),
            "p_reject": float(np.mean(~p_accept)),
            "p_reject_95ci": wilson(float(np.mean(~p_accept)), len(pion)),
            "e_valid_fraction": float(np.mean(e_valid)),
            "p_valid_fraction": float(np.mean(p_valid)),
        }
        if pi0 is not None:
            pi0_valid = pi0[:, path_idx] > 0.0
            pi0_accept = (~pi0_valid) | (pi0[:, dedx_idx] < threshold)
            pi0_keep = float(np.mean(pi0_accept))
            item.update({
                "pi0_keep": pi0_keep,
                "pi0_keep_95ci": wilson(pi0_keep, len(pi0)),
                "pi0_reject": 1.0 - pi0_keep,
                "pi0_reject_95ci": wilson(1.0 - pi0_keep, len(pi0)),
                "pi0_valid_fraction": float(np.mean(pi0_valid)),
            })
        points.append(item)
    return {
        "statistics": "all events (fixed threshold grid, no optimization)",
        "truth_path_usage": "none; path=(thRMax-thRMin)/sin(theta_BGO)",
        "points": points,
    }


def fixed_combined_response(electron, pion, names, pi0=None):
    required = {
        "tlcAnySegHitGe1_eff05", "thGeomPath_mm",
        "thMatchedDEdx_MeV_per_mm", "absDeltaZValid", "absDeltaZLt90",
    }
    if not required.issubset(names):
        return None
    tlc_idx = names.index("tlcAnySegHitGe1_eff05")
    path_idx = names.index("thGeomPath_mm")
    dedx_idx = names.index("thMatchedDEdx_MeV_per_mm")
    dz_valid_idx = names.index("absDeltaZValid")
    dz_cut_idx = names.index("absDeltaZLt90")

    def masks(data):
        tlc = data[:, tlc_idx] > 0.5
        path_valid = data[:, path_idx] > 0.0
        dz_full = ((data[:, dz_valid_idx] > 0.5) &
                   (data[:, dz_cut_idx] > 0.5))
        dedx03 = (~path_valid) | (data[:, dedx_idx] < 0.3)
        dedx04 = (~path_valid) | (data[:, dedx_idx] < 0.4)
        return {
            "tlc_any_eff05_npe1": tlc,
            "tlc_any_and_dedx_lt0p3": tlc & dedx03,
            "tlc_any_and_dedx_lt0p4": tlc & dedx04,
            "tlc_any_and_full_dz": tlc & dz_full,
            "tlc_any_and_dedx_lt0p3_and_full_dz": tlc & dedx03 & dz_full,
            "tlc_any_and_dedx_lt0p4_and_full_dz": tlc & dedx04 & dz_full,
        }

    e_masks, p_masks = masks(electron), masks(pion)
    pi0_masks = masks(pi0) if pi0 is not None else None
    points = []
    for name in e_masks:
        e_keep = float(np.mean(e_masks[name]))
        p_reject = float(np.mean(~p_masks[name]))
        item = {
            "selection": name,
            "e_keep": e_keep,
            "e_keep_95ci": wilson(e_keep, len(electron)),
            "p_reject": p_reject,
            "p_reject_95ci": wilson(p_reject, len(pion)),
        }
        if pi0_masks is not None:
            pi0_reject = float(np.mean(~pi0_masks[name]))
            item.update({
                "pi0_reject": pi0_reject,
                "pi0_reject_95ci": wilson(pi0_reject, len(pi0)),
            })
        points.append(item)
    return {
        "statistics": "all events; all thresholds fixed before evaluation",
        "definitions": {
            "tlc": "any of 30 segments Npe>=1 at collection*PDE=5%",
            "dedx": "electron-like pass dE/dx below cut; invalid geometry keeps",
            "full_dz": "require valid and abs(zReco-zPred)<90 mm, sigma_t=0",
        },
        "primary_operating_point": "tlc_any_eff05_npe1",
        "points": points,
    }


def split_mask(data, names):
    event_id = data[:, names.index("eventID")].astype(np.int64)
    # Stable 50/50 event-level split, independent within every geometry/sample.
    return (event_id & 1) == 0


def transform(data, names):
    out = np.asarray(data, dtype=np.float64).copy()
    for name in LOG_FEATURES:
        if name in names:
            idx = names.index(name)
            out[:, idx] = np.log1p(np.maximum(0.0, out[:, idx]))
    return out


def eval_cut(e_score, p_score, threshold, orientation):
    if orientation > 0:
        e_reject = e_score >= threshold
        p_reject = p_score >= threshold
    else:
        e_reject = e_score <= threshold
        p_reject = p_score <= threshold
    return {
        "e_keep": float(1.0 - e_reject.mean()),
        "p_reject": float(p_reject.mean()),
    }


def fit_cut(e_score, p_score, target):
    values = np.unique(p_score)
    ps = np.sort(p_score)
    es = np.sort(e_score)
    best = None
    for orientation in (1, -1):
        if orientation > 0:
            p_reject_count = len(ps) - np.searchsorted(
                ps, values, side="left")
            e_keep_count = np.searchsorted(es, values, side="left")
        else:
            p_reject_count = np.searchsorted(ps, values, side="right")
            e_keep_count = len(es) - np.searchsorted(
                es, values, side="right")
        ok = p_reject_count / len(ps) >= target
        if not np.any(ok):
            continue
        keep = e_keep_count / len(es)
        candidates = np.flatnonzero(ok)
        idx = candidates[np.argmax(keep[candidates])]
        model = {
            "threshold": float(values[idx]),
            "orientation": orientation,
            "e_keep_train": float(keep[idx]),
            "p_reject_train": float(p_reject_count[idx] / len(ps)),
        }
        if best is None or model["e_keep_train"] > best["e_keep_train"]:
            best = model
    return best


def wilson(value, n):
    z = 1.959963984540054
    denominator = 1.0 + z * z / n
    center = (value + z * z / (2.0 * n)) / denominator
    half = z * math.sqrt(
        value * (1.0 - value) / n + z * z / (4.0 * n * n)
    ) / denominator
    return [float(center - half), float(center + half)]


def add_intervals(metrics, ne, npion):
    metrics = dict(metrics)
    metrics["e_keep_95ci"] = wilson(metrics["e_keep"], ne)
    metrics["p_reject_95ci"] = wilson(metrics["p_reject"], npion)
    return metrics


def active_indices(names, requested, e_train, p_train):
    indices = [names.index(name) for name in requested]
    pooled = np.vstack([e_train[:, indices], p_train[:, indices]])
    active = [
        idx for idx, std in zip(indices, pooled.std(axis=0))
        if std > 1e-10
    ]
    if not active:
        raise RuntimeError("all requested features are constant")
    return active


def fit_fisher(e_train, p_train, indices, reg=0.05):
    pooled = np.vstack([e_train[:, indices], p_train[:, indices]])
    mean = pooled.mean(axis=0)
    scale = pooled.std(axis=0)
    scale[scale < 1e-10] = 1.0
    ez = (e_train[:, indices] - mean) / scale
    pz = (p_train[:, indices] - mean) / scale
    ce = np.atleast_2d(np.cov(ez, rowvar=False))
    cp = np.atleast_2d(np.cov(pz, rowvar=False))
    covariance = 0.5 * (ce + cp)
    covariance = ((1.0 - reg) * covariance +
                  reg * np.eye(covariance.shape[0]))
    weights = np.linalg.solve(
        covariance, pz.mean(axis=0) - ez.mean(axis=0))
    return mean, scale, weights


def fisher_score(data, indices, model):
    mean, scale, weights = model
    return ((data[:, indices] - mean) / scale) @ weights


def scan_pair(e_train, p_train, indices, target):
    pooled = np.vstack([e_train[:, indices], p_train[:, indices]])
    mean = pooled.mean(axis=0)
    scale = pooled.std(axis=0)
    scale[scale < 1e-10] = 1.0
    ez = (e_train[:, indices] - mean) / scale
    pz = (p_train[:, indices] - mean) / scale
    signs = np.sign(pz.mean(axis=0) - ez.mean(axis=0))
    signs[signs == 0] = 1.0
    ez *= signs
    pz *= signs
    angles = np.linspace(0.0, np.pi / 2.0, 25)
    best = None
    for ia in range(len(indices)):
        for ib in range(ia + 1, len(indices)):
            for angle in angles:
                ca, cb = np.cos(angle), np.sin(angle)
                es = ca * ez[:, ia] + cb * ez[:, ib]
                ps = ca * pz[:, ia] + cb * pz[:, ib]
                kth = max(0, int(np.floor((1.0 - target) * len(ps))))
                threshold = float(np.partition(ps, kth)[kth])
                candidate = {
                    "threshold": threshold,
                    "orientation": 1,
                    "e_keep_train": float(np.mean(es < threshold)),
                    "p_reject_train": float(np.mean(ps >= threshold)),
                    "feature_a_index": int(indices[ia]),
                    "feature_b_index": int(indices[ib]),
                    "coefficient_a": float(ca * signs[ia] / scale[ia]),
                    "coefficient_b": float(cb * signs[ib] / scale[ib]),
                    "center_a": float(mean[ia]),
                    "center_b": float(mean[ib]),
                }
                if (best is None or
                        candidate["e_keep_train"] > best["e_keep_train"]):
                    best = candidate
    return best


def pair_score(data, model):
    return (
        model["coefficient_a"] *
        (data[:, model["feature_a_index"]] - model["center_a"]) +
        model["coefficient_b"] *
        (data[:, model["feature_b_index"]] - model["center_b"])
    )


def optimize_group(group_name, requested, names, e_train, e_test,
                   p_train, p_test, targets, methods):
    indices = active_indices(names, requested, e_train, p_train)
    transformed_e_train = transform(e_train, names)
    transformed_e_test = transform(e_test, names)
    transformed_p_train = transform(p_train, names)
    transformed_p_test = transform(p_test, names)
    fisher = None
    if "fisher" in methods:
        fisher = fit_fisher(
            transformed_e_train, transformed_p_train, indices)
        fisher_e_train = fisher_score(transformed_e_train, indices, fisher)
        fisher_p_train = fisher_score(transformed_p_train, indices, fisher)
    result = {
        "group": group_name,
        "candidate_features": [names[i] for i in indices],
        "log1p_features": sorted(set(names[i] for i in indices) & LOG_FEATURES),
        "targets": {},
    }
    for target in targets:
        target_result = {}
        if "univariate" in methods:
            univariate = []
            for idx in indices:
                fit = fit_cut(e_train[:, idx], p_train[:, idx], target)
                test = eval_cut(
                    e_test[:, idx], p_test[:, idx],
                    fit["threshold"], fit["orientation"])
                univariate.append({
                    "feature": names[idx],
                    **fit,
                    "test": add_intervals(test, len(e_test), len(p_test)),
                })
            univariate.sort(
                key=lambda item: item["e_keep_train"], reverse=True)
            target_result["best_univariate"] = univariate[:10]

        if "pair" in methods:
            pair = scan_pair(
                transformed_e_train, transformed_p_train, indices, target)
            pair["feature_a"] = names[pair["feature_a_index"]]
            pair["feature_b"] = names[pair["feature_b_index"]]
            pair_test = eval_cut(
                pair_score(transformed_e_test, pair),
                pair_score(transformed_p_test, pair),
                pair["threshold"], pair["orientation"])
            pair["test"] = add_intervals(
                pair_test, len(e_test), len(p_test))
            target_result["pair_linear"] = pair

        if "fisher" in methods:
            fisher_cut = fit_cut(fisher_e_train, fisher_p_train, target)
            fisher_test = eval_cut(
                fisher_score(transformed_e_test, indices, fisher),
                fisher_score(transformed_p_test, indices, fisher),
                fisher_cut["threshold"], fisher_cut["orientation"])
            mean, scale, weights = fisher
            target_result["fisher"] = {
                **fisher_cut,
                "features": [names[i] for i in indices],
                "centers": [float(v) for v in mean],
                "scales": [float(v) for v in scale],
                "weights": [float(v) for v in weights],
                "regularization": 0.05,
                "test": add_intervals(
                    fisher_test, len(e_test), len(p_test)),
            }
        result["targets"][f"{100.0 * target:.3f}"] = target_result
    return result


def comparable_meta(e_meta, p_meta):
    keys = [
        "version", "ncol", "nLayer", "nSector", "segmentationMode",
        "physicsFlag", "nSegTH", "nSegTLC", "thetaMin_deg",
        "thetaMax_deg", "threshold_MeV",
    ]
    mismatch = {}
    for key in keys:
        if not np.isclose(e_meta[key], p_meta[key], rtol=0, atol=1e-9):
            mismatch[key] = [e_meta[key], p_meta[key]]
    return mismatch


def repo_relative_source(path):
    if path is None:
        return None
    value = str(path)
    marker = "beta_org/"
    position = value.find(marker)
    return value[position:] if position >= 0 else value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--electron", type=Path, required=True)
    parser.add_argument("--pion", type=Path, required=True)
    parser.add_argument(
        "--pi0", type=Path,
        help="optional strict-timing pi0 BGO2 file for fixed delta-z efficiency")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--targets", default="0.93,0.94,0.95",
        help="comma-separated pion rejection targets")
    parser.add_argument(
        "--groups", default="bgo",
        help=("comma-separated groups; legacy_response_scan uses TH dE plus "
              "fixed TLC binary, final_response_scan additionally requires "
              "strict measured TH dE/dx and delta-z, while "
              "diagnostic_response_scan permits free continuous Npe cuts"))
    parser.add_argument(
        "--methods", default="univariate,pair,fisher",
        help="comma-separated: univariate,pair,fisher")
    args = parser.parse_args()

    targets = [float(value) for value in args.targets.split(",")]
    if any(value <= 0 or value >= 1 for value in targets):
        raise ValueError("targets must lie between zero and one")
    electron, e_meta, e_manifest = load_bgo2(args.electron, copy=True)
    pion, p_meta, p_manifest = load_bgo2(args.pion, copy=True)
    pi0 = None
    pi0_manifest = None
    mismatch = comparable_meta(e_meta, p_meta)
    if mismatch:
        raise ValueError(f"electron/pion metadata mismatch: {mismatch}")
    if e_meta["physicsFlag"] != 4:
        raise ValueError("production optimizer accepts physicsFlag=4 only")
    if e_manifest["features"] != p_manifest["features"]:
        raise ValueError("electron/pion feature order mismatch")
    e_readout = e_manifest.get("readoutMode", "legacy-readout")
    p_readout = p_manifest.get("readoutMode", "legacy-readout")
    if e_readout != p_readout:
        raise ValueError(
            f"electron/pion readoutMode mismatch: {e_readout}/{p_readout}")
    if e_readout not in ("legacy-readout", "strict-timing"):
        raise ValueError(f"unknown readoutMode: {e_readout}")
    if args.pi0 is not None:
        pi0, pi0_meta, pi0_manifest = load_bgo2(args.pi0, copy=True)
        pi0_mismatch = comparable_meta(e_meta, pi0_meta)
        if pi0_mismatch:
            raise ValueError(f"electron/pi0 metadata mismatch: {pi0_mismatch}")
        if pi0_manifest["features"] != e_manifest["features"]:
            raise ValueError("electron/pi0 feature order mismatch")
        pi0_readout = pi0_manifest.get("readoutMode", "legacy-readout")
        if pi0_readout != e_readout:
            raise ValueError(
                f"electron/pi0 readoutMode mismatch: {e_readout}/{pi0_readout}")
    names = list(e_manifest["features"])
    electron, pion, names = ensure_fixed_binary_features(
        electron, pion, names)
    if pi0 is not None:
        # Electron is passed only to apply the identical deterministic column
        # derivation/order to pi0 when reading an older BGO2 schema.
        _, pi0, pi0_names = ensure_fixed_binary_features(
            electron, pi0, list(e_manifest["features"]))
        if pi0_names != names:
            raise ValueError("derived pi0 feature order mismatch")
    e_train_mask = split_mask(electron, names)
    p_train_mask = split_mask(pion, names)
    e_train, e_test = electron[e_train_mask], electron[~e_train_mask]
    p_train, p_test = pion[p_train_mask], pion[~p_train_mask]
    if min(len(e_train), len(e_test), len(p_train), len(p_test)) < 20:
        raise ValueError("too few events after train/test split")

    group_features = {
        "bgo": BGO_FEATURES,
        "th": TH_FEATURES,
        "tlc_expected": TLC_EXPECTED_FEATURES,
        "th_tlc_expected": HODO_FEATURES,
        "hodo": HODO_FEATURES,
        "bgo_hodo": BGO_FEATURES + HODO_FEATURES,
        "all": BGO_FEATURES + HODO_FEATURES + DIAGNOSTIC_HODO_FEATURES,
    }
    if e_readout == "strict-timing":
        group_features["th_final"] = FINAL_TH_FEATURES
    diagnostic_response_groups = []
    legacy_response_groups = []
    final_response_groups = []
    for tag in ("eff01", "eff02", "eff05", "eff10"):
        for npe_threshold in (1, 2, 3):
            tlc_name = f"diagnostic_tlc_free_{tag}_npe{npe_threshold}"
            hodo_name = f"diagnostic_th_tlc_free_{tag}_npe{npe_threshold}"
            combined_name = (
                f"diagnostic_bgo_th_tlc_free_{tag}_npe{npe_threshold}")
            tlc_features = poisson_tlc_features(tag, npe_threshold)
            features = poisson_hodo_features(tag, npe_threshold)
            group_features[tlc_name] = tlc_features
            group_features[hodo_name] = features
            group_features[combined_name] = BGO_FEATURES + features
            diagnostic_response_groups.extend(
                [tlc_name, hodo_name, combined_name])

            binary = fixed_tlc_binary_feature(tag, npe_threshold)
            matched_binary = fixed_tlc_binary_feature(
                tag, npe_threshold, matched=True)
            final_name = f"final_{tag}_npe{npe_threshold}"
            matched_name = f"final_matched_{tag}_npe{npe_threshold}"
            legacy_name = f"legacy_bgo_th_tlc_{tag}_npe{npe_threshold}"
            legacy_matched_name = (
                f"legacy_bgo_th_tlc_matched_{tag}_npe{npe_threshold}")
            group_features[legacy_name] = BGO_FEATURES + TH_FEATURES + [binary]
            group_features[legacy_matched_name] = (
                BGO_FEATURES + TH_FEATURES + [matched_binary])
            legacy_response_groups.extend(
                [legacy_name, legacy_matched_name])
            if (e_readout == "strict-timing" and
                    all(name in names for name in
                        FINAL_TH_FEATURES + [binary])):
                group_features[final_name] = (
                    BGO_FEATURES + FINAL_TH_FEATURES + [binary])
                final_response_groups.append(final_name)
            if (e_readout == "strict-timing" and
                    all(name in names for name in
                        FINAL_TH_FEATURES + [matched_binary])):
                group_features[matched_name] = (
                    BGO_FEATURES + FINAL_TH_FEATURES + [matched_binary])
                final_response_groups.append(matched_name)
    requested_groups = []
    for group in args.groups.split(","):
        if group == "final_response_scan":
            if e_readout != "strict-timing":
                raise ValueError(
                    "final_response_scan requires strict-timing inputs")
            requested_groups.extend(final_response_groups)
        elif group == "legacy_response_scan":
            requested_groups.extend(legacy_response_groups)
        elif group == "diagnostic_response_scan":
            requested_groups.extend(diagnostic_response_groups)
        else:
            requested_groups.append(group)
    unknown = sorted(set(requested_groups) - set(group_features))
    if unknown:
        raise ValueError(f"unknown group(s): {unknown}")
    n_hit = names.index("nHit")
    ncl4 = names.index("nCl4")
    ncl8 = names.index("nCl8")
    methods = args.methods.split(",")
    unknown_methods = sorted(
        set(methods) - {"univariate", "pair", "fisher"})
    if unknown_methods:
        raise ValueError(f"unknown method(s): {unknown_methods}")
    result = {
        "status": "complete",
        "electron": str(args.electron),
        "pion": str(args.pion),
        "pi0": str(args.pi0) if args.pi0 is not None else None,
        "feature_files_ephemeral": {
            "electron": str(args.electron),
            "pion": str(args.pion),
            "pi0": str(args.pi0) if args.pi0 is not None else None,
            "note": "temporary extracted BGO2 files; regenerate from source_root_files",
        },
        "source_root_files": {
            "electron": repo_relative_source(e_manifest.get("input")),
            "pion": repo_relative_source(p_manifest.get("input")),
            "pi0": repo_relative_source(
                pi0_manifest.get("input") if pi0_manifest else None),
        },
        "extractor_provenance": {
            "classifier_inputs": e_manifest.get("classifierInputs"),
            "forbidden_inputs_read": e_manifest.get("forbiddenInputsRead"),
            "event_join": e_manifest.get("eventJoin"),
            "poisson_response": e_manifest.get("poissonResponse"),
            "timing_response": e_manifest.get("timingResponse"),
        },
        "metadata": {
            key: value for key, value in e_meta.items()
            if key not in ("magic",)
        },
        "readout_mode": e_readout,
        "split": "even eventID train, odd eventID test, independently per sample",
        "counts": {
            "electron_train": len(e_train), "electron_test": len(e_test),
            "pion_train": len(p_train), "pion_test": len(p_test),
        },
        "classifier_truth_usage": "sample label only",
        "forbidden_inputs": [
            "target", "pdg", "pid", "tlc.dE_truth_MeV",
            "th.chargedPath_truth_mm", "tlc.chargedPath_truth_mm",
        ],
        "baseline_any_bgo_hit": {
            "train": {
                "e_keep": float(np.mean(e_train[:, n_hit] == 0)),
                "p_reject": float(np.mean(p_train[:, n_hit] > 0)),
            },
            "test": add_intervals({
                "e_keep": float(np.mean(e_test[:, n_hit] == 0)),
                "p_reject": float(np.mean(p_test[:, n_hit] > 0)),
            }, len(e_test), len(p_test)),
        },
        "baseline_ncluster_eq_1": {
            "four_neighbor": {
                "train": {
                    "e_keep": float(np.mean(e_train[:, ncl4] == 1)),
                    "p_reject": float(np.mean(p_train[:, ncl4] != 1)),
                },
                "test": add_intervals({
                    "e_keep": float(np.mean(e_test[:, ncl4] == 1)),
                    "p_reject": float(np.mean(p_test[:, ncl4] != 1)),
                }, len(e_test), len(p_test)),
            },
            "eight_neighbor": {
                "train": {
                    "e_keep": float(np.mean(e_train[:, ncl8] == 1)),
                    "p_reject": float(np.mean(p_train[:, ncl8] != 1)),
                },
                "test": add_intervals({
                    "e_keep": float(np.mean(e_test[:, ncl8] == 1)),
                    "p_reject": float(np.mean(p_test[:, ncl8] != 1)),
                }, len(e_test), len(p_test)),
            },
        },
        "baseline_abs_delta_z_lt90_mm": (
            abs_delta_z_baseline(
                electron, e_train, e_test, pion, p_train, p_test, names)
            if e_readout == "strict-timing" else None),
        "fixed_timing_response": (
            fixed_timing_response(electron, pion, names, pi0)
            if e_readout == "strict-timing" else None),
        "fixed_th_dedx_response": (
            fixed_th_dedx_response(electron, pion, names, pi0)
            if e_readout == "strict-timing" else None),
        "fixed_combined_response": (
            fixed_combined_response(electron, pion, names, pi0)
            if e_readout == "strict-timing" else None),
        "fixed_tlc_binary_response": fixed_tlc_response_table(
            electron, e_train, e_test, pion, p_train, p_test, names),
        "continuous_tlc_policy": (
            "tlcExpected*, tlcNpeTotal, tlcNpeMax, and freely optimized Npe "
            "cuts are diagnostic/systematic only; final groups use one fixed "
            "binary any-segment hit per response point"),
        "groups": [],
    }
    for group in requested_groups:
        result["groups"].append(optimize_group(
            group, group_features[group], names,
            e_train, e_test, p_train, p_test, targets, methods))
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"wrote {args.output}")
    print(json.dumps({
        "metadata": result["metadata"],
        "counts": result["counts"],
        "baseline_any_bgo_hit": result["baseline_any_bgo_hit"],
        "groups_completed": [group["group"] for group in result["groups"]],
    }, indent=2))


if __name__ == "__main__":
    main()
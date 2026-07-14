#!/usr/bin/env python3
import argparse
import json
import struct
from pathlib import Path

import numpy as np


MAGIC = 0x42474F32
VERSION = 2
HEADER_FORMAT = "<IIIIiiiiiiddd"
HEADER_SIZE = struct.calcsize(HEADER_FORMAT)
HEADER_KEYS = (
    "magic", "version", "nrow", "ncol", "nLayer", "nSector",
    "segmentationMode", "physicsFlag", "nSegTH", "nSegTLC",
    "thetaMin_deg", "thetaMax_deg", "threshold_MeV",
)


def load_bgo2(path, copy=False):
    path = Path(path)
    with path.open("rb") as stream:
        values = struct.unpack(HEADER_FORMAT, stream.read(HEADER_SIZE))
    meta = dict(zip(HEADER_KEYS, values))
    if meta["magic"] != MAGIC or meta["version"] != VERSION:
        raise ValueError(f"{path}: unsupported magic/version")
    expected_size = HEADER_SIZE + 4 * meta["nrow"] * meta["ncol"]
    if path.stat().st_size != expected_size:
        raise ValueError(
            f"{path}: size={path.stat().st_size}, expected={expected_size}")
    manifest_path = Path(str(path) + ".json")
    manifest = json.loads(manifest_path.read_text())
    for key in ("nrow", "ncol", "nLayer", "nSector", "segmentationMode",
                "physicsFlag", "nSegTH", "nSegTLC"):
        if manifest[key] != meta[key]:
            raise ValueError(f"{path}: header/manifest mismatch for {key}")
    if len(manifest["features"]) != meta["ncol"]:
        raise ValueError(f"{path}: feature-name count mismatch")
    data = np.memmap(path, dtype="<f4", mode="r", offset=HEADER_SIZE,
                     shape=(meta["nrow"], meta["ncol"]))
    if copy:
        data = np.array(data)
    return data, meta, manifest


def validate(path, expect_rows=None, expect_physics=4):
    data, meta, manifest = load_bgo2(path)
    names = manifest["features"]
    col = {name: i for i, name in enumerate(names)}
    required = {
        "eventID", "nHit", "thPresent", "thTotalE_MeV", "thMaxE_MeV",
        "thNHit", "thLeadingSeg", "tlcPresent", "tlcExpectedTotal",
        "tlcExpectedMax", "tlcNExpectedSeg", "matchedExpectedPhotons",
        "matchedExpectedFraction", "cherenkovDt_ns", "cherenkovDtValid",
    }
    missing = sorted(required - set(names))
    if missing:
        raise ValueError(f"missing output features: {missing}")
    if expect_rows is not None and meta["nrow"] != expect_rows:
        raise ValueError(f"rows={meta['nrow']}, expected={expect_rows}")
    if expect_physics is not None and meta["physicsFlag"] != expect_physics:
        raise ValueError(
            f"physicsFlag={meta['physicsFlag']}, expected={expect_physics}")
    if meta["segmentationMode"] not in (0, 1, 2):
        raise ValueError("invalid segmentationMode")
    if not np.all(np.isfinite(data)):
        raise ValueError("non-finite output feature")
    event_ids = data[:, col["eventID"]]
    if not np.all(event_ids == np.rint(event_ids)):
        raise ValueError("non-integral eventID")
    if len(np.unique(event_ids)) != len(event_ids):
        raise ValueError("duplicate output eventID")
    if np.any((data[:, col["nHit"]] < 0) |
              (data[:, col["nHit"]] > meta["nLayer"] * meta["nSector"])):
        raise ValueError("BGO nHit out of range")
    if np.any((data[:, col["thNHit"]] < 0) |
              (data[:, col["thNHit"]] > meta["nSegTH"])):
        raise ValueError("TH nHit out of range")
    if np.any((data[:, col["tlcNExpectedSeg"]] < 0) |
              (data[:, col["tlcNExpectedSeg"]] > meta["nSegTLC"])):
        raise ValueError("TLC nSeg out of range")
    for name, nseg in (("thLeadingSeg", meta["nSegTH"]),
                       ("tlcLeadingSeg", meta["nSegTLC"]),
                       ("matchedTlcLeadingSeg", meta["nSegTLC"])):
        values = data[:, col[name]]
        if np.any((values < -1) | (values >= nseg)):
            raise ValueError(f"{name} out of range")
    fraction = data[:, col["matchedExpectedFraction"]]
    if np.any((fraction < -1e-6) | (fraction > 1.0 + 1e-5)):
        raise ValueError("matchedExpectedFraction out of range")
    if np.any(data[:, col["matchedExpectedPhotons"]] >
              data[:, col["tlcExpectedTotal"]] + 1e-3):
        raise ValueError("matched expected photons exceeds total")
    pc_names = ("pcSumE_MeV", "pcDownE_MeV", "pcUpE_MeV")
    if any(name in col for name in pc_names):
        if not all(name in col for name in pc_names):
            raise ValueError("incomplete photon-counter feature set")
        for name in pc_names:
            if np.any(data[:, col[name]] < 0):
                raise ValueError(f"negative photon-counter energy: {name}")
        if not np.allclose(
            data[:, col["pcSumE_MeV"]],
            data[:, col["pcDownE_MeV"]] + data[:, col["pcUpE_MeV"]],
            rtol=1e-6, atol=1e-6,
        ):
            raise ValueError("photon-counter sum differs from side energies")
    pc_gamma_names = (
        "pcGammaN", "pcGammaDownN", "pcGammaUpN",
        "pcGammaEnergy_MeV", "pcGammaDownEnergy_MeV",
        "pcGammaUpEnergy_MeV", "pcGammaMaxEnergy_MeV",
    )
    if any(name in col for name in pc_gamma_names):
        if not all(name in col for name in pc_gamma_names):
            raise ValueError("incomplete photon-counter gamma-entrance feature set")
        for name in pc_gamma_names:
            if np.any(data[:, col[name]] < 0):
                raise ValueError(f"negative photon-counter gamma value: {name}")
        for name in ("pcGammaN", "pcGammaDownN", "pcGammaUpN"):
            if np.any(data[:, col[name]] != np.rint(data[:, col[name]])):
                raise ValueError(f"non-integral photon-counter gamma count: {name}")
        if not np.allclose(
            data[:, col["pcGammaN"]],
            data[:, col["pcGammaDownN"]] + data[:, col["pcGammaUpN"]],
            rtol=0, atol=0,
        ):
            raise ValueError("photon-counter gamma count differs from side counts")
        if not np.allclose(
            data[:, col["pcGammaEnergy_MeV"]],
            data[:, col["pcGammaDownEnergy_MeV"]]
            + data[:, col["pcGammaUpEnergy_MeV"]],
            rtol=1e-6, atol=1e-6,
        ):
            raise ValueError("photon-counter gamma energy differs from side energies")
    valid = data[:, col["cherenkovDtValid"]]
    if np.any((valid != 0) & (valid != 1)):
        raise ValueError("cherenkovDtValid is not binary")
    if not np.all(data[:, col["thPresent"]] == 1):
        raise ValueError("TH event missing in strict output")
    if not np.all(data[:, col["tlcPresent"]] == 1):
        raise ValueError("TLC event missing in strict output")
    if any(manifest[key] != 0 for key in
           ("missingTH", "missingTLC", "extraTH", "extraTLC")):
        raise ValueError("strict output records missing/extra hodoscope events")
    if manifest["eventJoin"] != "eventID":
        raise ValueError("manifest does not declare eventID join")
    forbidden = {"target", "pdg", "pid", "dE_truth_MeV",
                 "chargedPath_truth_mm"}
    inputs = set(manifest["classifierInputs"])
    bad = sorted(inputs & forbidden)
    if bad or manifest["forbiddenInputsRead"]:
        raise ValueError(f"forbidden classifier input recorded: {bad}")

    if manifest.get("finalReadoutMetadata", False):
        final_required = {
            "bgoLeadingCell", "bgoLeadingTheta_deg", "bgoLeadingPhi_deg",
            "finalReadoutPresent", "thBgoMatchedSeg",
            "thBgoMatchedE_MeV", "thGeomPath_mm",
            "thMatchedDEdx_MeV_per_mm", "thMatchedZReco_mm",
            "thMatchedZRecoValid", "zPred_mm", "absDeltaZ_mm",
            "absDeltaZValid", "absDeltaZLt90",
            "absDeltaZ_sigT0p2ns_mm", "absDeltaZLt90_sigT0p2ns",
            "absDeltaZ_sigT0p5ns_mm", "absDeltaZLt90_sigT0p5ns",
            "absDeltaZ_sigT1p0ns_mm", "absDeltaZLt90_sigT1p0ns",
        }
        final_missing = sorted(final_required - set(names))
        if final_missing:
            raise ValueError(f"missing final-readout features: {final_missing}")
        if not np.all(data[:, col["finalReadoutPresent"]] == 1):
            raise ValueError("strict final output has missing readout rows")
        for name in ("thMatchedZRecoValid", "absDeltaZValid",
                     "absDeltaZLt90"):
            values = data[:, col[name]]
            if np.any((values != 0) & (values != 1)):
                raise ValueError(f"{name} is not binary")
        valid = data[:, col["absDeltaZValid"]] > 0.5
        if np.any(data[valid, col["absDeltaZ_mm"]] < 0):
            raise ValueError("valid absDeltaZ is negative")
        expected_lt90 = data[valid, col["absDeltaZ_mm"]] < 90.0
        if np.any((data[valid, col["absDeltaZLt90"]] > 0.5) !=
                  expected_lt90):
            raise ValueError("absDeltaZLt90 is inconsistent with absDeltaZ")
        if np.any(data[~valid, col["absDeltaZLt90"]] != 0):
            raise ValueError("invalid absDeltaZ passes the 90 mm cut")
        for abs_name, cut_name in (
                ("absDeltaZ_sigT0p2ns_mm", "absDeltaZLt90_sigT0p2ns"),
                ("absDeltaZ_sigT0p5ns_mm", "absDeltaZLt90_sigT0p5ns"),
                ("absDeltaZ_sigT1p0ns_mm", "absDeltaZLt90_sigT1p0ns")):
            cut_values = data[:, col[cut_name]]
            if np.any((cut_values != 0) & (cut_values != 1)):
                raise ValueError(f"{cut_name} is not binary")
            if np.any((cut_values[valid] > 0.5) !=
                      (data[valid, col[abs_name]] < 90.0)):
                raise ValueError(f"{cut_name} inconsistent with {abs_name}")
            if np.any(cut_values[~valid] != 0):
                raise ValueError(f"invalid delta-z passes {cut_name}")
        for name, limit in (("bgoLeadingCell",
                             meta["nLayer"] * meta["nSector"]),
                            ("thBgoMatchedSeg", meta["nSegTH"])):
            values = data[:, col[name]]
            if np.any((values < -1) | (values >= limit)):
                raise ValueError(f"{name} out of range")
        if manifest.get("timingResponse") is None:
            raise ValueError("strict timing output lacks response metadata")
        if manifest.get("thTimingSmearingAppliedInSource") != 0:
            raise ValueError("strict timing source was already smeared")
    else:
        if manifest.get("readoutMode", "legacy-readout") != "legacy-readout":
            raise ValueError("non-final manifest has inconsistent readoutMode")
        if "th.zReco_mm" in manifest["classifierInputs"]:
            raise ValueError("legacy manifest claims unavailable th.zReco_mm")
        if manifest.get("timingResponse") is not None:
            raise ValueError("legacy manifest has fake timing response")

    for tag in ("eff01", "eff02", "eff05", "eff10"):
        for threshold in (1, 2, 3):
            for prefix, count_prefix in (
                    ("tlcAnySegHit", "tlcNpeNSeg"),
                    ("matchedAnySegHit", "matchedNpeNSeg")):
                binary_name = f"{prefix}Ge{threshold}_{tag}"
                count_name = f"{count_prefix}Ge{threshold}_{tag}"
                if binary_name not in col:
                    continue
                values = data[:, col[binary_name]]
                if np.any((values != 0) & (values != 1)):
                    raise ValueError(f"{binary_name} is not binary")
                if np.any((values > 0.5) !=
                          (data[:, col[count_name]] > 0)):
                    raise ValueError(
                        f"{binary_name} inconsistent with {count_name}")
    summary_names = [
        "maxE", "sumE", "nHit", "thTotalE_MeV", "thMaxE_MeV", "thNHit",
        "tlcExpectedTotal", "tlcExpectedMax", "tlcNExpectedSeg",
        "matchedExpectedPhotons", "matchedExpectedFraction",
        "cherenkovDtValid",
    ]
    if manifest.get("finalReadoutMetadata", False):
        summary_names += [
            "thMatchedDEdx_MeV_per_mm", "absDeltaZ_mm",
            "absDeltaZValid", "absDeltaZLt90",
        ]
    if all(name in col for name in pc_names):
        summary_names += list(pc_names)
    if all(name in col for name in pc_gamma_names):
        summary_names += list(pc_gamma_names)
    summary = {
        name: {
            "min": float(np.min(data[:, col[name]])),
            "max": float(np.max(data[:, col[name]])),
            "mean": float(np.mean(data[:, col[name]])),
        }
        for name in summary_names
    }
    return meta, manifest, summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("feature_file", type=Path)
    parser.add_argument("--expect-rows", type=int)
    parser.add_argument("--expect-physics", type=int, default=4)
    args = parser.parse_args()
    meta, manifest, summary = validate(
        args.feature_file, args.expect_rows, args.expect_physics)
    print(json.dumps({
        "status": "ok",
        "file": str(args.feature_file),
        "metadata": {key: value for key, value in meta.items()
                     if key not in ("magic",)},
        "missingTH": manifest["missingTH"],
        "missingTLC": manifest["missingTLC"],
        "extraTH": manifest["extraTH"],
        "extraTLC": manifest["extraTLC"],
        "summary": summary,
    }, indent=2))


if __name__ == "__main__":
    main()

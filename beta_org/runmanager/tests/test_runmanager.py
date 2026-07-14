from __future__ import annotations

import importlib.util
import contextlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest import mock
from types import SimpleNamespace

import yaml


MODULE_PATH = Path(__file__).resolve().parents[1] / "runmanager.py"
SPEC = importlib.util.spec_from_file_location("beta_runmanager", MODULE_PATH)
assert SPEC and SPEC.loader
rm = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rm)


def base_manifest() -> dict:
    return {
        "schema": rm.EXPECTED_SCHEMA,
        "tag": "unit_scan",
        "build_dir": "build_opt",
        "events": 100000,
        "seed": 6302026,
        "threads": 8,
        "geometries": [
            {
                "name": "u15x15",
                "n_layer": 15,
                "n_sector": 15,
                "segmentation": "uniform_theta",
                "geometry_mode": "current",
                "photon_counter": "none",
            },
            {
                "name": "a10x20",
                "n_layer": 10,
                "n_sector": 20,
                "segmentation": "equal_solid_angle",
                "geometry_mode": "bgoegg_envelope",
                "photon_counter": "two_sided",
            },
        ],
        "primaries": ["e", "pim"],
    }


class RunManagerTests(unittest.TestCase):
    def write_manifest(self, directory: str, content: dict) -> Path:
        path = Path(directory) / "manifest.yml"
        path.write_text(yaml.safe_dump(content, sort_keys=False), encoding="utf-8")
        return path

    def test_cartesian_product_and_output_names(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = rm.load_manifest(self.write_manifest(directory, base_manifest()))
            jobs = rm.expand_jobs(manifest, Path(directory) / "state.json")
        self.assertEqual(len(jobs), 4)
        self.assertEqual(
            {job["key"] for job in jobs},
            {"u15x15:e", "u15x15:pim", "a10x20:e", "a10x20:pim"},
        )
        self.assertTrue(jobs[0]["output_root"].endswith("unit_scan_u15x15_e.root"))

    def test_rejects_different_schema(self):
        content = base_manifest()
        content["schema"] = "beta-bgo-v0"
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_manifest(directory, content)
            with self.assertRaisesRegex(rm.RunManagerError, "schema mismatch"):
                rm.load_manifest(path)

    def test_rejects_non_100k_events(self):
        content = base_manifest()
        content["events"] = 10
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_manifest(directory, content)
            with self.assertRaisesRegex(rm.RunManagerError, "events must be 100000"):
                rm.load_manifest(path)

    def test_rejects_duplicate_geometry_name(self):
        content = base_manifest()
        content["geometries"][1]["name"] = content["geometries"][0]["name"]
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_manifest(directory, content)
            with self.assertRaisesRegex(rm.RunManagerError, "duplicate geometry name"):
                rm.load_manifest(path)

    def test_rejects_unknown_geometry_and_counter_modes(self):
        for field, value in (
            ("geometry_mode", "exact_bgoegg"),
            ("photon_counter", "upstream_only"),
        ):
            content = base_manifest()
            content["geometries"][0][field] = value
            with tempfile.TemporaryDirectory() as directory:
                path = self.write_manifest(directory, content)
                with self.assertRaisesRegex(rm.RunManagerError, field):
                    rm.load_manifest(path)

    def test_current_geometry_rejects_photon_counter(self):
        content = base_manifest()
        content["geometries"][0]["photon_counter"] = "downstream"
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_manifest(directory, content)
            with self.assertRaisesRegex(rm.RunManagerError, "BGOegg geometry"):
                rm.load_manifest(path)

    def test_v3_offset_manifest_and_command(self):
        content = base_manifest()
        content["schema"] = rm.SCHEMA_V3
        content["geometries"][0]["bgo_z_offset_cm"] = 0
        content["geometries"][1]["bgo_z_offset_cm"] = -5
        with tempfile.TemporaryDirectory() as directory:
            manifest = rm.load_manifest(self.write_manifest(directory, content))
            jobs = rm.expand_jobs(manifest, Path(directory) / "state.json")
        egg_job = next(job for job in jobs if job["key"] == "a10x20:e")
        command = rm.bsub_command(manifest, egg_job)
        self.assertEqual(command[-1], "-5.0")
        self.assertEqual(egg_job["geometry"]["bgo_z_offset_cm"], -5.0)

    def test_v3_requires_offset_and_current_nonzero_is_rejected(self):
        content = base_manifest()
        content["schema"] = rm.SCHEMA_V3
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_manifest(directory, content)
            with self.assertRaisesRegex(rm.RunManagerError, "bgo_z_offset_cm"):
                rm.load_manifest(path)
        content["geometries"][0]["bgo_z_offset_cm"] = 1
        content["geometries"][1]["bgo_z_offset_cm"] = 0
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_manifest(directory, content)
            with self.assertRaisesRegex(rm.RunManagerError, "BGOegg geometry"):
                rm.load_manifest(path)

    def test_v3_published_bgoegg_frustum31(self):
        content = base_manifest()
        content["schema"] = rm.SCHEMA_V3
        content["geometries"] = [
            {
                "name": "published31",
                "n_layer": 31,
                "n_sector": 60,
                "segmentation": "bgoegg_published",
                "geometry_mode": "bgoegg_frustum",
                "photon_counter": "none",
                "bgo_z_offset_cm": -10,
            }
        ]
        with tempfile.TemporaryDirectory() as directory:
            manifest = rm.load_manifest(self.write_manifest(directory, content))
            jobs = rm.expand_jobs(manifest, Path(directory) / "state.json")
        geometry = jobs[0]["geometry"]
        self.assertAlmostEqual(geometry["theta_min_deg"], 5.336032242257286)
        self.assertEqual(geometry["theta_max_deg"], 168.0)
        command = rm.bsub_command(manifest, jobs[0])
        self.assertEqual(command[-5:], [
            "bgoegg_published", "e", "bgoegg_frustum", "none", "-10.0",
        ])

    def test_published_bgoegg_frustum_rejects_photon_counter(self):
        content = base_manifest()
        content["schema"] = rm.SCHEMA_V3
        content["geometries"] = [
            {
                "name": "published31",
                "n_layer": 31,
                "n_sector": 60,
                "segmentation": "bgoegg_published",
                "geometry_mode": "bgoegg_frustum",
                "photon_counter": "downstream",
                "bgo_z_offset_cm": 0,
            }
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_manifest(directory, content)
            with self.assertRaisesRegex(rm.RunManagerError, "overlap"):
                rm.load_manifest(path)

    def test_v4_coverage_manifest_and_command(self):
        content = base_manifest()
        content["schema"] = rm.SCHEMA_V4
        content["geometries"] = [content["geometries"][1]]
        content["geometries"][0].update(
            bgo_z_offset_cm=0, theta_min_deg=18, theta_max_deg=150,
        )
        with tempfile.TemporaryDirectory() as directory:
            manifest = rm.load_manifest(self.write_manifest(directory, content))
            jobs = rm.expand_jobs(manifest, Path(directory) / "state.json")
        command = rm.bsub_command(manifest, jobs[0])
        self.assertEqual(command[-3:], ["0.0", "18.0", "150.0"])
        self.assertEqual(jobs[0]["geometry"]["theta_min_deg"], 18.0)

    def test_v4_requires_ordered_egg_theta_range(self):
        content = base_manifest()
        content["schema"] = rm.SCHEMA_V4
        content["geometries"] = [content["geometries"][1]]
        content["geometries"][0].update(
            bgo_z_offset_cm=0, theta_min_deg=150, theta_max_deg=18,
        )
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_manifest(directory, content)
            with self.assertRaisesRegex(rm.RunManagerError, "less than"):
                rm.load_manifest(path)

    def test_schema_and_runmeta_validation(self):
        job = {
            "events": 100000,
            "seed": 6302026,
            "primary": "pim",
            "output_stem": "output/scan/unit__u15x15_pim",
            "geometry": {
                "n_layer": 15,
                "n_sector": 15,
                "segmentation": "uniform_theta",
                "geometry_mode": "current",
                "photon_counter": "none",
            },
        }
        inspection = {
            "file": {"zombie": "0", "recovered": "0"},
            "trees": {name: "0" for name in rm.EXPECTED_BRANCHES},
            "branches": {name: set(branches) for name, branches in rm.EXPECTED_BRANCHES.items()},
            "meta": {
                "nLayer": "15",
                "nSector": "15",
                "segmentationMode": "0",
                "physicsFlag": "4",
                "writeCalHit": "0",
                "neutronScale": "2",
                "inelasticBias": "3",
                "pionInelasticXSScale": "1.65",
                "seed": "6302026",
                "segmentation": "uniform_theta",
                "primary": "pim",
                "output": job["output_stem"],
                "nSegTH": "30",
                "nSegTLC": "30",
                "geometryMode": "0",
                "photonCounterMode": "0",
                "geometry": "current",
                "geometryModel": "spherical_shell_current",
                "photonCounter": "none",
                "pcNLayers": "8",
                "thetaMin_deg": "5.666",
                "thetaMax_deg": "170.302",
                "rMin_cm": "30",
                "thickness_cm": "20",
                "pcPbThickness_mm": "1",
                "pcScintiThickness_mm": "5",
                "pcZFront_cm": "52",
                "pcDownThetaInner_deg": "9.698",
                "pcDownThetaOuter_deg": "24",
                "pcUpThetaInner_deg": "5.666",
                "pcUpThetaOuter_deg": "36",
            },
            "vectors": {
                tree: {
                    branch: (225 if size == "n_cells" else size)
                    for branch, size in sizes.items()
                }
                for tree, sizes in rm.EXPECTED_VECTOR_SIZES.items()
            },
        }
        for name in ("evt", "calarr", "th", "tlc"):
            inspection["trees"][name] = "100000"
        inspection["trees"]["runmeta"] = "1"
        self.assertEqual(rm.validate_inspection(job, inspection), [])
        inspection["meta"]["physicsFlag"] = "3"
        errors = rm.validate_inspection(job, inspection)
        self.assertTrue(any("physicsFlag" in error for error in errors))

    def test_retry_predicate_is_documented_by_status_values(self):
        self.assertEqual(rm.RETRYABLE_VALIDATION, {"failed", "missing"})
        self.assertNotIn("valid", rm.RETRYABLE_VALIDATION)
        self.assertTrue(
            rm.is_retryable_job(
                {"lsf_state": "DONE", "validation": {"status": "missing"}}
            )
        )
        self.assertFalse(
            rm.is_retryable_job(
                {"lsf_state": "EXIT", "validation": {"status": "valid"}}
            )
        )

    def test_eligible_retry_jobs_selects_each_job_once(self):
        jobs = [
            {
                "key": "a10x20:e",
                "lsf_state": "DONE",
                "validation": {"status": "missing"},
            },
            {
                "key": "a10x20:pim",
                "lsf_state": "DONE",
                "validation": {"status": "valid"},
            },
        ]
        selected = rm.eligible_retry_jobs(jobs)
        self.assertEqual([job["key"] for job in selected], ["a10x20:e"])

    def test_existing_state_plain_dry_run_lists_all_jobs(self):
        with tempfile.TemporaryDirectory() as directory:
            build_dir = Path(directory)
            executable = build_dir / "beta"
            executable.touch(mode=0o755)
            manifest = {
                "tag": "unit_scan",
                "build_dir": str(build_dir),
                "geometries": [],
                "primaries": [],
            }
            jobs = [{"key": "g:e"}, {"key": "g:pim"}]
            args = SimpleNamespace(
                manifest="unused.yml",
                state_dir=directory,
                retry=False,
                dry_run=True,
            )
            output = io.StringIO()
            with (
                mock.patch.object(rm, "load_manifest", return_value=manifest),
                mock.patch.object(rm, "load_state", return_value={"jobs": []}),
                mock.patch.object(rm, "ensure_state_matches"),
                mock.patch.object(rm, "expand_jobs", return_value=jobs),
                mock.patch.object(rm, "bsub_command", side_effect=lambda _m, j: ["bsub", j["key"]]),
                contextlib.redirect_stdout(output),
            ):
                self.assertEqual(rm.submit(args), 0)
            self.assertEqual(output.getvalue().count("DRY-RUN"), 2)

    def test_binary_hash_drift_is_rejected(self):
        manifest = {
            "schema": rm.EXPECTED_SCHEMA,
            "tag": "unit_scan",
            "digest": "manifest-hash",
            "build_dir": "/unused",
        }
        state = {
            "schema": rm.EXPECTED_SCHEMA,
            "tag": "unit_scan",
            "manifest_sha256": "manifest-hash",
            "provenance": {
                "git_commit": "abc",
                "git_dirty": False,
                "build_executable": "/build/beta",
                "build_sha256": "old-hash",
            },
        }
        current = dict(state["provenance"], build_sha256="new-hash")
        with mock.patch.object(rm, "current_provenance", return_value=current):
            with self.assertRaisesRegex(rm.RunManagerError, "build_sha256"):
                rm.ensure_state_matches(manifest, state)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Small, schema-strict LSF run manager for beta_org BGO scans."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
from typing import Any

import yaml


RUNMANAGER_DIR = Path(__file__).resolve().parent
PROJECT_DIR = RUNMANAGER_DIR.parent
RUNNER = PROJECT_DIR / "scripts" / "run_bgo_sample.sh"
EXPECTED_SCHEMA = "beta-bgo-th-tlc-v1"
SUPPORTED_EVENTS = 100_000
ACTIVE_LSF_STATES = {"PEND", "RUN", "PROV", "WAIT", "SSUSP", "USUSP", "PSUSP"}
RETRYABLE_VALIDATION = {"failed", "missing"}
SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")

EXPECTED_BRANCHES = {
    "evt": {"eventID", "EdepCell_MeV", "EdepTarget_MeV"},
    "calhit": {
        "eventID", "copyNo", "t_ns", "dE_MeV", "pdg", "px_MeV",
        "py_MeV", "pz_MeV", "creator", "originType",
    },
    "target": {"eventID", "t_ns", "dE_MeV", "pdg", "px_MeV", "py_MeV", "pz_MeV"},
    "calarr": {"eventID", "dE_MeV", "pid"},
    "runmeta": {
        "nLayer", "nSector", "segmentationMode", "physicsFlag", "writeCalHit",
        "thetaMin_deg", "thetaMax_deg", "rMin_cm", "thickness_cm",
        "neutronScale", "inelasticBias", "pionInelasticXSScale", "seed",
        "segmentation", "primary", "output", "nSegTH", "nSegTLC",
        "tlcRefractiveIndex", "tlcLambdaMin_nm", "tlcLambdaMax_nm",
        "tlcResponseModel", "tlcCollectionEfficiencyApplied", "thBarZMin_mm",
        "thBarZMax_mm", "thEffectiveLightSpeed_mm_per_ns",
        "thTimingSmearingApplied", "thTimingModel", "thRMin_mm", "thRMax_mm",
    },
    "th": {
        "eventID", "dE_MeV", "time_ns", "timeLeft_ns", "timeRight_ns",
        "timeLeftMinusRight_ns", "zReco_mm", "chargedPath_truth_mm",
    },
    "tlc": {
        "eventID", "dE_truth_MeV", "cherenkovTime_ns", "chargedPath_truth_mm",
        "cherenkovPath_mm", "cherenkovExpectedPhotons",
    },
}

EXPECTED_VECTOR_SIZES = {
    "calarr": {"dE_MeV": "n_cells", "pid": "n_cells"},
    "th": {
        "dE_MeV": 30,
        "time_ns": 30,
        "timeLeft_ns": 30,
        "timeRight_ns": 30,
        "timeLeftMinusRight_ns": 30,
        "zReco_mm": 30,
        "chargedPath_truth_mm": 30,
    },
    "tlc": {
        "dE_truth_MeV": 30,
        "cherenkovTime_ns": 30,
        "chargedPath_truth_mm": 30,
        "cherenkovPath_mm": 30,
        "cherenkovExpectedPhotons": 30,
    },
}


class RunManagerError(RuntimeError):
    pass


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def require_int(value: Any, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise RunManagerError(f"{name} must be an integer")
    if not minimum <= value <= maximum:
        raise RunManagerError(f"{name} must be in [{minimum}, {maximum}]")
    return value


def require_name(value: Any, name: str) -> str:
    if not isinstance(value, str) or not SAFE_NAME.fullmatch(value):
        raise RunManagerError(
            f"{name} must match {SAFE_NAME.pattern!r}; got {value!r}"
        )
    return value


def manifest_digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_output(*arguments: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(PROJECT_DIR.parent), *arguments],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        raise RunManagerError(
            f"git {' '.join(arguments)} failed: {(result.stderr or result.stdout).strip()}"
        )
    return result.stdout.strip()


def current_provenance(manifest: dict[str, Any]) -> dict[str, Any]:
    executable = Path(manifest["build_dir"]) / "beta"
    if not executable.is_file():
        raise RunManagerError(f"beta executable is missing: {executable}")
    return {
        "git_commit": git_output("rev-parse", "HEAD"),
        "git_dirty": bool(git_output("status", "--porcelain")),
        "build_executable": str(executable.resolve()),
        "build_sha256": file_digest(executable),
    }


def load_manifest(path: Path) -> dict[str, Any]:
    path = path.resolve()
    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RunManagerError(f"cannot read manifest {path}: {exc}") from exc
    if not isinstance(raw, dict):
        raise RunManagerError("manifest root must be a mapping")

    allowed = {
        "schema", "tag", "build_dir", "events", "seed", "threads",
        "geometries", "primaries", "lsf",
    }
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise RunManagerError(f"unknown manifest fields: {', '.join(unknown)}")
    missing = sorted(allowed - {"lsf"} - set(raw))
    if missing:
        raise RunManagerError(f"missing manifest fields: {', '.join(missing)}")

    if raw["schema"] != EXPECTED_SCHEMA:
        raise RunManagerError(
            f"schema mismatch: expected {EXPECTED_SCHEMA!r}, got {raw['schema']!r}"
        )
    tag = require_name(raw["tag"], "tag")
    events = require_int(raw["events"], "events", 1, 2_000_000_000)
    if events != SUPPORTED_EVENTS:
        raise RunManagerError(
            f"events must be {SUPPORTED_EVENTS}; scripts/run_bgo_sample.sh uses the fixed 100k macro"
        )
    seed = require_int(raw["seed"], "seed", 1, 2_147_483_647)
    threads = require_int(raw["threads"], "threads", 1, 256)

    build_dir_raw = raw["build_dir"]
    if not isinstance(build_dir_raw, str) or not build_dir_raw:
        raise RunManagerError("build_dir must be a non-empty path string")
    build_dir = Path(build_dir_raw).expanduser()
    if not build_dir.is_absolute():
        build_dir = (PROJECT_DIR / build_dir).resolve()
    else:
        build_dir = build_dir.resolve()

    primary_raw = raw["primaries"]
    if not isinstance(primary_raw, list) or not primary_raw:
        raise RunManagerError("primaries must be a non-empty list")
    primaries: list[str] = []
    for item in primary_raw:
        if item not in {"e", "pim", "pi0"}:
            raise RunManagerError(f"unsupported primary: {item!r}")
        if item in primaries:
            raise RunManagerError(f"duplicate primary: {item}")
        primaries.append(item)

    geometry_raw = raw["geometries"]
    if not isinstance(geometry_raw, list) or not geometry_raw:
        raise RunManagerError("geometries must be a non-empty list")
    geometries: list[dict[str, Any]] = []
    geometry_names: set[str] = set()
    for index, item in enumerate(geometry_raw):
        if not isinstance(item, dict):
            raise RunManagerError(f"geometries[{index}] must be a mapping")
        unknown_geometry = sorted(
            set(item) - {"name", "n_layer", "n_sector", "segmentation"}
        )
        if unknown_geometry:
            raise RunManagerError(
                f"unknown fields in geometries[{index}]: {', '.join(unknown_geometry)}"
            )
        missing_geometry = sorted(
            {"name", "n_layer", "n_sector", "segmentation"} - set(item)
        )
        if missing_geometry:
            raise RunManagerError(
                f"missing fields in geometries[{index}]: {', '.join(missing_geometry)}"
            )
        name = require_name(item["name"], f"geometries[{index}].name")
        if name in geometry_names:
            raise RunManagerError(f"duplicate geometry name: {name}")
        geometry_names.add(name)
        n_layer = require_int(item["n_layer"], f"{name}.n_layer", 1, 200)
        n_sector = require_int(item["n_sector"], f"{name}.n_sector", 1, 360)
        if n_layer * n_sector > 20_000:
            raise RunManagerError(f"{name}: n_layer*n_sector exceeds 20000")
        segmentation = item["segmentation"]
        if segmentation not in {"uniform_theta", "equal_solid_angle"}:
            raise RunManagerError(
                f"{name}.segmentation must be uniform_theta or equal_solid_angle"
            )
        geometries.append(
            {
                "name": name,
                "n_layer": n_layer,
                "n_sector": n_sector,
                "segmentation": segmentation,
            }
        )

    lsf_raw = raw.get("lsf", {})
    if lsf_raw is None:
        lsf_raw = {}
    if not isinstance(lsf_raw, dict):
        raise RunManagerError("lsf must be a mapping")
    unknown_lsf = sorted(set(lsf_raw) - {"queue", "wall_minutes", "memory_mb"})
    if unknown_lsf:
        raise RunManagerError(f"unknown lsf fields: {', '.join(unknown_lsf)}")
    lsf: dict[str, Any] = {}
    if "queue" in lsf_raw:
        lsf["queue"] = require_name(lsf_raw["queue"], "lsf.queue")
    if "wall_minutes" in lsf_raw:
        lsf["wall_minutes"] = require_int(
            lsf_raw["wall_minutes"], "lsf.wall_minutes", 1, 7 * 24 * 60
        )
    if "memory_mb" in lsf_raw:
        lsf["memory_mb"] = require_int(
            lsf_raw["memory_mb"], "lsf.memory_mb", 128, 1_048_576
        )

    return {
        "path": str(path),
        "digest": manifest_digest(path),
        "schema": EXPECTED_SCHEMA,
        "tag": tag,
        "build_dir": str(build_dir),
        "events": events,
        "seed": seed,
        "threads": threads,
        "geometries": geometries,
        "primaries": primaries,
        "lsf": lsf,
    }


def state_path(manifest: dict[str, Any], override: str | None) -> Path:
    directory = Path(override).expanduser().resolve() if override else RUNMANAGER_DIR / ".state"
    return directory / f"{manifest['tag']}.json"


def load_state(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        state = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RunManagerError(f"cannot read state {path}: {exc}") from exc
    if not isinstance(state, dict) or not isinstance(state.get("jobs"), list):
        raise RunManagerError(f"invalid state format: {path}")
    return state


def write_state(path: Path, state: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    state["updated_at"] = now_iso()
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(state, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def expand_jobs(manifest: dict[str, Any], state_file: Path) -> list[dict[str, Any]]:
    build_dir = Path(manifest["build_dir"])
    logs = state_file.parent / "logs" / manifest["tag"]
    jobs: list[dict[str, Any]] = []
    for geometry in manifest["geometries"]:
        for primary in manifest["primaries"]:
            run_tag = f"{manifest['tag']}_{geometry['name']}"
            output_stem = f"output/scan/{run_tag}_{primary}"
            key = f"{geometry['name']}:{primary}"
            job_name = f"bgo_{manifest['tag']}_{geometry['name']}_{primary}"
            if len(job_name) > 120:
                digest = hashlib.sha1(job_name.encode()).hexdigest()[:10]
                job_name = f"{job_name[:109]}_{digest}"
            jobs.append(
                {
                    "key": key,
                    "job_name": job_name,
                    "run_tag": run_tag,
                    "geometry": geometry,
                    "primary": primary,
                    "events": manifest["events"],
                    "seed": manifest["seed"],
                    "threads": manifest["threads"],
                    "output_stem": output_stem,
                    "output_root": str(build_dir / f"{output_stem}.root"),
                    "stdout": str(logs / f"{geometry['name']}_{primary}.out"),
                    "stderr": str(logs / f"{geometry['name']}_{primary}.err"),
                    "lsf_id": None,
                    "lsf_state": "NOT_SUBMITTED",
                    "submission": "pending",
                    "validation": {"status": "not_run", "errors": []},
                }
            )
    return jobs


def make_state(manifest: dict[str, Any], jobs: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "schema": manifest["schema"],
        "tag": manifest["tag"],
        "manifest": manifest["path"],
        "manifest_sha256": manifest["digest"],
        "provenance": current_provenance(manifest),
        "created_at": now_iso(),
        "updated_at": now_iso(),
        "jobs": jobs,
    }


def ensure_state_identity(manifest: dict[str, Any], state: dict[str, Any]) -> None:
    if state.get("schema") != EXPECTED_SCHEMA:
        raise RunManagerError(
            f"state schema mismatch: expected {EXPECTED_SCHEMA!r}, got {state.get('schema')!r}"
        )
    if state.get("tag") != manifest["tag"]:
        raise RunManagerError("state tag does not match manifest")
    if state.get("manifest_sha256") != manifest["digest"]:
        raise RunManagerError(
            "manifest changed after state creation; use a new tag instead of mixing schemas/configurations"
        )


def ensure_state_matches(manifest: dict[str, Any], state: dict[str, Any]) -> None:
    ensure_state_identity(manifest, state)
    recorded = state.get("provenance")
    if not isinstance(recorded, dict):
        raise RunManagerError("state lacks build/source provenance; recreate or re-adopt it")
    current = current_provenance(manifest)
    for key in ("git_commit", "git_dirty", "build_executable", "build_sha256"):
        if recorded.get(key) != current.get(key):
            raise RunManagerError(
                f"provenance drift for {key}: state={recorded.get(key)!r}, current={current.get(key)!r}; "
                "use a new tag instead of mixing binaries/source states"
            )


def lsf_wall(minutes: int) -> str:
    hours, minute = divmod(minutes, 60)
    return f"{hours}:{minute:02d}"


def bsub_command(manifest: dict[str, Any], job: dict[str, Any]) -> list[str]:
    command = [
        "bsub",
        "-J", job["job_name"],
        "-n", str(job["threads"]),
        "-R", "span[hosts=1]",
        "-oo", job["stdout"],
        "-eo", job["stderr"],
    ]
    lsf = manifest["lsf"]
    if "queue" in lsf:
        command.extend(["-q", lsf["queue"]])
    if "wall_minutes" in lsf:
        command.extend(["-W", lsf_wall(lsf["wall_minutes"])])
    if "memory_mb" in lsf:
        command.extend(["-M", str(lsf["memory_mb"])])
        command.extend(["-R", f"rusage[mem={lsf['memory_mb']}]"])
    geometry = job["geometry"]
    command.extend(
        [
            "env",
            f"BETA_BUILD_DIR={manifest['build_dir']}",
            f"BETA_THREADS={job['threads']}",
            f"BETA_SEED={job['seed']}",
            str(RUNNER),
            job["run_tag"],
            str(geometry["n_layer"]),
            str(geometry["n_sector"]),
            geometry["segmentation"],
            job["primary"],
        ]
    )
    return command


def parse_bsub_id(output: str) -> int:
    match = re.search(r"Job <(\d+)>", output)
    if not match:
        raise RunManagerError(f"could not parse bsub job ID from: {output.strip()!r}")
    return int(match.group(1))


def query_lsf_state(job_id: int) -> str:
    result = subprocess.run(
        ["bjobs", "-a", "-noheader", "-o", "stat", str(job_id)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        return "UNKNOWN"
    tokens = result.stdout.strip().split()
    return tokens[0].upper() if tokens else "UNKNOWN"


def is_retryable_job(job: dict[str, Any]) -> bool:
    validation_status = job.get("validation", {}).get("status", "not_run")
    if validation_status == "valid":
        return False
    return (
        validation_status in RETRYABLE_VALIDATION
        or job.get("lsf_state") == "EXIT"
        or job.get("submission") == "failed"
    )


def eligible_retry_jobs(jobs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    seen: set[str] = set()
    for job in jobs:
        key = str(job.get("key", ""))
        if not key or key in seen:
            raise RunManagerError(f"duplicate or empty job key in state: {key!r}")
        seen.add(key)
        if job.get("lsf_id"):
            job["lsf_state"] = query_lsf_state(int(job["lsf_id"]))
        if is_retryable_job(job) and job.get("lsf_state") not in ACTIVE_LSF_STATES:
            selected.append(job)
    return selected


def submit(args: argparse.Namespace) -> int:
    manifest = load_manifest(Path(args.manifest))
    if not RUNNER.is_file() or not os.access(RUNNER, os.X_OK):
        raise RunManagerError(f"runner is missing or not executable: {RUNNER}")
    executable = Path(manifest["build_dir"]) / "beta"
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise RunManagerError(f"beta executable is missing: {executable}")

    state_file = state_path(manifest, args.state_dir)
    state = load_state(state_file)
    if state is None:
        if args.retry:
            raise RunManagerError("--retry requires existing state")
        jobs = expand_jobs(manifest, state_file)
        existing = [job["output_root"] for job in jobs if Path(job["output_root"]).exists()]
        if existing and not args.dry_run:
            raise RunManagerError(
                "refusing to overwrite outputs without state: " + ", ".join(existing[:4])
            )
        state = make_state(manifest, jobs)
        selected = state["jobs"]
    else:
        ensure_state_matches(manifest, state)
        if args.dry_run and not args.retry:
            # Inspection mode always shows the complete manifest expansion,
            # even when state already exists. It never mutates state.
            selected = expand_jobs(manifest, state_file)
        elif not args.retry:
            raise RunManagerError("state already exists; use --retry for failed/missing jobs only")
        else:
            selected = eligible_retry_jobs(state["jobs"])
            if not selected:
                print("no failed/missing jobs are eligible for retry")
                if not args.dry_run:
                    write_state(state_file, state)
                return 0

    if args.dry_run:
        for job in selected:
            print(f"DRY-RUN {job['key']}: {shlex.join(bsub_command(manifest, job))}")
        return 0

    state_file.parent.mkdir(parents=True, exist_ok=True)
    (state_file.parent / "logs" / manifest["tag"]).mkdir(parents=True, exist_ok=True)
    write_state(state_file, state)
    any_failure = False
    for job in selected:
        command = bsub_command(manifest, job)
        result = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if result.returncode == 0:
            job_id = parse_bsub_id(result.stdout + result.stderr)
            job["lsf_id"] = job_id
            job["lsf_state"] = "SUBMITTED"
            job["submission"] = "submitted"
            job["submitted_at"] = now_iso()
            job["validation"] = {"status": "not_run", "errors": []}
            print(f"submitted {job['key']}: LSF job {job_id}")
        else:
            any_failure = True
            job["submission"] = "failed"
            job["lsf_state"] = "SUBMIT_FAILED"
            job["submit_error"] = (result.stderr or result.stdout).strip()
            print(f"submit failed {job['key']}: {job['submit_error']}", file=sys.stderr)
        write_state(state_file, state)
    return 1 if any_failure else 0


def status(args: argparse.Namespace) -> int:
    manifest = load_manifest(Path(args.manifest))
    state_file = state_path(manifest, args.state_dir)
    state = load_state(state_file)
    if state is None:
        raise RunManagerError(f"state does not exist: {state_file}")
    ensure_state_matches(manifest, state)
    for job in state["jobs"]:
        if job.get("lsf_id"):
            job["lsf_state"] = query_lsf_state(int(job["lsf_id"]))
        validation = job.get("validation", {}).get("status", "not_run")
        output = "present" if Path(job["output_root"]).is_file() else "missing"
        print(
            f"{job['key']:<32} id={str(job.get('lsf_id') or '-'):>8} "
            f"lsf={job.get('lsf_state', 'UNKNOWN'):<13} output={output:<7} validate={validation}"
        )
    write_state(state_file, state)
    return 0


ROOT_DUMP_CODE = r'''
const char *path = gSystem->Getenv("BETA_VALIDATE_ROOT");
TFile *f = path ? TFile::Open(path, "READ") : nullptr;
if (!f || f->IsZombie()) {
  std::cout << "FILE\tzombie\t1" << std::endl;
  gSystem->Exit(20);
}
std::cout << "FILE\tzombie\t0" << std::endl;
std::cout << "FILE\trecovered\t" << (f->TestBit(TFile::kRecovered) ? 1 : 0) << std::endl;
const char *treeNames[] = {"evt", "calhit", "target", "calarr", "runmeta", "th", "tlc"};
for (const char *name : treeNames) {
  TTree *tree = dynamic_cast<TTree *>(f->Get(name));
  if (!tree) {
    std::cout << "TREE\t" << name << "\tMISSING" << std::endl;
    continue;
  }
  std::cout << "TREE\t" << name << "\t" << tree->GetEntries() << std::endl;
  TObjArray *branches = tree->GetListOfBranches();
  for (int i = 0; i < branches->GetEntries(); ++i) {
    TBranch *branch = static_cast<TBranch *>(branches->At(i));
    std::cout << "BRANCH\t" << name << "\t" << branch->GetName() << std::endl;
  }
}
TTree *meta = dynamic_cast<TTree *>(f->Get("runmeta"));
if (meta && meta->GetEntries() == 1) {
  int nLayer=-1,nSector=-1,segmentationMode=-1,physicsFlag=-1,writeCalHit=-1,nSegTH=-1,nSegTLC=-1;
  double neutronScale=-1,inelasticBias=-1,pionInelasticXSScale=-1,seed=-1;
  char segmentation[128]={0},primary[128]={0},output[1024]={0};
  if(meta->GetBranch("nLayer")) meta->SetBranchAddress("nLayer",&nLayer);
  if(meta->GetBranch("nSector")) meta->SetBranchAddress("nSector",&nSector);
  if(meta->GetBranch("segmentationMode")) meta->SetBranchAddress("segmentationMode",&segmentationMode);
  if(meta->GetBranch("physicsFlag")) meta->SetBranchAddress("physicsFlag",&physicsFlag);
  if(meta->GetBranch("writeCalHit")) meta->SetBranchAddress("writeCalHit",&writeCalHit);
  if(meta->GetBranch("neutronScale")) meta->SetBranchAddress("neutronScale",&neutronScale);
  if(meta->GetBranch("inelasticBias")) meta->SetBranchAddress("inelasticBias",&inelasticBias);
  if(meta->GetBranch("pionInelasticXSScale")) meta->SetBranchAddress("pionInelasticXSScale",&pionInelasticXSScale);
  if(meta->GetBranch("seed")) meta->SetBranchAddress("seed",&seed);
  if(meta->GetBranch("segmentation")) meta->SetBranchAddress("segmentation",&segmentation);
  if(meta->GetBranch("primary")) meta->SetBranchAddress("primary",&primary);
  if(meta->GetBranch("output")) meta->SetBranchAddress("output",&output);
  if(meta->GetBranch("nSegTH")) meta->SetBranchAddress("nSegTH",&nSegTH);
  if(meta->GetBranch("nSegTLC")) meta->SetBranchAddress("nSegTLC",&nSegTLC);
  meta->GetEntry(0);
  std::cout << std::setprecision(17);
  std::cout << "META\tnLayer\t" << nLayer << std::endl;
  std::cout << "META\tnSector\t" << nSector << std::endl;
  std::cout << "META\tsegmentationMode\t" << segmentationMode << std::endl;
  std::cout << "META\tphysicsFlag\t" << physicsFlag << std::endl;
  std::cout << "META\twriteCalHit\t" << writeCalHit << std::endl;
  std::cout << "META\tneutronScale\t" << neutronScale << std::endl;
  std::cout << "META\tinelasticBias\t" << inelasticBias << std::endl;
  std::cout << "META\tpionInelasticXSScale\t" << pionInelasticXSScale << std::endl;
  std::cout << "META\tseed\t" << seed << std::endl;
  std::cout << "META\tsegmentation\t" << segmentation << std::endl;
  std::cout << "META\tprimary\t" << primary << std::endl;
  std::cout << "META\toutput\t" << output << std::endl;
  std::cout << "META\tnSegTH\t" << nSegTH << std::endl;
  std::cout << "META\tnSegTLC\t" << nSegTLC << std::endl;
}
auto dumpVectorSizes = [](TTree *tree, const char *treeName, std::initializer_list<const char *> names) {
  if (!tree || tree->GetEntries() < 1) return;
  for (const char *name : names) {
    if (!tree->GetBranch(name)) continue;
    if (std::string(name) == "pid") {
      std::vector<int> *value=nullptr;
      tree->SetBranchAddress(name,&value);
      tree->GetEntry(0);
      std::cout << "VECTOR\t" << treeName << "\t" << name << "\t" << (value ? value->size() : 0) << std::endl;
      tree->ResetBranchAddresses();
    } else {
      std::vector<double> *value=nullptr;
      tree->SetBranchAddress(name,&value);
      tree->GetEntry(0);
      std::cout << "VECTOR\t" << treeName << "\t" << name << "\t" << (value ? value->size() : 0) << std::endl;
      tree->ResetBranchAddresses();
    }
  }
};
dumpVectorSizes(dynamic_cast<TTree *>(f->Get("calarr")),"calarr",{"dE_MeV","pid"});
dumpVectorSizes(dynamic_cast<TTree *>(f->Get("th")),"th",{"dE_MeV","time_ns","timeLeft_ns","timeRight_ns","timeLeftMinusRight_ns","zReco_mm","chargedPath_truth_mm"});
dumpVectorSizes(dynamic_cast<TTree *>(f->Get("tlc")),"tlc",{"dE_truth_MeV","cherenkovTime_ns","chargedPath_truth_mm","cherenkovPath_mm","cherenkovExpectedPhotons"});
f->Close();
gSystem->Exit(0);
'''


def inspect_root(path: Path, root_command: str = "root") -> dict[str, Any]:
    environment = os.environ.copy()
    environment["BETA_VALIDATE_ROOT"] = str(path)
    try:
        result = subprocess.run(
            [root_command, "-l", "-b", "-q", "-e", ROOT_DUMP_CODE],
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=120,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RunManagerError(f"ROOT inspection failed for {path}: {exc}") from exc
    if result.returncode != 0:
        detail = (result.stderr or result.stdout).strip()
        raise RunManagerError(f"ROOT inspection failed for {path}: {detail}")
    inspection: dict[str, Any] = {
        "file": {}, "trees": {}, "branches": {}, "meta": {}, "vectors": {}
    }
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) == 3 and fields[0] == "FILE":
            inspection["file"][fields[1]] = fields[2]
        elif len(fields) == 3 and fields[0] == "TREE":
            inspection["trees"][fields[1]] = fields[2]
        elif len(fields) == 3 and fields[0] == "BRANCH":
            inspection["branches"].setdefault(fields[1], set()).add(fields[2])
        elif len(fields) == 3 and fields[0] == "META":
            inspection["meta"][fields[1]] = fields[2]
        elif len(fields) == 4 and fields[0] == "VECTOR":
            inspection["vectors"].setdefault(fields[1], {})[fields[2]] = int(fields[3])
    return inspection


def float_equal(raw: str | None, expected: float, tolerance: float = 1e-12) -> bool:
    try:
        return abs(float(raw) - expected) <= tolerance * max(1.0, abs(expected))
    except (TypeError, ValueError):
        return False


def validate_inspection(job: dict[str, Any], inspection: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    if inspection["file"].get("zombie") != "0":
        errors.append("ROOT file is zombie")
    if inspection["file"].get("recovered") != "0":
        errors.append("ROOT file has the recovered bit")

    for tree, expected_branches in EXPECTED_BRANCHES.items():
        entry_raw = inspection["trees"].get(tree)
        if entry_raw in {None, "MISSING"}:
            errors.append(f"missing tree: {tree}")
            continue
        actual_branches = inspection["branches"].get(tree, set())
        if actual_branches != expected_branches:
            missing = sorted(expected_branches - actual_branches)
            extra = sorted(actual_branches - expected_branches)
            errors.append(f"schema mismatch in {tree}: missing={missing}, extra={extra}")

    for tree in ("evt", "calarr", "th", "tlc"):
        raw = inspection["trees"].get(tree)
        if raw not in {None, "MISSING"} and int(raw) != job["events"]:
            errors.append(f"{tree} entries={raw}, expected={job['events']}")
    runmeta_entries = inspection["trees"].get("runmeta")
    if runmeta_entries not in {None, "MISSING"} and int(runmeta_entries) != 1:
        errors.append(f"runmeta entries={runmeta_entries}, expected=1")

    meta = inspection["meta"]
    geometry = job["geometry"]
    expected_meta = {
        "nLayer": str(geometry["n_layer"]),
        "nSector": str(geometry["n_sector"]),
        "segmentationMode": "1" if geometry["segmentation"] == "equal_solid_angle" else "0",
        "physicsFlag": "4",
        "writeCalHit": "0",
        "segmentation": geometry["segmentation"],
        "primary": job["primary"],
        "output": job["output_stem"],
        "nSegTH": "30",
        "nSegTLC": "30",
    }
    for key, expected in expected_meta.items():
        if meta.get(key) != expected:
            errors.append(f"runmeta {key}={meta.get(key)!r}, expected={expected!r}")
    for key, expected in {
        "neutronScale": 2.0,
        "inelasticBias": 3.0,
        "pionInelasticXSScale": 1.65,
        "seed": float(job["seed"]),
    }.items():
        if not float_equal(meta.get(key), expected):
            errors.append(f"runmeta {key}={meta.get(key)!r}, expected={expected}")

    n_cells = geometry["n_layer"] * geometry["n_sector"]
    for tree, branch_sizes in EXPECTED_VECTOR_SIZES.items():
        for branch, expected_raw in branch_sizes.items():
            expected = n_cells if expected_raw == "n_cells" else int(expected_raw)
            actual = inspection["vectors"].get(tree, {}).get(branch)
            if actual != expected:
                errors.append(f"{tree}.{branch} vector size={actual}, expected={expected}")
    return errors


def validate(args: argparse.Namespace) -> int:
    manifest = load_manifest(Path(args.manifest))
    state_file = state_path(manifest, args.state_dir)
    state = load_state(state_file)
    needs_adoption = state is None or (
        args.adopt_existing and not isinstance(state.get("provenance"), dict)
    )
    if needs_adoption:
        if not args.adopt_existing:
            raise RunManagerError(
                f"state does not exist: {state_file}; use --adopt-existing for completed legacy jobs"
            )
        if state is not None:
            ensure_state_identity(manifest, state)
        jobs = expand_jobs(manifest, state_file)
        missing = [job["output_root"] for job in jobs if not Path(job["output_root"]).is_file()]
        if missing:
            raise RunManagerError(
                "--adopt-existing requires every declared output; missing: "
                + ", ".join(missing[:4])
            )
        for job in jobs:
            job["submission"] = "adopted"
            job["lsf_state"] = "UNKNOWN"
        state = make_state(manifest, jobs)
    ensure_state_matches(manifest, state)
    failures = 0
    for job in state["jobs"]:
        path = Path(job["output_root"])
        if not path.is_file():
            errors = [f"missing output: {path}"]
            status_name = "missing"
        else:
            try:
                inspection = inspect_root(path, args.root_command)
                errors = validate_inspection(job, inspection)
            except RunManagerError as exc:
                errors = [str(exc)]
            status_name = "failed" if errors else "valid"
        job["validation"] = {
            "status": status_name,
            "errors": errors,
            "checked_at": now_iso(),
        }
        if errors:
            failures += 1
            print(f"{status_name.upper()} {job['key']}")
            for error in errors:
                print(f"  - {error}")
        else:
            print(f"VALID {job['key']}: {path}")
    write_state(state_file, state)
    return 1 if failures else 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    submit_parser = subparsers.add_parser("submit", help="submit geometry x primary jobs")
    submit_parser.add_argument("manifest")
    submit_parser.add_argument("--dry-run", action="store_true")
    submit_parser.add_argument(
        "--retry", action="store_true", help="resubmit failed/missing jobs only"
    )
    submit_parser.add_argument("--state-dir", help=argparse.SUPPRESS)
    submit_parser.set_defaults(function=submit)

    status_parser = subparsers.add_parser("status", help="query bjobs without killing jobs")
    status_parser.add_argument("manifest")
    status_parser.add_argument("--state-dir", help=argparse.SUPPRESS)
    status_parser.set_defaults(function=status)

    validate_parser = subparsers.add_parser("validate", help="validate ROOT outputs and runmeta")
    validate_parser.add_argument("manifest")
    validate_parser.add_argument(
        "--adopt-existing",
        action="store_true",
        help="create state for a complete pre-runmanager output set",
    )
    validate_parser.add_argument("--root-command", default="root", help=argparse.SUPPRESS)
    validate_parser.add_argument("--state-dir", help=argparse.SUPPRESS)
    validate_parser.set_defaults(function=validate)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.function(args))
    except RunManagerError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

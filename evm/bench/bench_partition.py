#!/usr/bin/env python3
"""Benchmark harness for the evidence_modeler per-partition step.

Runs the Rust (and optionally Perl) evidence_modeler binary on a selected
partition directory, reports wall-clock time and peak RSS, and optionally
compares outputs.

Example:
    python3 bench/bench_partition.py \
        --partition example/Rhodotorula_sphaerocarpa/EVM/scaffold_1/scaffold_1_262149-1061113 \
        --weights example/Rhodotorula_sphaerocarpa/EVM/weights.txt \
        --iterations 3
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def find_repo_root() -> Path:
    p = Path(__file__).resolve().parent.parent.parent
    return p


def build_release(repo: Path) -> Path:
    print("Building release binary...")
    subprocess.run(["cargo", "build", "--release", "-p", "evm-utils"], cwd=repo / "evm", check=True)
    return repo / "evm" / "target" / "release" / "evidence_modeler"


def run_rust(binary: Path, partition: Path, weights: Path) -> dict:
    genome = (partition / "genome.softmasked.fa").resolve()
    gene_pred = (partition / "gene_predictions.gff3").resolve()
    protein = (partition / "protein_alignments.gff3").resolve()
    transcript = (partition / "transcript_alignments.gff3").resolve()
    weights = weights.resolve()

    cmd = [
        str(binary),
        "-G", str(genome),
        "-g", str(gene_pred),
        "-p", str(protein),
        "-w", str(weights),
        "--min_intron_length", "20",
        "--terminal_intergenic_re_search", "10000",
        "-o", "/dev/null",
    ]
    if transcript.exists():
        cmd.extend(["-e", str(transcript)])

    env = os.environ.copy()
    env["RUST_LOG"] = "warn"

    start = time.perf_counter()
    proc = subprocess.run(
        ["/usr/bin/time", "-v"] + cmd,
        capture_output=True,
        text=True,
        cwd=partition,
        env=env,
    )
    elapsed = time.perf_counter() - start

    if proc.returncode != 0:
        print("STDOUT:", proc.stdout[-2000:] if len(proc.stdout) > 2000 else proc.stdout)
        print("STDERR:", proc.stderr[-2000:] if len(proc.stderr) > 2000 else proc.stderr)
        raise RuntimeError(f"evidence_modeler failed with code {proc.returncode}")

    rss_kb = None
    for line in proc.stderr.splitlines():
        if "Maximum resident set size" in line:
            rss_kb = int(line.split(":")[-1].strip())
            break

    return {"elapsed": elapsed, "rss_kb": rss_kb}


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark evidence_modeler on a partition")
    parser.add_argument("--partition", required=True, type=Path, help="Path to partition directory")
    parser.add_argument("--weights", required=True, type=Path, help="Path to weights file")
    parser.add_argument("--iterations", type=int, default=3, help="Number of iterations")
    parser.add_argument("--skip-build", action="store_true", help="Skip cargo build --release")
    args = parser.parse_args()

    repo = find_repo_root()
    binary = repo / "evm" / "target" / "release" / "evidence_modeler"
    if not args.skip_build:
        binary = build_release(repo)

    print(f"Benchmarking {args.partition}")
    print(f"Weights: {args.weights}")
    print(f"Binary:  {binary}")
    print(f"Iterations: {args.iterations}\n")

    results = []
    for i in range(args.iterations):
        print(f"Run {i + 1}/{args.iterations} ...", flush=True)
        r = run_rust(binary, args.partition, args.weights)
        results.append(r)
        print(f"  wall={r['elapsed']:.3f}s  rss={r['rss_kb'] / 1024:.1f}MB")

    times = [r["elapsed"] for r in results]
    rss = [r["rss_kb"] for r in results]
    summary = {
        "partition": str(args.partition),
        "binary": str(binary),
        "iterations": args.iterations,
        "wall_seconds": {
            "min": min(times),
            "max": max(times),
            "mean": sum(times) / len(times),
        },
        "peak_rss_mb": {
            "min": min(rss) / 1024,
            "max": max(rss) / 1024,
            "mean": sum(rss) / len(rss) / 1024,
        },
    }

    print("\nSummary:")
    print(json.dumps(summary, indent=2))

    out = repo / "evm" / "bench" / "latest_bench.json"
    with open(out, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

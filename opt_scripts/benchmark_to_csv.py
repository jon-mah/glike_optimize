#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path


def parse_log(filepath):
    """Extract NUM_TREES, runtime, logp, and true from a log file."""

    with open(filepath, "r") as f:
        text = f.read()

    # NUM_TREES = 1
    num_trees_match = re.search(r"NUM_TREES\s*=\s*(\d+)", text)

    # Finished in 1h 14m 58.7s
    runtime_match = re.search(r"Finished in\s+(.+)", text)

    # logp = -11452.473967799393  (true = -11486.13318118242)
    logp_match = re.search(
        r"logp\s*=\s*([+-]?\d+(?:\.\d+)?)\s*\(true\s*=\s*([+-]?\d+(?:\.\d+)?)\)",
        text
    )

    if num_trees_match is None:
        raise ValueError(f"Could not find NUM_TREES in {filepath}")

    if runtime_match is None:
        raise ValueError(f"Could not find runtime in {filepath}")

    if logp_match is None:
        raise ValueError(f"Could not find logp/true in {filepath}")

    return {
        "file": filepath.name,
        "NUM_TREES": int(num_trees_match.group(1)),
        "runtime": runtime_match.group(1).strip(),
        "logp": float(logp_match.group(1)),
        "true": float(logp_match.group(2)),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Extract information from *_simulate_ARG.log files."
    )

    parser.add_argument(
        "--indir",
        required=True,
        help="Directory containing *_simulate_ARG.log files"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV file"
    )

    args = parser.parse_args()

    indir = Path(args.indir)

    files = sorted(indir.glob("*_simulate_ARG.log"))

    if not files:
        print(f"No *_simulate_ARG.log files found in {indir}")
        return

    results = []

    for filepath in files:
        print(f"Processing: {filepath}")
        try:
            results.append(parse_log(filepath))
        except ValueError as e:
            print(f"WARNING: {e}")

    # Write CSV
    fieldnames = [
        "file",
        "NUM_TREES",
        "runtime",
        "logp",
        "true",
    ]

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"\nWrote {len(results)} rows to {args.output}")


if __name__ == "__main__":
    main()

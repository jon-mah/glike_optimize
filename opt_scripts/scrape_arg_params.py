#!/usr/bin/env python3

import argparse
import ast
import csv
import re
from pathlib import Path


def parse_log(filepath):
    """Extract NUM_TREES and parameter estimates from a log file."""

    with open(filepath, "r") as f:
        text = f.read()

    # Extract NUM_TREES
    num_trees_match = re.search(r"NUM_TREES\s*=\s*(\d+)", text)

    if num_trees_match is None:
        raise ValueError(f"Could not find NUM_TREES in {filepath}")

    num_trees = int(num_trees_match.group(1))

    # Extract the Estimated dictionary
    estimated_match = re.search(
        r"Estimated:\s*(\{.*?\})",
        text
    )

    if estimated_match is None:
        raise ValueError(f"Could not find Estimated parameters in {filepath}")

    # Safely convert the string representation of the dictionary
    estimates = ast.literal_eval(estimated_match.group(1))

    return {
        "file": filepath.name,
        "NUM_TREES": num_trees,
        **estimates
    }


def main():
    parser = argparse.ArgumentParser(
        description="Extract parameter estimates from *_simulate_ARG.log files."
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

    if not results:
        print("No valid log files found.")
        return

    # Get all parameter names
    parameter_names = set()

    for result in results:
        parameter_names.update(
            key for key in result.keys()
            if key not in {"file", "NUM_TREES"}
        )

    # Sort parameters alphabetically
    parameter_names = sorted(parameter_names)

    # Make sure file and NUM_TREES are the first two columns
    fieldnames = [
        "file",
        "NUM_TREES",
        *parameter_names
    ]

    # Write CSV
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(results)

    print(f"\nWrote {len(results)} rows to {args.output}")


if __name__ == "__main__":
    main()
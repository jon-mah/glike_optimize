#!/usr/bin/env python3
"""
Parse CMA-ES optimization generations from a gLike log and write them to CSV.

Expected log records:
    Generation   1  Likelihood = -266715.702819
    {'t1': 14.7, 't2': 49.3, ...}

Usage:
    python parse_generations.py input.log
    python parse_generations.py input.log -o generations.csv
    python parse_generations.py input.log -o generations.csv --include-initial
"""

import argparse
import ast
import csv
import re
import sys
from pathlib import Path

GENERATION_RE = re.compile(
    r"^\s*Generation\s+(\d+)\s+Likelihood\s*=\s*([-+0-9.eE]+)\s*$"
)
PARAMETER_DICT_RE = re.compile(r"^\s*(\{.*\})\s*$")


def parse_log(log_path, include_initial=False):
    """Return parsed generation records and a stable list of parameter names."""
    records = []
    parameter_names = []
    known_parameters = set()
    pending_generation = None

    with open(log_path, "r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            generation_match = GENERATION_RE.match(line)
            if generation_match:
                pending_generation = {
                    "generation": int(generation_match.group(1)),
                    "log_likelihood": float(generation_match.group(2)),
                    "_line": line_number,
                }
                continue

            # A parameter dictionary normally appears immediately after a
            # generation's likelihood line.
            dict_match = PARAMETER_DICT_RE.match(line)
            if not dict_match:
                continue

            try:
                parameters = ast.literal_eval(dict_match.group(1))
            except (ValueError, SyntaxError):
                continue

            if not isinstance(parameters, dict):
                continue

            if pending_generation is not None:
                record = {
                    "generation": pending_generation["generation"],
                    "log_likelihood": pending_generation["log_likelihood"],
                }
                record.update(parameters)
                records.append(record)

                for name in parameters:
                    if name not in known_parameters:
                        known_parameters.add(name)
                        parameter_names.append(name)

                pending_generation = None
            elif include_initial and not records:
                # Capture the dictionary after "Initial parameters:" only.
                # Avoid accidentally treating unrelated dictionaries as initial
                # parameters by checking the preceding line.
                # The preceding-line check is handled in the main loop below
                # in typical logs; this fallback is intentionally not used.
                pass

    if include_initial:
        initial = extract_initial_parameters(log_path)
        if initial is not None:
            initial_record = {
                "generation": 0,
                "log_likelihood": "",
                **initial,
            }
            records.insert(0, initial_record)
            for name in initial:
                if name not in known_parameters:
                    known_parameters.add(name)
                    parameter_names.append(name)

    return records, parameter_names


def extract_initial_parameters(log_path):
    """Extract the dictionary immediately following 'Initial parameters:'."""
    expecting_dict = False
    with open(log_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip() == "Initial parameters:":
                expecting_dict = True
                continue
            if expecting_dict:
                match = PARAMETER_DICT_RE.match(line)
                if match:
                    try:
                        value = ast.literal_eval(match.group(1))
                        return value if isinstance(value, dict) else None
                    except (ValueError, SyntaxError):
                        return None
                if line.strip():
                    return None
    return None


def write_csv(records, parameter_names, output_path):
    """Write records with a consistent schema and blank cells for missing data."""
    columns = ["generation", "log_likelihood"] + parameter_names
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(
        description="Extract CMA-ES generation parameters from a gLike log."
    )
    parser.add_argument("input_log", help="Path to the optimization log")
    parser.add_argument(
        "-o", "--output",
        help="Output CSV path (default: <input_log_stem>_generations.csv)",
    )
    parser.add_argument(
        "--include-initial",
        action="store_true",
        help="Include initial parameters as generation 0 (likelihood left blank)",
    )
    args = parser.parse_args()

    input_path = Path(args.input_log)
    if not input_path.is_file():
        parser.error("Input log does not exist or is not a file: {}".format(input_path))

    output_path = Path(args.output) if args.output else input_path.with_name(
        input_path.stem + "_generations.csv"
    )

    records, parameter_names = parse_log(
        input_path, include_initial=args.include_initial
    )
    if not records:
        print(
            "No generation records found. Check that the log contains lines like "
            "'Generation 1 Likelihood = ...' followed by a parameter dictionary.",
            file=sys.stderr,
        )
        return 1

    write_csv(records, parameter_names, output_path)
    print("Wrote {} generation records to {}".format(len(records), output_path))
    print("Parameters: {}".format(", ".join(parameter_names)))
    return 0


if __name__ == "__main__":
    sys.exit(main())

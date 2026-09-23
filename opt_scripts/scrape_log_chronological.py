#!/usr/bin/env python3

import argparse
import ast
import csv
import re
from pathlib import Path


def get_replicate(log_file):
    """
    Extract the replicate number from the input filename.

    For example:

        gLike_3G09_12273717_1.log  -> 1
        gLike_3G09_12273717_10.log -> 10

    The replicate is defined as the final underscore-delimited
    component of the filename before '.log'.
    """

    stem = log_file.stem

    try:
        replicate = stem.rsplit("_", 1)[1]
    except IndexError:
        raise ValueError(
            f"Could not determine replicate from filename: {log_file}"
        )

    return replicate


def parse_iteration(line):
    """
    Attempt to extract an optimizer iteration number from a log line.

    Recognizes forms such as:

        Iteration 1
        iteration 1
        Iteration: 1
        iteration: 1
        Iteration = 1
        iteration = 1

    Returns the iteration number if found, otherwise None.
    """

    match = re.search(
        r"\biteration\s*(?:[:=]|\s)\s*(\d+)\b",
        line,
        flags=re.IGNORECASE
    )

    if match:
        return int(match.group(1))

    return None


def parse_log(log_file):
    """
    Parse a glike log file.

    Captures only lines whose first non-whitespace character is '{',
    followed by a Python dictionary and a log-likelihood value.

    Each captured record is additionally assigned:

        iteration
            The optimizer iteration associated with the record.

        replicate
            The replicate number extracted from the input filename.

    Example captured line:

        {'t1': 2.7675, 't2': 12.94295, 't3': 100} -75619.74213766128

    Examples ignored:

        x_up: {'t1': 3.65125, 't2': 12.94295, 't3': 100} -75986.61
        x: {'t1': 2.7675, 't2': 12.94295, 't3': 100} -75619.74
        x_down: {'t1': 1.88375, 't2': 12.94295, 't3': 100} -75207.96

    The returned records remain in exactly the order in which they
    appeared in the log file.
    """

    results = []

    replicate = get_replicate(log_file)

    # Current optimizer iteration.
    iteration = None

    with open(log_file, "r") as f:

        for line_number, raw_line in enumerate(f, start=1):

            line = raw_line.strip()

            # Ignore empty lines.
            if not line:
                continue

            # ---------------------------------------------------------
            # Check whether this line announces a new iteration.
            #
            # This is done BEFORE checking for parameter dictionaries,
            # so the iteration number is available for subsequent
            # parameter records.
            # ---------------------------------------------------------
            detected_iteration = parse_iteration(line)

            if detected_iteration is not None:
                iteration = detected_iteration
                continue

            # ---------------------------------------------------------
            # Only process lines whose first character is "{"
            #
            # This prevents x_up, x, x_down, etc. from being captured.
            # ---------------------------------------------------------
            if not line.startswith("{"):
                continue

            # ---------------------------------------------------------
            # Separate the dictionary from the final likelihood.
            #
            # Using rsplit() means that the last whitespace-separated
            # field is treated as the likelihood.
            # ---------------------------------------------------------
            try:
                parameter_string, likelihood_string = line.rsplit(
                    None, 1
                )
            except ValueError:
                print(
                    f"Warning: could not split line {line_number}:"
                )
                print(f"    {line}")
                continue

            # ---------------------------------------------------------
            # Parse the parameter dictionary.
            # ---------------------------------------------------------
            try:
                params = ast.literal_eval(parameter_string)
            except (ValueError, SyntaxError) as e:
                print(
                    f"Warning: could not parse dictionary "
                    f"on line {line_number}: {e}"
                )
                print(f"    {line}")
                continue

            # Make sure what we parsed really is a dictionary.
            if not isinstance(params, dict):
                print(
                    f"Warning: object on line {line_number} "
                    f"is not a dictionary."
                )
                continue

            # ---------------------------------------------------------
            # Parse the likelihood.
            # ---------------------------------------------------------
            try:
                log_likelihood = float(likelihood_string)
            except ValueError:
                print(
                    f"Warning: could not parse likelihood "
                    f"on line {line_number}: "
                    f"{likelihood_string}"
                )
                continue

            # ---------------------------------------------------------
            # Store the result.
            # ---------------------------------------------------------
            results.append(
                {
                    "params": params,
                    "log_likelihood": log_likelihood,
                    "iteration": iteration,
                    "replicate": replicate,
                }
            )

    return results


def write_csv(results, output_file):
    """
    Write parsed results to CSV while preserving log-file order.
    """

    if not results:
        raise RuntimeError(
            "No parameter/log-likelihood pairs were found."
        )

    # Use the parameter order from the first dictionary.
    parameter_names = list(results[0]["params"].keys())

    # Put iteration and replicate first, followed by the parameters
    # and finally the likelihood.
    fieldnames = (
        ["replicate", "iteration"]
        + parameter_names
        + ["log_likelihood"]
    )

    with open(output_file, "w", newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            extrasaction="ignore"
        )

        writer.writeheader()

        for result in results:

            row = result["params"].copy()

            row["replicate"] = result["replicate"]
            row["iteration"] = result["iteration"]
            row["log_likelihood"] = result["log_likelihood"]

            writer.writerow(row)


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Extract parameter dictionaries, optimizer iterations, "
            "replicate numbers, and log likelihoods from a glike log file."
        )
    )

    parser.add_argument(
        "log_file",
        type=Path,
        help="Path to the glike .log file."
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to the output CSV file."
    )

    args = parser.parse_args()

    if not args.log_file.is_file():
        raise FileNotFoundError(
            f"Log file not found: {args.log_file}"
        )

    print(f"Reading: {args.log_file}")

    replicate = get_replicate(args.log_file)

    print(f"Replicate: {replicate}")

    results = parse_log(args.log_file)

    print(
        f"Found {len(results):,} "
        "parameter/log-likelihood records."
    )

    # Warn if parameter records were found before an iteration marker.
    records_without_iteration = sum(
        result["iteration"] is None
        for result in results
    )

    if records_without_iteration:
        print(
            f"Warning: {records_without_iteration:,} records "
            "had no detected optimizer iteration."
        )

    write_csv(results, args.output)

    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()

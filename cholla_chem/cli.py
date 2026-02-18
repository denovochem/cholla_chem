import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from cholla_chem.main import resolve_compounds_to_smiles


def _infer_input_format(path: str) -> str:
    ext = Path(path).suffix.lower()
    if ext == ".csv":
        return "csv"
    if ext == ".tsv":
        return "tsv"
    if ext == ".smi":
        return "smi"
    return "txt"


def _infer_output_format(path: str) -> str:
    ext = Path(path).suffix.lower()
    if ext == ".json":
        return "json"
    if ext == ".csv":
        return "csv"
    if ext == ".tsv":
        return "tsv"
    raise ValueError(f"Unsupported output extension: {ext} (use .json, .csv, .tsv)")


def read_names_from_file(
    path: str,
    *,
    input_format: str | None = None,
    input_column: str = "name",
    encoding: str = "utf-8",
) -> List[str]:
    fmt = (input_format or _infer_input_format(path)).lower()

    if fmt == "txt":
        with open(path, "r", encoding=encoding, newline="") as f:
            return [line.strip() for line in f if line.strip()]

    out: List[str] = []

    if fmt in ("csv", "tsv"):
        delimiter = "," if fmt == "csv" else "\t"
        with open(path, "r", encoding=encoding, newline="") as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            if not reader.fieldnames:
                raise ValueError(
                    f"{path} has no header row; expected a column like '{input_column}'"
                )
            if input_column not in reader.fieldnames:
                raise ValueError(
                    f"{path} is missing column '{input_column}'. Columns: {reader.fieldnames}"
                )
            for row in reader:
                val = (row.get(input_column) or "").strip()
                if val:
                    out.append(val)
            return out

    if fmt == "smi":
        with open(path, "r", encoding=encoding, newline="") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    out.append(" ".join(parts[1:]).strip())
                else:
                    out.append(parts[0])
        return out

    raise ValueError(f"Unsupported input format: {fmt} (use txt, csv, tsv, smi)")


def write_results(
    results: Dict,
    *,
    output_path: str | None,
    output_format: str | None = None,
    encoding: str = "utf-8",
) -> None:
    if not output_path:
        out_text = "\n".join(f"{k}\t{v}" for k, v in results.items())
        print(out_text)
        return

    fmt = (output_format or _infer_output_format(output_path)).lower()

    if fmt == "json":
        with open(output_path, "w", encoding=encoding) as f:
            json.dump(results, f, ensure_ascii=False, indent=2)
            f.write("\n")
        return

    if fmt in ("csv", "tsv"):
        delimiter = "," if fmt == "csv" else "\t"
        with open(output_path, "w", encoding=encoding, newline="") as f:
            w = csv.DictWriter(f, fieldnames=["name", "smiles"], delimiter=delimiter)
            w.writeheader()
            for name, smiles in results.items():
                w.writerow({"name": name, "smiles": smiles})
        return

    raise ValueError(f"Unsupported output format: {fmt}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="cholla-chem")
    p.add_argument("names", nargs="*", type=str, help="Chemical names to resolve")
    p.add_argument("--input", "-i", type=str, help="Text file with one name per line")
    p.add_argument(
        "--input-format",
        type=str,
        default=None,
        choices=["txt", "csv", "tsv", "smi"],
        help="Optional override for input format (otherwise inferred from file extension)",
    )
    p.add_argument(
        "--input-column",
        type=str,
        default="name",
        help="Column name to read from CSV/TSV inputs (default: name)",
    )
    p.add_argument(
        "--output", "-o", type=str, help="Write results to a file (default: stdout)"
    )
    p.add_argument(
        "--output-format",
        type=str,
        default=None,
        choices=["json", "csv", "tsv"],
        help="Optional override for output format (otherwise inferred from output file extension)",
    )
    p.add_argument(
        "--smiles-selection-mode",
        default="weighted",
        type=str,
        help="Smiles selection mode",
    )
    p.add_argument(
        "--detailed-name-dict",
        default=False,
        type=bool,
        help="Whether to return a detailed name dictionary",
    )
    p.add_argument("--batch-size", default=500, type=int, help="Batch size")
    p.add_argument(
        "--normalize-unicode",
        default=True,
        type=bool,
        help="Whether to normalize unicode",
    )
    p.add_argument(
        "--split-names-to-solve",
        default=True,
        type=bool,
        help="Whether to split names to solve",
    )
    p.add_argument(
        "--resolve-peptide-shorthand",
        default=True,
        type=bool,
        help="Whether to resolve peptide shorthand",
    )
    p.add_argument(
        "--attempt-name-correction",
        default=True,
        type=bool,
        help="Whether to attempt name correction",
    )
    p.add_argument(
        "--internet-connection-available",
        default=True,
        type=bool,
        help="Whether an internet connection is available",
    )
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    names = list(args.names)
    if args.input:
        names.extend(
            read_names_from_file(
                args.input,
                input_format=args.input_format,
                input_column=args.input_column,
            )
        )

    if not names:
        parser.error("Provide names as arguments or via --input")

    results = resolve_compounds_to_smiles(
        names,
        smiles_selection_mode=args.smiles_selection_mode,
        detailed_name_dict=args.detailed_name_dict,
        batch_size=args.batch_size,
        normalize_unicode=args.normalize_unicode,
        split_names_to_solve=args.split_names_to_solve,
        resolve_peptide_shorthand=args.resolve_peptide_shorthand,
        attempt_name_correction=args.attempt_name_correction,
        internet_connection_available=args.internet_connection_available,
    )

    write_results(
        results,
        output_path=args.output,
        output_format=args.output_format,
    )

    return 0
